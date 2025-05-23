suppressPackageStartupMessages({
  library(sommer)
  library(tidyverse)
  library(furrr)

  source('quant_gen/prepare_pedigrees.R')
  source('sequencing/genomics_helpers.R')
  source('selection_decisions/compile_decisions.R')
  source('sequencing/gwas/peak_viz/peak_viz_tools.R')
  source('sequencing/gwas_sommer/two_step_GWAS.R')
})

reference <- 'male'

cat('Loading phenotypes and pedigree\n')

# load photodatabase with fish ids etc
pd <- data.table::fread('photo_database.csv') |>
  dplyr::select(replicate, generation, fish_id, unique_id, facing_direction)

# load ornament data, summarise to individual level and add relevant metadata
ornaments <- bind_rows(
  read_rds('ornament_analysis/car_ornaments.rds'),
  read_rds('ornament_analysis/mel_ornaments.rds')
) |>
  pivot_wider(id_cols = unique_id, names_from = ornament, values_from = present_10) |>
  left_join(pd, join_by(unique_id)) |>
  summarise(across(c(car_1:car_7, mel_1:mel_8), mean), .by = c(fish_id, facing_direction)) |>
  summarise(across(c(car_1:car_7, mel_1:mel_8), mean), .by = fish_id) |>
  left_join(
    selection |> dplyr::select(fish_id, replicate, selection2) |> mutate(fish_id = tolower(fish_id)),
    join_by(fish_id)
  ) |> add_patriline(ped_df) |>
  mutate(fish_idA = fish_id, fish_idX = fish_id, fish_idY = fish_id)

# get the genomic relationship matrix
ss <- get_sampling_structure()
grm <- get_GRM(reference)
nm <- ss$fish_id[match(rownames(grm), ss$sample_name)]
dimnames(grm) <- list(nm, nm)

# adjust the scale of the GRM, so that it matches the scale of A
grm_rescaled <- grm * (mean(diag(A)) / mean(diag(grm)))

# Make H, the ssGBLUP combined matrix of A and the GRM
H <- H.mat(A, grm_rescaled)
samp <- sample(rownames(A), 100); plot(H[samp, samp], A[samp, samp]); abline(0, 1)

cat('Loading genotypes into RAM\n')

geno <- load_vcf(reference = reference, min_maf = 0.1)
#geno <- load_vcf(region = 'NC_024344.1:26587000-26588000', min_maf = 0.1) # for testing
rownames(geno$mat) <- ss$fish_id[match(rownames(geno$mat), ss$sample_name)]

# subset all matrices to only the males.
H2 <- H[ornaments$fish_id, ornaments$fish_id]
X2 <- X[ornaments$fish_id, ornaments$fish_id]
Y2 <- Y[ornaments$fish_id, ornaments$fish_id]
# and make a subset with just the genotyped males
geno_id <- rownames(geno$mat)
H3 <- H2[geno_id, geno_id] # note that this is just the GRM again
X3 <- X2[geno_id, geno_id]
Y3 <- Y2[geno_id, geno_id]

fitGWAS <- function(
    trait, reference, df, geno,
    outfile = glue::glue('sequencing/gwas_sommer/sommer_output/{reference}_{trait}.csv.gz'),
    modelfile = glue::glue('sequencing/gwas_sommer/sommer_models/{reference}_{trait}.rds')
  ) {
  df_geno <- filter(df, fish_id %in% geno_id)

  if (!file.exists(modelfile)) {
    pre <- mmer(
      #nIters = 1, tolParConvLL = 1,

      fixed = reformulate(c('selection2'), trait, intercept = FALSE),
      random = ~ vsr(fish_idA, Gu = H2) + vsr(fish_idX, Gu = X2) + vsr(fish_idY, Gu = Y2),
      rcov = ~ units,
      data = df
    )
    write_rds(pre, modelfile, compress = "xz")
  } else {
    pre <- read_rds(modelfile)
  }
  # we have to assign to the global environment, because {sommer} uses eval(parse()) :( :( :(
  Sigma <<- pre$sigma_scaled

  random_formula <- ~ vsr(fish_idA, Gu = H3, Gti = Sigma[['u:fish_idA']], Gtc = matrix(3)) +
    vsr(fish_idX, Gu = X3, Gti = Sigma[['u:fish_idX']], Gtc = matrix(3)) +
    vsr(fish_idY, Gu = Y3, Gti = Sigma[['u:fish_idY']], Gtc = matrix(3))

  m <- GWAS(
    fixed = reformulate(c('selection2'), trait, intercept = FALSE),
    random = random_formula,
    rcov = ~ units,
    data = df_geno,
    M = geno$mat, gTerm = 'u:fish_idA'
  )
  out <- geno$df
  out$p.value <- m$pvals[, 1]
  out$score <- m$scores[, 1]
  out$beta <- m$effects[, 1]
  out$beta.se <- m$effects.se[, 1]
  data.table::fwrite(out, outfile)
}

cat('Starting GWAS\n')

# DO NOT RUN IN PARALLEL - THE RANDOM EFFECT FORMULA WILL NOT WORK
walk(
  c(paste0('car_', 1:7), paste0('mel_', 1:8)),
  \(t) fitGWAS(trait = t, reference = reference, df = ornaments, geno = geno)
)
