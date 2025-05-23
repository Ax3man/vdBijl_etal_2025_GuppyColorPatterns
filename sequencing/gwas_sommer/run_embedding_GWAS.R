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

# load embedding data, summarise to individual level and add relevant metadata
embeddings <- full_join(
  read_rds(
    'dimension_reduction/triplet_loss_encoders/embeddings/car_model_ped_fullcolor_comparison/embed_dim_5.rds'
  ) |> as.data.frame() |> setNames(paste0('car_V', 1:5)) |> rownames_to_column('unique_id'),
  read_rds(
    'dimension_reduction/triplet_loss_encoders/embeddings/mel_model_ped_fullcolor_comparison/embed_dim_5.rds'
  ) |> as.data.frame() |> setNames(paste0('mel_V', 1:5)) |> rownames_to_column('unique_id'),
  join_by(unique_id)
) |>
  left_join(pd, join_by(unique_id)) |>
  summarise(across(car_V1:mel_V5, mean), .by = c(fish_id, facing_direction)) |>
  summarise(across(car_V1:mel_V5, mean), .by = fish_id) |>
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
#samp <- sample(rownames(A), 100); plot(H[samp, samp], A[samp, samp]); abline(0, 1)

cat('Loading genotypes into RAM\n')

geno <- load_vcf(reference = reference, min_maf = 0.1)
#geno <- load_vcf(region = 'NC_024344.1:26587000-26588000', min_maf = 0.1) # for testing
rownames(geno$mat) <- ss$fish_id[match(rownames(geno$mat), ss$sample_name)]

# subset all matrices to only the males.
H2 <- H[embeddings$fish_id, embeddings$fish_id]
X2 <- X[embeddings$fish_id, embeddings$fish_id]
Y2 <- Y[embeddings$fish_id, embeddings$fish_id]
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
  cat(glue::glue('Analyzing {trait}\n\n'))
  df_geno <- filter(df, fish_id %in% geno_id)

  if (!file.exists(modelfile)) {
    cat('Obtaining estimates of variance components by REML\n')
    pre <- mmer(
      #nIters = 1, tolParConvLL = 1e4, # for rapid testing

      fixed = reformulate(c('selection2'), trait, intercept = FALSE),
      random = ~ vsr(fish_idA, Gu = H2) + vsr(fish_idX, Gu = X2) + vsr(fish_idY, Gu = Y2),
      rcov = ~ units,
      data = df
    )
    write_rds(pre, modelfile, compress = "xz")
  } else {
    cat('Loading estimates of variance components from disk\n')
    pre <- read_rds(modelfile)
  }
  # we have to assign to the global environment, because {sommer} uses eval(parse()) :( :( :(
  Sigma <<- pre$sigma_scaled

  random_formula <- ~ vsr(fish_idA, Gu = H3, Gti = Sigma[['u:fish_idA']], Gtc = matrix(3)) +
    vsr(fish_idX, Gu = X3, Gti = Sigma[['u:fish_idX']], Gtc = matrix(3)) +
    vsr(fish_idY, Gu = Y3, Gti = Sigma[['u:fish_idY']], Gtc = matrix(3))

  cat('Obtaining marker-wise coeffients and scores\n')
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

fitGWAS(trait = 'car_V4', reference = 'male', df = embeddings, geno = geno)
fitGWAS(trait = 'car_V5', reference = 'male', df = embeddings, geno = geno)

# DO NOT RUN IN PARALLEL - THE RANDOM EFFECT FORMULA WILL NOT WORK

# walk(
#   paste0('car_V', 1:5),
#   \(t) fitGWAS(trait = t, reference = reference, df = embeddings, geno = geno)
# )
# walk(
#   paste0('mel_V', 1:5),
#   \(t) fitGWAS(trait = t, reference = reference, df = embeddings, geno = geno)
# )

#options(future.globals.maxSize = 1e12)
#plan(multicore, workers = 15)
# future_walk(
#   c(paste0('car_', 1:7), paste0('mel_', 1:8)),
#   \(t) fitGWAS(trait = t, reference = reference, df = embeddings, geno = geno)#,
#   .options = furrr_options(seed = TRUE)
# )
# plan(sequential)
