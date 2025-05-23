library(imager)
library(tidyverse)
library(sommer)
source('sequencing/gwas_sommer/gwas_sommer_tools.R')
source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('sequencing/genomics_helpers.R')
source('quant_gen/prepare_pedigrees.R')
source('selection_decisions/compile_decisions.R')

reference <- 'female'
trait <- 'orange_ornaments'
pval_column <- 'p.value'
fdr.level <- 0.05
min_dist <- 1e5 / 2

# Load GWAS results --------------------------------------------------------------------------------
reduce_peaks <- function(d, min_dist = 1e4, pval_column = 'p_lrt') {
  # sort SNPs by significance, and assign windows around them
  d <- arrange(d, !!sym(pval_column)) %>%
    mutate(.start = pos - min_dist, .end = pos + min_dist)

  out <- dplyr::slice(d, 1)
  while (TRUE) {
    a <- anti_join(dplyr::slice(d, -1), out, join_by(chr, between(x$pos, y$.start, y$.end)))
    if (nrow(a) == nrow(d)) break
    d <- a
    out <- bind_rows(out, dplyr::slice(d, 1))
  }
  out <- dplyr::select(out, -.start, -.end)
  return(out)
}

gwas <- prep_sommer_gwas_table(
  str_glue('car_{c(1:3, 5:7)}'), name = 'orange_ornaments', pval_column = pval_column
) |>
  filter(
    str_starts(chr, 'NC'),
    p.value.adj < fdr.level
  )

gwas_selected <- gwas |>
  filter(!(source %in% c('X', 'Y'))) |>
  group_by(trait) |>
  group_modify(\(dat, key) reduce_peaks(dat, min_dist = min_dist, pval_column = pval_column)) |>
  ungroup() |>
  bind_rows(
    gwas |> filter(source %in% c('X', 'Y')) |> arrange(ref) |>
      group_by(chr == 'NC_024342.1', source, trait) |>
      slice_min(.data[[pval_column]], with_ties = FALSE)
  )

# Load genotypes -----------------------------------------------------------------------------------
ss <- get_sampling_structure()
geno <- load_vcf(region = str_glue_data(gwas_selected, '{chr}:{pos}-{pos}'))
rownames(geno$mat) <- ss$fish_id[match(rownames(geno$mat), ss$sample_name)]

# Load the VCV matrices ----------------------------------------------------------------------------
# get the genomic relationship matrix
grm <- get_GRM(reference)
nm <- ss$fish_id[match(rownames(grm), ss$sample_name)]
dimnames(grm) <- list(nm, nm)

# adjust the scale of the GRM, so that it matches the scale of A
grm_rescaled <- grm * (mean(diag(A)) / mean(diag(grm)))

# Make H, the ssGBLUP combined matrix of A and the GRM
H <- H.mat(A, grm_rescaled)

# Load the phenotypes ------------------------------------------------------------------------------
img_files <- list.files('data/carotenoid_coloration_warped', recursive = TRUE, full.names = TRUE) |>
  (\(x) setNames(x, basename(x) |> tools::file_path_sans_ext()))()

# first, match the sequenced individuals with all their color pattern images
phenotypes <- data.table::fread('photo_database.csv') |>
  dplyr::select(fish_id, facing_direction, unique_id) |>
  inner_join(ss, join_by(fish_id)) |>
  # load those images
  filter(file.exists(img_files[unique_id])) |>
  mutate(image = map(
    img_files[unique_id],
    \(p) load.image(p) |> channel(4) |> resize_halfXY() |> threshold(thr = 0.5)
  )) |>
  # reduce to 1 measurment per individual, by averaging the images first per side, then per fish.
  group_by(fish_id, facing_direction) |>
  summarise(image = list(average(image)), .groups = 'drop_last') |>
  summarise(image = list(average(image)), .groups = 'drop') |>
  mutate(image = map(image, as.data.frame)) |>
  unnest(cols = image) |>
  # since we have averaged a lot of binary values, we need to round again to binary values
  # mutate(value = round(value)) %>%
  # finally, drop pixels where < 1% of sequenced males have color
  group_by(x, y) |>
  filter(mean(value) > 0.01) |>
  left_join(
    dplyr::select(selection, fish_id, selection2) |> mutate(fish_id = tolower(fish_id)),
    join_by(fish_id)
  )

# Filter the VCV matrices --------------------------------------------------------------------------
# subset all matrices to only the males.
f <- unique(phenotypes$fish_id)
H2 <- H[f, f]
X2 <- X[f, f]
Y2 <- Y[f, f]
# and make a subset with just the genotyped males
geno_id <- rownames(geno$mat)
H3 <- H2[geno_id, geno_id] # note that this is just the GRM again
X3 <- X2[geno_id, geno_id]
Y3 <- Y2[geno_id, geno_id]


cat('Analyzing', n_groups(phenotypes), 'pixels for', ncol(geno$mat), 'SNPs...\n')

# Run pixelwise SNP effect estimation --------------------------------------------------------------

#trait <- 'value'; reference = 'female'; df = group_split(phenotypes)[[1]];
fitGWAS <- function(
    trait, reference, df, geno,
    name = str_glue('pixel_{str_pad(df$x[1], 4, pad = "0")}_{str_pad(df$y[1], 4, pad = "0")}'),
    outfile = str_glue('sequencing/gwas_sommer/snp_effects/pixel_results/orange_ornaments/{reference}_{name}.csv.gz')
) {
  df_geno <- filter(df, fish_id %in% geno_id)

  brms_file <- str_glue('quant_gen/cluster_brms_carotenoid/posterior_summaries_{1:2}/{name}.rds')
  if (!any(file.exists(brms_file))) return(NULL)
  brms_post_summ <- read_rds(brms_file[file.exists(brms_file)]) |>
    filter(variable %in% str_glue('sd_{c(1:3, 5)}[1]')) |>
    pull(mean) |> setNames(c('res', 'A', 'X', 'Y'))
  # first rescale the variances to unity
  brms_post_summ <- brms_post_summ / sum(brms_post_summ)
  # then scale to sum to the variance of y. Note that we rescale including the residual variance,
  # which is the individual level in the brms model, but this is still estimated in the GWAS fit.
  brms_post_summ <- brms_post_summ / var(df_geno$value)

  # random_formula <- ~
  #   vsr(fish_idA, Gu = H3, Gti = matrix(brms_post_summ['A']), Gtc = matrix(3)) +
  #   vsr(fish_idX, Gu = X3, Gti = matrix(brms_post_summ['X']), Gtc = matrix(3)) +
  #   vsr(fish_idY, Gu = Y3, Gti = matrix(brms_post_summ['Y']), Gtc = matrix(3))

  # we are forced to do global shenanigans, since there is some internal eval(parse) in sommer
  assign(name, brms_post_summ, envir = .GlobalEnv)
  random_formula <- reformulate(c(
    str_glue("vsr(fish_idA, Gu = H3, Gti = matrix({name}['A']), Gtc = matrix(3))"),
    str_glue("vsr(fish_idX, Gu = X3, Gti = matrix({name}['X']), Gtc = matrix(3))"),
    str_glue("vsr(fish_idY, Gu = Y3, Gti = matrix({name}['Y']), Gtc = matrix(3))")
  ))

  m <- GWAS(
    fixed = reformulate(c('selection2'), trait, intercept = FALSE),
    random = random_formula,
    rcov = ~ units,
    data = df_geno |> mutate(fish_idA = fish_id, fish_idX = fish_id, fish_idY = fish_id),
    M = geno$mat, gTerm = 'u:fish_idA'
  )
  out <- geno$df
  out$p.value <- m$pvals[, 1]
  out$score <- m$scores[, 1]
  out$beta <- m$effects[, 1]
  out$beta.se <- m$effects.se[, 1]
  data.table::fwrite(out, outfile)
}

l <- split(phenotypes, interaction(phenotypes$x, phenotypes$y))
for (i in 1:length(l)) fitGWAS('value', reference, l[[i]], geno)
# walk(
#   split(phenotypes, interaction(phenotypes$x, phenotypes$y)),
#   \(x) fitGWAS('value', reference, x, geno),
# )

# library(furrr)
# plan(multisession, workers = 24)
# future_walk(
#   split(phenotypes, interaction(phenotypes$x, phenotypes$y)),
#   \(x) fitGWAS('value', reference, x, geno),
#   .options = furrr_options(seed = TRUE), .progress = TRUE
# )
# plan(sequential)
