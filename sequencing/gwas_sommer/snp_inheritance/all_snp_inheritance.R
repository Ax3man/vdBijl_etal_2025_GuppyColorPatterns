library(tidyverse, warn.conflicts = FALSE)

## Define helper functions
source('sequencing/genomics_helpers.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')

source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('sequencing/gwas/peak_viz/panel_functions/make_linkage_panel.R')

get_linkage <- function(geno, workers = 25, reference = 'female') {
  suppressPackageStartupMessages(require(furrr))

  workers <- pmin(parallelly::availableCores(logical = FALSE) - 1, workers)

  cat('Preparing data...\n')

  source('quant_gen/prepare_pedigrees.R')
  sampling <- get_sampling_structure() |> add_patriline(ped_df)
  A_small <- A[sampling$fish_id, sampling$fish_id]
  X_small <- X[sampling$fish_id, sampling$fish_id]

  sample_pats <- sampling$patriline[match(rownames(geno$mat), sampling$sample_name)]
  rownames(geno$mat) <- sampling$fish_id[match(rownames(geno$mat), sampling$sample_name)]

  if (workers > 1) plan(multisession, workers = workers)
  options(future.globals.maxSize = Inf)

  library(progressr)
  progressr::handlers(progressr::handler_progress(intrusiveness = 100.0))

  linkage_wrapper <- function(g) {
    df <- data.frame(dose = g, fish_id = names(g), fish_idX = names(g), patriline = sample_pats)
    l <- safe_snp_sex_linkage_AIC(df, A_small, X_small)
    case_when(
      nrow(l) == 0 ~ NA_character_,
      l$AIC_weight[2] > 0.8 ~ "X",
      l$AIC_weight[3] > 0.8 ~ "Y",
      TRUE ~ 'auto'
    )
  }

  with_progress({
    p <- progressor(ncol(geno$mat))
    geno$df$source <- future.apply::future_apply(
      geno$mat,
      2,
      \(g) { p('Fitting linkage models'); linkage_wrapper(g) },
      simplify = TRUE, future.seed = TRUE, future.chunk.size = 100
    )
  })

  if (workers > 1) plan(sequential)

  cat('Done!\n')
  return(geno$df)
}

find_snp_inheritance <- function(
  reference = 'female', chr, region, overwrite = FALSE
) {
  outfile <- glue::glue('sequencing/gwas_sommer/snp_inheritance/sommer_snp_inheritance/{reference}/{chr}.csv.gz')
  if (!overwrite && file.exists(outfile)) {
    message('File exists and overwrite is FALSE, so skipping')
    return(invisible(NULL))
  }

  geno <- load_vcf(reference, region = region, min_maf = 0.1)

  link <- get_linkage(geno, workers = 24)

  data.table::fwrite(link, outfile)
}

for (i in 1:nrow(scaff_labs)) {
  cat(str_glue('Fitting {scaff_labs$chr[i]}...\n\n'))
  find_snp_inheritance(
    reference = 'female',
    chr = scaff_labs$chr[i],
    region = str_glue_data(scaff_labs[i, ], '{chr}:1-{scaff_len}')
  )
}

for (i in 1:nrow(male_scaff_labs)) {
  cat(str_glue('Fitting {male_scaff_labs$chr[i]}...\n\n'))
  find_snp_inheritance(
    reference = 'male',
    chr = male_scaff_labs$chr[i],
    region = str_glue_data(male_scaff_labs[i, ], 'chr{chr}:1-{scaff_len}')
  )
}

