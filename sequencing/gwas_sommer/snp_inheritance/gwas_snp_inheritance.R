suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(tidyverse)
})

## Define helper functions
source('sequencing/genomics_helpers.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')

source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('sequencing/gwas/peak_viz/panel_functions/make_linkage_panel.R')

get_linkage <- function(geno, trait, workers = 25, min_MAF = 0.1, reference = 'female', use_previous = TRUE) {
  suppressPackageStartupMessages(require(furrr))

  workers <- pmin(parallelly::availableCores(logical = FALSE) - 1, workers)

  cat('Preparing data...\n')

  source('quant_gen/prepare_pedigrees.R')
  sampling <- get_sampling_structure() |> add_patriline(ped_df)
  A_small <- A[sampling$fish_id, sampling$fish_id]
  X_small <- X[sampling$fish_id, sampling$fish_id]

  # load previous results (only significant variants)
  if (use_previous) {
    prev_file <- glue::glue(
      'sequencing/gwas_sommer/snp_inheritance/gwas_snp_inheritance/{reference}_{trait}.csv.gz'
    )

    if (file.exists(prev_file)) {
      prev <- data.table::fread(prev_file)
    } else {
      prev <- data.frame(chr = character(), start = integer(), end = integer(), alt = character())
    }
  }

  dat <- geno %>%
    as.data.frame() %>%
    mutate(MAF = ifelse(AF1 < 0.5, AF1, 1 - AF1)) %>%
    filter(MAF >= min_MAF) %>%
    dplyr::select(chr = seqnames, start, end, alt, sampleNames, GT) %>%
    inner_join(sampling, join_by(sampleNames == sample_name)) %>%
    mutate(
      dose = case_when(GT == '0/0' ~ 0, GT == '0/1' ~ 1, GT == '1/0' ~ 1, GT == '1/1' ~ 2, TRUE ~ NA_integer_),
      fish_idX = fish_id
    )

  if (use_previous) {
    done <- semi_join(prev, dat, join_by(chr, start, end, alt))
    to_do <- anti_join(dat, prev, join_by(chr, start, end, alt))
  } else {
    to_do <- dat
  }

  cat('Fitting models...\n')
  if (workers > 1) plan(multisession, workers = workers)
  link <- to_do %>%
    group_by(chr, start, end, alt) %>%
    group_nest()

  link$out <- future_map(
    link$data,
    \(x) safe_snp_sex_linkage_AIC(x, A_small, X_small),
    .options = furrr_options(
      seed = TRUE,
      globals = c('A_small', 'X_small', 'safe_snp_sex_linkage_AIC', 'snp_sex_linkage_AIC')
    )
  )
  if (workers > 1) plan(sequential)

  cat('Post-processing...\n')
  link <- link %>%
    unnest(cols = out) %>%
    dplyr::select(-data)

  if (use_previous) {
    link <- bind_rows(link, done)
  }
  cat('Done!\n')
  return(link)
}

find_snp_inheritance <- function(
    traits, name = traits, pval_column, reference = 'female', overwrite = TRUE
  ) {
  outfile <- glue::glue('sequencing/gwas_sommer/snp_inheritance/sommer_snp_inheritance/{reference}_{name}.csv.gz')
  if (!overwrite && file.exists(outfile)) {
    message('File exists and overwrite is FALSE, so skipping')
    return(invisible(NULL))
  }

  cat('Loading GWAS results...\n')
  if (length(traits) == 1) {
    g <- prep_sommer_gwas_table(traits, reference, pval_column, FALSE, FALSE) |>
      mutate(
        p.value.adj = p.adjust(.data[[pval_column]], method = 'BH'),
        significant = p.value.adj < 0.05
      ) |>
      filter(significant)
  } else {
    g <- map_dfr(
      setNames(traits, traits),
      \(tr) prep_sommer_gwas_table(tr, reference, pval_column, FALSE, FALSE),
      .id = 'trait'
    ) |>
    # the joint-set q-value
      mutate(
        p.value.adj = p.adjust(.data[[pval_column]], method = 'BH'),
        significant = p.value.adj < 0.05
      ) |>
      filter(significant)
  }

  cat('Loading genotypes from VCF...\n')
  if (reference == 'female') {
    vcf <- 'sequencing/gwas/filtered.vcf.gz'
  } else {
    stop('male vcf still needs to be specified')
  }

  geno <- load_vcf(reference, region = glue::glue_data(g, '{chr}:{pos}-{pos}'))

  cat('Estimating linkage...\n')
  link <- get_linkage(geno, trait = traits, workers = 24, min_MAF = 0.1, use_previous = FALSE)

  data.table::fwrite(link, outfile)
}

find_snp_inheritance('car_PIE', pval_column = 'p_SHet', reference = 'female')
find_snp_inheritance('mel_PIE', pval_column = 'p_SHet', reference = 'female')
find_snp_inheritance(
  glue::glue('car_{c(1:3, 5:7)}'), name = 'orange_ornaments',
  pval_column = 'p.value', reference = 'female'
)
find_snp_inheritance(
  glue::glue('mel_{c(1, 3, 5:8)}'), name = 'black_ornaments',
  pval_column = 'p.value', reference = 'female'
)

