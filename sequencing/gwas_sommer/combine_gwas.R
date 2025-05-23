combine_gwas_sHet <- function(name, input, workers = 42, reference = 'female', overwrite = FALSE) {
  require(tidyverse)
  require(future)
  source('sequencing/gwas/CPASSOC/FunctionSet.R')

  outfile <- glue::glue('sequencing/gwas_sommer/sommer_output/{reference}_{name}.csv.gz')
  if (!overwrite && file.exists(outfile)) return(invisible(NULL))

  # Find the correct files
  files <- glue::glue('sequencing/gwas_sommer/sommer_output/{reference}_{input}.csv.gz')
  #if (!any(file.exists(files))) files <- str_remove(files, '\\.gz')
  assoc <- combine_assoc(files)

  # Obtain the genetic correlations
  cor_mat <- assoc %>%
    # filter out significant SNPs as recommended by CPASSOC
    filter(if_all(starts_with('Z'), \(x) abs(x) < 1.96)) %>%
    # remove MT, sex chromosomes and unplaced scaffolds as LD may be very high
    filter(!chr %in% c('NC_024238.1', 'NC_024342.1'), str_starts(chr, 'NW', negate = TRUE)) %>%
    group_by(chr) %>%
    # Reduce LD by taking every 100th SNP
    dplyr::slice(seq(1, n(), 100)) %>%
    ungroup() %>%
    dplyr::select(starts_with('Z')) %>%
    as.matrix() %>%
    cor()

  print(cor_mat)

  if (length(files) != nrow(cor_mat)) {
    stop('Files and correlations have mismatched dimensions.')
  }

  Z <- dplyr::select(assoc, starts_with('Z'))
  message('NOTE: Assuming a sample size of 297 individuals.')
  S <- rep(297, ncol(Z))
  B <- cor_mat

  cat('Estimating gamma...\n')
  para <- EstimateGamma(N = 1e6, SampleSize = S, CorrMatrix = B, correct = 1, isAllpossible = TRUE)

  cat('Calculating SHet statistics...\n')
  old_plan <- plan(multisession, workers = workers)
  options(future.globals.maxSize = 1e12)
  x <- SHet(X = Z, SampleSize = S, CorrMatrix = B, correct = 1, isAllpossible = TRUE)
  plan(old_plan)

  cat('Calculating p-values...\n')
  assoc$p_SHet <- pgamma(q = x - para[3], shape = para[1], scale = para[2], lower.tail = F)
  assoc$log_p_SHet <- pgamma(q = x - para[3], shape = para[1], scale = para[2], lower.tail = F, log.p = TRUE)

  data.table::fwrite(assoc, outfile)

  return(invisible(NULL))
}

#combine_gwas_sHet('car_PIE', paste0('car_V', 1:5), reference = 'female', overwrite = TRUE)
#combine_gwas_sHet('mel_PIE', paste0('mel_V', 1:5), reference = 'female', overwrite = TRUE)
combine_gwas_sHet('car_PIE', paste0('car_V', 1:5), reference = 'male', overwrite = TRUE)
combine_gwas_sHet('mel_PIE', paste0('mel_V', 1:5), reference = 'male', overwrite = TRUE)
