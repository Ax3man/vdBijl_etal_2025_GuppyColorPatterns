get_gwas_results <- function(
    trait, chr = NULL, pos, viz_range, fdr.level, pval_column = 'p_SHet', reference = 'female'
) {
  if (reference == 'female') {
    gwas <- data.table::fread(glue::glue('sequencing/gwas/gemma_output/{trait}.assoc.txt'))
  } else {
    gwas <- data.table::fread(glue::glue('sequencing/male_reference_gwas/gemma_output/{trait}.assoc.txt'))
  }

  if (!is.null(chr)) {
    gwas <- dplyr::filter(
      gwas,
      .data$chr == .env$chr,
      between(ps, viz_range[1], viz_range[2])
    )
  }
  return(gwas)
}

get_genotypes <- function(vcf, chr, range) {
  message('Consider using sequencing/genomics_helpers/load_vcf instead.')

  if (chr %in% as.character(1:23)) chr <- paste0('chr', chr)
  #require(VariantAnnotation)
  if (length(range) == 1) {
    regions <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = range, width = 1)
    )
  } else {
    regions <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = range[1], width = diff(range))
    )
  }
  VariantAnnotation::readVcfAsVRanges(
    vcf,
    param = VariantAnnotation::ScanVcfParam(which = regions)
  ) %>%
    as.data.frame() %>%
    filter(
      !(sampleNames %in% c(
        'NS.2125.002.IDT_i7_111---IDT_i5_111.280',
        'NS.2145.001.IDT_i7_89---IDT_i5_89.355'
      ))
    )
}

get_phenotypes <- function(trait) {
  require(dplyr)
  if (trait %in% c('car_perc', 'mel_perc_v2')) {

  } else {
    data.table::fread(glue::glue('sequencing/kmer_gwas/config/phenotypes/{trait}.pheno')) %>%
      inner_join(
        get_sampling_structure() %>% dplyr::select(-sample_id),
        join_by(accession_id == fish_id)
      ) %>% dplyr::rename(fish_id = accession_id)
  }
}

get_ornaments <- function() {
  bind_rows(
    read_rds('ornament_analysis/car_ornaments.rds'),
    read_rds('ornament_analysis/mel_ornaments.rds')
  ) %>%
    pivot_wider(id_cols = unique_id, names_from = 'ornament', values_from = 'present_10') %>%
    left_join(
      data.table::fread('photo_database.csv') %>% dplyr::select(fish_id, unique_id, facing_direction),
      join_by(unique_id)
    ) %>%
    group_by(fish_id, facing_direction) %>%
    summarise(across(c(car_1:car_7, mel_1:mel_8), mean), .groups = 'drop_last') %>%
    summarise(across(c(car_1:car_7, mel_1:mel_8), mean), .groups = 'drop') %>%
    mutate(across(c(car_1:car_7, mel_1:mel_8), round)) %>%
    right_join(get_sampling_structure(), join_by(fish_id))
}

get_phenotypes2 <- function(type = 'present_10') {
  out <- bind_rows(
    read_rds('ornament_analysis/car_ornaments.rds'),
    read_rds('ornament_analysis/mel_ornaments.rds')
  ) %>%
    pivot_wider(id_cols = unique_id, names_from = 'ornament', values_from = all_of(type)) %>%
    left_join(
      data.table::fread('photo_database.csv') %>%
        dplyr::select(fish_id, unique_id, facing_direction, car_perc, mel_perc_v2),
      join_by(unique_id)
    ) %>%
    group_by(fish_id, facing_direction) %>%
    summarise(across(c(car_1:car_7, mel_1:mel_8, car_perc, mel_perc_v2), mean), .groups = 'drop_last') %>%
    summarise(across(c(car_1:car_7, mel_1:mel_8, car_perc, mel_perc_v2), mean), .groups = 'drop')

  if (type %in% c('present', 'present_10')) {
    out <- mutate(out, across(c(car_1:car_7, mel_1:mel_8), round))
  }

  right_join(out, get_sampling_structure(), join_by(fish_id))
}

get_sampling_structure <- function() {
  require(dplyr)
  sampling1 <- data.table::fread('sequencing/sample_lists/DNA_WGS_300samples_long_list.csv')
  sampling2 <- data.table::fread('data/NovaSeq_sample_list.csv') %>%
    dplyr::select(sample_id = Name, sample_name = `Filename Prefix`)
  inner_join(sampling1, sampling2, by = 'sample_id') %>%
    dplyr::select(fish_id, sample_id, sample_name) %>%
    mutate(
      fish_id = tolower(fish_id),
      # adjust these names, for samples that were sequenced twice.
      sample_name = case_when(
        sample_name == 'NS.2145.001.IDT_i7_77---IDT_i5_77.289' ~ 'sample_289',
        sample_name == 'NS.2145.001.IDT_i7_144---IDT_i5_144.323' ~ 'sample_323',
        TRUE ~ sample_name
      )
    ) %>%
    filter(!(sample_name %in% c(
      # remove duplicated entries due to double sequencing
      'NS.2118.002.IDT_i7_155---IDT_i5_155.323', 'NS.2118.002.IDT_i7_19---IDT_i5_19.289',
      # remove dropped sample
      'NS.2125.002.IDT_i7_111---IDT_i5_111.280', 'NS.2145.001.IDT_i7_89---IDT_i5_89.355'
    )))
}

get_mean_coverage <- function(reference = 'female') {
  if (reference == 'female') {
    cov <- read.csv('sequencing/gwas/coverage_summary.csv')
  } else {
    cov <- read.csv('sequencing/male_reference_gwas/coverage_summary.csv')
  }
  return(cov)
}

get_snp_effect_heatmap <- function(trait, chr, pos) {
  files <- list.files(glue::glue('sequencing/gwas/snp_effects/by_pixel/{trait}'), f = TRUE)
  files <- files[str_detect(files, chr)]
  file <- files[str_detect(files, paste0('-', pos, '-'))]
  if (length(file) != 1) {
    warning('Multiple snp-effect heatmaps found for this position.')
    return(ggplot())
  }
  img <- magick::image_read(file)
  if (magick::image_info(img)$width < 10) return(ggplot())
  cowplot::ggdraw() + cowplot::draw_image(img)
}

get_region_coverage <- function(chr, pos, viz_range, reference = 'female') {
  require(tidyverse)
  if (reference == 'female') {
    f <- list.files('sequencing/coverage_gwas_intervals', f = TRUE) %>% str_subset(chr)
  } else {
    f <- list.files('sequencing/male_reference_gwas/coverage_gwas_intervals', f = TRUE) %>%
      str_subset(paste0('chr', chr))
  }
  d <- map_dfr(f, data.table::fread)

  if (nrow(d) == 0) return(data.frame())

  x <- inner_join(
    data.frame(viz_min = viz_range[1], viz_max = viz_range[2]),
    d,
    join_by(overlaps(viz_min, viz_max, window_start, window_end))
  ) %>%
    dplyr::select(-viz_min, -viz_max)

  if (reference == 'male') {
    # fix 2 sample names
    x <- dplyr::rename(x, sample_323 = `323_merged`, sample_289 = `289_merged`)
  }

  x %>%
    pivot_longer(-(chr:window_end), names_to = 'sample_name', values_to = 'coverage') %>%
    left_join(get_mean_coverage(), join_by(sample_name)) %>%
    mutate(rel_coverage = coverage / mean_coverage)
}

# yuying parents
yuying_parents <- data.frame(
  fish_id = c('P433', 'P481', 'P409', 'P363', 'P271', 'FR26'),
  sex = c('female', 'female', 'male', 'male', 'male', 'female'),
  sample_name = c(
    'NS.1357.002.IDT_i7_9---IDT_i5_9.P433', 'NS.1357.002.IDT_i7_82---IDT_i5_82.P481',
    'NS.1357.002.IDT_i7_43---IDT_i5_43.P409', 'NS.1357.002.IDT_i7_55---IDT_i5_55.P363',
    'NS.1357.002.IDT_i7_31---IDT_i5_31.P271', 'NS.1357.002.IDT_i7_67---IDT_i5_67.FR26'
  )
)
get_mean_coverage_yuying <- function() {
  list.files('sequencing/yuying_parents/coverage/', f = TRUE) %>%
    setNames(., basename(.) %>% tools::file_path_sans_ext() %>% str_remove('_average_coverage')) %>%
    map_dfr(data.table::fread, .id = 'sample_name') %>%
    dplyr::rename(mean_coverage = V1)
}
get_region_coverage_yuying <- function(chr, pos, viz_range) {
  f <- list.files('sequencing/yuying_parents/coverage_gwas_intervals', f = TRUE) %>%
    str_subset(chr)
  d <- map_dfr(f, data.table::fread)

  if (nrow(d) == 0) return(data.frame())

  inner_join(
    data.frame(viz_min = viz_range[1], viz_max = viz_range[2]),
    d,
    join_by(overlaps(viz_min, viz_max, window_start, window_end))
  ) %>%
    dplyr::select(-viz_min, -viz_max) %>%
    pivot_longer(-(chr:window_end), names_to = 'sample_name', values_to = 'coverage') %>%
    left_join(get_mean_coverage_yuying(), join_by(sample_name)) %>%
    mutate(rel_coverage = coverage / mean_coverage) %>%
    left_join(yuying_parents, join_by(sample_name))
}
