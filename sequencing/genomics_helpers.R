load_vcf <- function(
    reference = 'female', region = '', min_maf = NULL, format = 'GT', vartype = 'snps', pass = TRUE,
    drop_invariable = FALSE
) {
  if (reference == 'female') {
    #message("Loading extra filtered vcf...")
    vcf <- 'sequencing/gwas_sommer/filtered2.vcf.gz'
  } else {
    #message("Loading extra filtered vcf...")
    vcf <- 'sequencing/gwas_sommer/male_filtered2.vcf.gz'
  }

  load_vcf_region <- function(r, min_maf) {
    g <- vcfppR::vcftable(
      vcf,
      region = r,
      samples = '^NS.2125.002.IDT_i7_111---IDT_i5_111.280,NS.2145.001.IDT_i7_89---IDT_i5_89.355',
      vartype = vartype,
      format = format,
      pass = pass
    )
    if (format == 'GT') {
      names(g)[names(g) == 'gt'] <- 'GT'
      g[[format]] <- g[[format]] - 1
    } else {
      if (!is.null(min_maf) || drop_invariable) {
        gt <- vcfppR::vcftable(
          vcf,
          region = r,
          samples = '^NS.2125.002.IDT_i7_111---IDT_i5_111.280,NS.2145.001.IDT_i7_89---IDT_i5_89.355',
          vartype = vartype,
          format = 'GT',
          pass = pass
        )
        AF <- Rfast::rowmeans(gt$gt) / 2
        Vars <- Rfast::rowVars(gt$gt)
        if (!is.null(min_maf)) excl <- which(AF < min_maf | AF > (1 - min_maf))
        if (drop_invariable) excl <- which(Vars == 0)
        if (!is.null(min_maf) && drop_invariable) excl <- which(AF < min_maf | AF > (1 - min_maf) | Vars == 0)
        g[-1] <- lapply(g[-1], \(x) {
          if (length(dim(x)) == 0) return(x[-excl])
          if (length(dim(x)) == 2) return(x[-excl, ])
        })
      }
    }

    if (!is.matrix(g[[format]])) g[[format]] <- matrix(g[[format]], nrow = 1)

    if ((!is.null(min_maf) || drop_invariable) && format == 'GT') {
      AF <- Rfast::rowmeans(g$GT + 1) / 2
      Vars <- Rfast::rowVars(g$GT)
      if (!is.null(min_maf)) excl <- which(AF < min_maf | AF > (1 - min_maf))
      if (drop_invariable) excl <- which(Vars == 0)
      if (!is.null(min_maf) && drop_invariable) excl <- which(AF < min_maf | AF > (1 - min_maf) | Vars == 0)
      df <- as.data.frame(g[c('chr', 'pos', 'ref', 'alt')])[-excl, ]
      mat <- g[[format]][-excl, ] |> magrittr::set_colnames(g$samples) |> t()
    } else {
      df <- as.data.frame(g[c('chr', 'pos', 'ref', 'alt')])
      mat <- g[[format]] |> magrittr::set_colnames(g$samples) |> t()
    }
    return(list(df = df, mat = mat))
  }

  if (length(region) == 1) return(load_vcf_region(region, min_maf = min_maf))

  # support supplying a vector of regions
  l <- lapply(region, load_vcf_region, min_maf = min_maf) |> transpose()
    return(
    list(
      df = list_rbind(l$df),
      mat = do.call(cbind, l$mat)
    )
  )
}

get_GRM <- function(reference = 'female') {
  vcf <- 'sequencing/gwas/filtered.vcf.gz'

  # read a random region of the VCF to get the sample names encoded in the VCF
  sample_names <- VariantAnnotation::readVcf(
    vcf,
    param = VariantAnnotation::ScanVcfParam(
      which = GenomicRanges::GRanges("NC_024331.1", IRanges::IRanges(24115677))
    )
  ) %>%
    SummarizedExperiment::colData() %>% rownames()

  grm_file <- 'sequencing/gwas/gemma_output/for_kinship_comparison.cXX.txt'

  GRM <- data.table::fread(grm_file) %>%
    as.data.frame() %>%
    `dimnames<-`(list(sample_names, sample_names)) %>%
    as.matrix()

  to_drop <- c('NS.2125.002.IDT_i7_111---IDT_i5_111.280', 'NS.2145.001.IDT_i7_89---IDT_i5_89.355')
  GRM[!(sample_names %in% to_drop), !(sample_names %in% to_drop)]
}

prep_gwas_table <- function(
    trait, produce_plots = FALSE, pval_column, fdr_level = 0.05, calc_q = TRUE, reference = 'female',
    # whether to apply the Fraser remapping of LG12 coordinates, don't use if e.g. matching to vcf
    adjust_chr12 = TRUE
) {
  require(qvalue)
  if (reference == 'female') {
    gwas <- data.table::fread(glue::glue('sequencing/gwas/gemma_output/{trait}.assoc.txt'), data.table = FALSE)
    gwas$chr2 <- fct_relabel(
      factor(gwas$chr),
      \(l) case_when(l == 'NC_024238.1' ~ 'MT',
                     str_starts(l, 'NC') ~ (str_sub(l, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                     TRUE ~ 'Un'
      ))
    scaff_sizes <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
      dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
      dplyr::select(-V4, -V5)
  } else {
    if (reference != 'male') stop('`reference` should be "male" or "female".')
    source('sequencing/male_reference_gwas/update_chr12_liftover.R')

    gwas <- data.table::fread(glue::glue('sequencing/male_reference_gwas/gemma_output/{trait}.assoc.txt'), data.table = FALSE)
    gwas$chr2 <- fct_relabel(factor(gwas$chr), \(l) ifelse(str_starts(l, '0'), 'Un', l))
    # Note that something in the pipeline (plink?) has assumed that chr23 is the X (assumed human).
    gwas$chr2 <- factor(gwas$chr2, levels = c(1:22, 'X', 'Un'), labels = c(1:23, 'Un'))
    gwas$chr <- ifelse(gwas$chr == 'X', '23', gwas$chr)

    gwas <- mutate(gwas, chr = fct_relevel(factor(chr), as.character(1:23)))

    # for the male reference, we update the chr12 coordinates, following Fraser et al. and relevel the chr
    if (adjust_chr12) {
      gwas <- mutate(gwas, ps = ifelse(chr == '12', update_chr12(ps), ps))
    }

    scaff_sizes <- data.table::fread('sequencing/male_reference_gwas/male_reference.fna.fai') %>%
      dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
      dplyr::select(-V4, -V5) |>
      mutate(
        chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr),
        chr = fct_relevel(factor(chr), as.character(1:23))
      ) |>
      arrange(chr) |>
      mutate(cum_scaff_start = cumsum(lag(scaff_len, default = 0)))
  }

  # make qqplot
  if (produce_plots) {
    opa <- par(mfrow = c(2, 1))
    hist(gwas[[pval_column]])
    GWASTools::qqPlot(gwas[[pval_column]], thinThreshold = 4)
    par(opa)
  }

  if (calc_q) {
    q <- qvalue(gwas[[pval_column]], fdr.level = fdr_level)
    summary(q)
    gwas$qvalue <- q$qvalues
    gwas$significant <- q$significant
  }

  gwas %>%
    left_join(scaff_sizes, join_by(chr)) %>%
    arrange(chr2) |>
    mutate(cum_ps = ps + cum_scaff_start)
}

read_kmer_blast <- function(trait, parse_pvals = FALSE) {
  file <- glue::glue('sequencing/kmer_gwas/kmer_blast2ref/{trait}.tsv')
  d <- data.table::fread(file) %>% as_tibble()
  colnames(d) <- c(
    'kmer', 'chr', 'perc_ident', 'align_length', 'nr_mismatch', 'nr_gaps', 'kmer_start', 'kmer_end',
    'start', 'end', 'e_value', 'bit_score'
  )
  if (parse_pvals) {
    d <- d %>%
      separate(kmer, into = c('kmer', 'p_value'), sep = '_') %>%
      mutate(p_value = parse_number(p_value))
  }

  scaff_sizes <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
    dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
    dplyr::select(-V4, -V5) %>%
    mutate(chr2 = case_when(chr == 'NC_024238.1' ~ 'MT',
                            str_starts(chr, 'NC') ~ (str_sub(chr, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                            TRUE ~ 'Un'))

  d %>%
    left_join(scaff_sizes, join_by(chr)) %>%
    mutate(cum_start = start + cum_scaff_start)
}

scaff_labs <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
  dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
  mutate(
    scaff_mid = cum_scaff_start + 0.5 * scaff_len,
    chr2 = case_when(chr == 'NC_024238.1' ~ 'MT',
                     str_starts(chr, 'NC') ~ (str_sub(chr, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                     TRUE ~ 'Un')
  ) %>%
  filter(!(chr2 %in% c('Un', 'MT'))) %>%
  dplyr::select(-V4, -V5)

male_scaff_labs <- data.table::fread('sequencing/male_reference_gwas/male_reference.fna.fai') %>%
  dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
  dplyr::select(-V4, -V5) |>
  mutate(
    chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr),
    chr = fct_relevel(factor(chr), as.character(1:23))
  ) |>
  arrange(chr) |>
  mutate(
    cum_scaff_start = cumsum(lag(scaff_len, default = 0)),
    scaff_mid = cum_scaff_start + 0.5 * scaff_len,
    chr2 = chr
  ) |>
  filter(!str_starts(chr, '0'))

reduce_peaks <- function(d, min_dist = 1e4, pval_column = 'p_lrt') {
  # sort SNPs by significance, and assign windows around them
  d <- arrange(d, !!sym(pval_column)) %>%
    mutate(.start = ps - min_dist, .end = ps + min_dist)

  out <- dplyr::slice(d, 1)
  while (TRUE) {
    a <- anti_join(dplyr::slice(d, -1), out, join_by(chr, between(x$ps, y$.start, y$.end)))
    if (nrow(a) == nrow(d)) break
    d <- a
    out <- bind_rows(out, dplyr::slice(d, 1))
  }
  out <- dplyr::select(out, -.start, -.end)
  return(out)
}
