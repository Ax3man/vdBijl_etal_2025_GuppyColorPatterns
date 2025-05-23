prep_sommer_gwas_table <- function(
    trait, name = trait, reference = 'female', pval_column = 'p.value', qqplots = FALSE,
    # whether to apply the Fraser remapping of LG12 coordinates, don't use if e.g. matching to vcf
    adjust_chr12 = TRUE
) {
  if (length(trait) == 1) {
    gwas <- data.table::fread(
      glue::glue('sequencing/gwas_sommer/sommer_output/{reference}_{trait}.csv.gz'),
      data.table = FALSE
    )
  } else {
    gwas <- map(
      setNames(glue::glue('sequencing/gwas_sommer/sommer_output/{reference}_{trait}.csv.gz'), trait),
      \(x) data.table::fread(x, data.table = FALSE)
    ) |> list_rbind(names_to = 'trait')
  }

  if (reference == 'female') {
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

    gwas$chr2 <- fct_relabel(factor(gwas$chr), \(l) ifelse(str_starts(l, '0'), 'Un', l))
    gwas$chr <- str_remove(gwas$chr, 'chr')
    gwas$chr2 <- str_remove(gwas$chr2, 'chr')

    gwas <- mutate(gwas, chr = fct_relevel(factor(chr), as.character(1:23)))

    # for the male reference, we update the chr12 coordinates, following Fraser et al. and relevel the chr
    if (adjust_chr12) {
      gwas <- mutate(gwas, pos = ifelse(chr == '12', update_chr12(pos), pos))
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
  if (qqplots) {
    opa <- par(mfrow = c(2, 1))
    hist(gwas[[pval_column]])
    GWASTools::qqPlot(gwas[[pval_column]], thinThreshold = 4)
    par(opa)
  }

  gwas <- gwas %>%
    left_join(scaff_sizes, join_by(chr)) %>%
    arrange(chr2) |>
    mutate(
      cum_pos = pos + cum_scaff_start,
      p.value.adj = p.adjust(.data[[pval_column]], method = 'fdr')
    )

  # add linkage info for gwas hits, and convenient colors etc for manhattan
  linkage <- list.files(
    str_glue('sequencing/gwas_sommer/snp_inheritance/sommer_snp_inheritance/{reference}/'), full.names = TRUE
  ) |>
    map(data.table::fread) |>
    list_rbind() |>
    mutate(
      source = ifelse(source == '', 'auto', source),
      chr = str_remove(chr, 'chr')
    )

  gwas <- gwas |>
    #SNPs only
    filter(nchar(ref) == 1, nchar(alt) == 1) |>
    mutate(alt = toupper(alt)) |>
    left_join(
      linkage,
      join_by(chr, pos, ref, alt)
    ) |>
    mutate(
      significant = p.value.adj < 0.05,
      source = ifelse(source == '', 'auto', source),

      facet = case_when(
        significant & chr != 'NC_024342.1' & source %in% c('X', 'Y') ~ 'inset',
        TRUE ~ 'standard'
      ) |> factor(levels = c('standard', 'inset')),
      alt_cols = case_when(
        significant & chr != 'NC_024342.1' & source == 'X' ~ 'color1',
        significant & chr != 'NC_024342.1' & source == 'Y' ~ 'color2',
        as.numeric(factor(chr2)) %% 2 == 1 ~ 'color2',
        TRUE ~ 'color1'
      ),
      category = case_when(
        significant ~ source,
        as.numeric(factor(chr2)) %% 2 == 1 ~ 'unsig1',
        TRUE ~ 'unsig2'
      ),
      cum_pos_XYsep = case_when(
        significant & chr != 'NC_024342.1' & source == 'X' ~ 750e6 + rank(.data[[pval_column]])/1e5,
        significant & chr != 'NC_024342.1' & source == 'Y' ~ 765e6 + rank(.data[[pval_column]])/1e5,
        .default = cum_pos,
      )
    ) |>
    arrange(-.data[[pval_column]])

  return(gwas)
}

get_genes <- function(chr, start, end, reference = 'female') {
  require(VariantAnnotation)
  if (reference == 'female') {
    regions <- GRanges(seqnames = chr, ranges = IRanges(start, end, width = end - start + 1))
    ann <- import('sequencing/GCF_000633615.1_Guppy_female_1.0_MT_genomic.gff.gz', which = regions) %>%
      as.data.frame()

    genes <- ann %>% filter(type == 'gene')
    exons <- ann %>% filter(type == 'exon') %>%
      dplyr::select(substart = start, subend = end, gene) %>%
      left_join(dplyr::select(genes, start, end, gene, strand), 'gene')
    genes$gene <- ifelse(genes$gene == 'LOC103476393', 'texim', genes$gene)
    exons$gene <- ifelse(exons$gene == 'LOC103476393', 'texim', exons$gene)
  } else {
    regions <- GRanges(seqnames = paste0('chr', chr), ranges = IRanges(start = range[1], width = diff(range)))
    ann <- import('sequencing/PRET-male-geneID.annotations.gff3', which = regions) %>%
      as.data.frame()
    genes <- ann %>% filter(type == 'gene') %>% rename(gene = Name)
    exons <- ann %>% filter(type == 'exon') %>%
      dplyr::select(substart = start, subend = end, gene = Parent) %>%
      mutate(gene = str_split_i(gene, '\\.', 1)) |>
      left_join(dplyr::select(genes, start, end, gene, strand), 'gene')
  }
  genes <- genes |>
    mutate(start = ifelse(strand == '+', start, end), end = ifelse(strand == '+', end, end - width - 1))
  exons <- exons |>
    mutate(
      width = end - start + 1,
      start = ifelse(strand == '+', start, end), end = ifelse(strand == '+', end, end - width - 1),
      substart = ifelse(strand == '+', substart, subend),
      end = ifelse(strand == '+', subend, subend - width - 1)
    )
  return(list(genes = genes, exons = exons))
}

get_closest_genes <- function(chr, start, end, reference = 'female') {
  require(VariantAnnotation)
  if (reference == 'female') {
    regions <- GRanges(seqnames = chr, ranges = IRanges(start, end, width = end - start + 1))
    ann <- import('sequencing/GCF_000633615.1_Guppy_female_1.0_MT_genomic.gff.gz', which = regions) %>%
      as.data.frame()

    genes <- ann %>% filter(type == 'gene')
    exons <- ann %>% filter(type == 'exon') %>%
      dplyr::select(substart = start, subend = end, gene) %>%
      left_join(dplyr::select(genes, start, end, gene, strand), 'gene')
    genes$gene <- ifelse(genes$gene == 'LOC103476393', 'texim', genes$gene)
    exons$gene <- ifelse(exons$gene == 'LOC103476393', 'texim', exons$gene)
  } else {
    regions <- GRanges(seqnames = paste0('chr', chr), ranges = IRanges(start = range[1], width = diff(range)))
    ann <- import('sequencing/PRET-male-geneID.annotations.gff3', which = regions) %>%
      as.data.frame()
    genes <- ann %>% filter(type == 'gene') %>% rename(gene = Name)
    exons <- ann %>% filter(type == 'exon') %>%
      dplyr::select(substart = start, subend = end, gene = Parent) %>%
      mutate(gene = str_split_i(gene, '\\.', 1)) |>
      left_join(dplyr::select(genes, start, end, gene, strand), 'gene')
  }
  genes <- genes |>
    mutate(start = ifelse(strand == '+', start, end), end = ifelse(strand == '+', end, end - width - 1))
  exons <- exons |>
    mutate(
      width = end - start + 1,
      start = ifelse(strand == '+', start, end), end = ifelse(strand == '+', end, end - width - 1),
      substart = ifelse(strand == '+', substart, subend),
      end = ifelse(strand == '+', subend, subend - width - 1)
    )
  return(list(genes = genes, exons = exons))
}

get_ranges <- function(gwa, width = 1e4) {
  ranges <- data.frame(chr = gwa$chr[1], start = gwa$pos[1] - width, end = gwa$pos[1] + width)
  pointer <- 1
  for (locus in 2:nrow(gwa)) {
    if (gwa$chr[locus] == ranges$chr[pointer]) {
      if (gwa$pos[locus] - width < ranges$end[pointer]) {
        ranges$end[pointer] <- gwa$pos[locus] + width
        next
      }
    }
    # otherwise add new range to ranges
    ranges <- rbind(
      ranges,
      data.frame(chr = gwa$chr[locus], start = gwa$pos[locus] - width, end = gwa$pos[locus] + width)
    )
    pointer <- pointer + 1
  }
  ranges$start <- pmax(ranges$start, 0L)
  return(ranges)
}

reduce_peaks <- function(d, pval_column = 'p_SHet', min_dist = 1e4) {
  # sort the table by p-value, and define the window
  d <- dplyr::arrange(d, .data[[pval_column]]) |>
    dplyr::mutate(.start = pos - min_dist/2, .end = pos + min_dist/2)

  out <- dplyr::slice(d, 1)
  while (TRUE) {
    a <- dplyr::anti_join(dplyr::slice(d, -1), out, join_by(chr, between(x$pos, y$.start, y$.end)))
    if (nrow(a) == nrow(d)) break
    d <- a
    out <- dplyr::bind_rows(out, dplyr::slice(d, 1))
  }

  return(out)
}

make_sommer_manhattan <- function(
    gwas_table = NULL, pval_column = 'p.value', AIC_limit = 0.8, reference = 'female'
) {
  if (!is.null(trait)) {
    if (reference == 'female') {
      linkage_file <- glue::glue('sequencing/gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv')
    } else {
      linkage_file <- glue::glue('sequencing/male_reference_gwas/snp_inheritance/gwas_snp_inheritance/{trait}.csv')
    }
    if (any(gwas_table$significant) && file.exists(linkage_file)) {
      linkage_table <- data.table::fread(linkage_file) %>%
        group_by(chr, start, end, allele1) %>%
        slice_max(AIC_weight) %>%
        ungroup() %>%
        mutate(
          source = ifelse(AIC_weight < AIC_limit, 'Not called', source),
          chr = ifelse(str_starts(chr, 'chr'), parse_number(chr), chr)
        )
      if (reference == 'male') {
        # for the male reference, we update the chr12 coordinates, following Fraser et al. and relevel the chr
        source('sequencing/male_reference_gwas/update_chr12_liftover.R')
        linkage_table <- linkage_table %>% mutate(
          across(c(start, end), \(pos) ifelse(chr == '12', update_chr12(pos), pos))
        )
      }
      gwas_table <- left_join(gwas_table, linkage_table, join_by(chr, pos == start, allele1)) |>
        # relevel for the male reference
        mutate(chr = fct_relevel(factor(chr), as.character(1:23))) |>
        # drop the MT for the female reference
        filter(chr != 'NC_024238.1')
    } else {
      gwas_table$significant <- gwas_table$source <- FALSE
    }
  }

  if (reference == 'female') {
    scaff_labs <- data.table::fread('sequencing/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai') %>%
      dplyr::rename(chr = V1, scaff_len = V2, cum_scaff_start = V3) %>%
      mutate(
        scaff_mid = cum_scaff_start + 0.5 * scaff_len,
        chr2 = case_when(chr == 'NC_024238.1' ~ 'MT',
                         str_starts(chr, 'NC') ~ (str_sub(chr, 8, 9) %>% as.numeric() - 30) %>% as.character(),
                         TRUE ~ 'Un')
      ) %>% filter(!(chr2 %in% c('Un', 'MT')))
  } else {
    scaff_labs <- data.table::fread('sequencing/male_reference_gwas/male_reference.fna.fai') %>%
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
  }

  xscale <- scale_x_continuous(expand = expansion(c(0.01, 0.01)),
                               breaks = scaff_labs$scaff_mid, labels = scaff_labs$chr2)
  yscale <- scale_y_continuous(expand = expansion(c(0.005, 0.05)))

  # add two categories for the non-signficant values (for two tone chr display)
  gwas_table <- mutate(gwas_table, category = case_when(
    significant ~ source,
    !significant & as.numeric(factor(chr2)) %% 2 == 1 ~ 'unsig1',
    TRUE ~ 'unsig2'
  )) |>
    arrange(-.data[[pval_column]])

  ggplot(gwas_table, aes(x = cum_pos, y = -log10(.data[[pval_column]]), color = category)) +
    geom_point(size = 0.6, stroke = 0) +
    xscale + yscale +
    scale_color_manual(
      values = c(
        'auto' = 'black', 'X' = 'firebrick', 'Y' = 'blue3', 'Not called' = 'grey40',
        'unsig1' = 'grey80', 'unsig2' = 'grey90'
      ),
      breaks = c('auto', 'X', 'Y', 'Not called'),
      labels = c('Autosomal', 'X-linked', 'Y-linked', 'Not called'),
      name = 'Pattern of\ninheritance\n(>80% support)',
      limits = c('auto', 'X', 'Y', 'Not called', 'unsig1', 'unsig2')
    ) +
    coord_cartesian(clip = 'off') +
    labs(y = expression(-log[10](italic(p))), x = NULL) +
    theme(
      axis.ticks.x = element_blank(),
      strip.background = element_blank(), strip.placement = 'outside',
      #strip.text.y = element_blank(),
      panel.spacing = unit(2, 'points'),
      panel.background = element_rect(fill = NA) # to help avoid clipping of points
    )
}
