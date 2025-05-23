library(tidyverse)
library(patchwork)
source('sequencing/genomics_helpers.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')
source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('paper_figures/theme.R')
source('sequencing/gwas_sommer/ornament_peak_figure/make_minis.R')

pheno <- get_phenotypes2()

car_gwas <- prep_sommer_gwas_table('car_PIE', pval_column = 'p_SHet')
mel_gwas <- prep_sommer_gwas_table('mel_PIE', pval_column = 'p_SHet')

if (file.exists('sequencing/gwas_sommer/snp_effects/orange_ornaments.csv.gz')) {
  orange_snp_effects <- data.table::fread('sequencing/gwas_sommer/snp_effects/orange_ornaments.csv.gz')
} else {
  orange_snp_effects <- prep_for_minis('orange_ornaments')
  data.table::fwrite(orange_snp_effects, 'sequencing/gwas_sommer/snp_effects/orange_ornaments.csv.gz')
}
if (file.exists('sequencing/gwas_sommer/snp_effects/black_ornaments.csv.gz')) {
  black_snp_effects <- data.table::fread('sequencing/gwas_sommer/snp_effects/black_ornaments.csv.gz')
} else {
  black_snp_effects <- prep_for_minis('black_ornaments')
  data.table::fwrite(black_snp_effects, 'sequencing/gwas_sommer/snp_effects/black_ornaments.csv.gz')
}

triplet <- function(ornament, chosen_chr, title, snp_effects, gwas,
                    viz_width = 6e4, skip = 0, present_color = 'orange3') {
  d <- gwas |> filter(chr2 != 'Un')
  d_sig <- d |> filter(significant)

  ggplot(d_sig, aes(cum_pos, -log10(p_SHet), color = source)) + geom_point() +
    facet_wrap(vars(chr2), scales = 'free_x')
  d_sig |> filter(chr2 == chosen_chr) |>
    ggplot(aes(pos / 1e6, -log10(p_SHet), color = source)) + geom_point()

  if (skip == 0) {
    choice <- slice_min(d_sig |> filter(chr2 == chosen_chr), p_SHet, with_ties = FALSE)
  } else {
    choice <- slice_min(d_sig |> filter(chr2 == chosen_chr), p_SHet, n = skip + 1) |>
      mutate(p_SHet = ifelse(nchar(ref) + nchar(alt) == 2, p_SHet, p_SHet - 1e-60)) |>
      slice_max(p_SHet, with_ties = FALSE)
  }
  genes <- get_genes(choice$chr, choice$pos - 0.5 * viz_width, choice$pos + 0.5 * viz_width)

  d2 <- d |>
    filter(chr2 == chosen_chr, between(pos, choice$pos - viz_width/2, choice$pos + viz_width/2))

  p1 <- ggplot(d2, aes(pos / 1e6, -log10(p_SHet), color = significant)) +
    geom_segment(
      aes(
        start / 1e6, xend = end / 1e6,
        max(-log10(d2$p_SHet)) * (1.1 + 0.05 * seq_len(nrow(genes$genes)) %% 3), yend = after_scale(y)
      ),
      genes$genes,
      color = 'black'#,
      #arrow = grid::arrow(length = unit(0.05, 'inch')), lineend = 'round', linejoin = 'mitre'
    ) +
    ggpp::geom_text_s(
      aes(
        (start + end) / 2e6,
        max(-log10(d2$p_SHet)) * (1.1 + 0.05 * seq_len(nrow(genes$genes)) %% 3),
        label = ifelse(str_starts(gene, 'LOC'), '', gene)
      ),
      genes$genes,
      position = ggpp::position_nudge_to(
        x = (range(d2$pos) + diff(range(d2$pos)) * c(0.1, -0.2)) / 1e6,
        y = max(-log10(d2$p_SHet)) * 1.25,
        x.action = 'spread'
      ),
      default.color = 'black', color = 'grey80', color.target = 'segment',
      size = 1.5, angle = 30, hjust = 0, segment.linewidth = 0.1
    ) +
    geom_point(size = 0.3) +
    geom_point(data = choice, color = 'black', shape = 23, size = 0.7) +
    scale_color_manual(values = c('grey60', 'black'), guide = 'none') +
    scale_x_continuous(
      expand = c(0, 0),
      name = str_glue('Chr {choice$chr2} (Mb)'),
      breaks = scales::extended_breaks(3)
    ) +
    facet_wrap(vars('a'), labeller = \(x) title) +
    scale_y_continuous(expand = expansion(c(0, 0.2)), name = expression(-log[10](italic(p)))) +
    coord_cartesian(xlim = range(d2$pos) / 1e6) +
    theme(strip.background = element_rect(fill = 'grey90', color = NA))

  cov <- with(choice, get_region_coverage(chr, pos, c(pos - viz_width, pos + viz_width)))
  femcov <- with(choice, get_region_coverage_yuying(chr, pos, c(pos - viz_width, pos + viz_width))) |>
    filter(sex == 'female')

  p2 <- ggplot(cov, aes((window_start + window_end) / 2e6, rel_coverage, group = sample_name)) +
    geom_line(linewidth = 0.1, alpha = 0.1) +
    geom_line(data = femcov, linewidth = 0.3, alpha = 1, color = 'firebrick') +
    scale_x_continuous(
      expand = c(0, 0),
      name = str_glue('Chr {choice$chr2} (Mb)'),
      breaks = scales::extended_breaks(3)
    ) +
    scale_y_continuous(expand = expansion(c(0, 0.01)), name = 'Relative depth') +
    coord_cartesian(xlim = range(d2$pos) / 1e6)

  g <- load_vcf(region = str_glue_data(choice, '{chr}:{pos}-{pos}'))$mat |>
    as.data.frame() |> rownames_to_column('sample_name') |>
    left_join(pheno, join_by(sample_name)) |>
    dplyr::rename(GT = V1)

  gd <- get_genotypes('sequencing/gwas_sommer/filtered2.vcf.gz', choice$chr, c(choice$pos - 1, choice$pos + 1)) |>
    dplyr::select(sample_name = sampleNames, refDepth, altDepth) |>
    left_join(get_mean_coverage(), join_by(sample_name)) |>
    mutate(ref_depth = refDepth / mean_coverage, alt_depth = altDepth / mean_coverage) |>
    left_join(g, join_by(sample_name))

  m <- glm(reformulate('GT', ornament), data = g, family = 'binomial')

  r <- choice$ref; a <- choice$alt
  labs <- paste(c(r, r, a), c(r, a, a), sep = '/')

  p3 <- ggplot(g, aes(
    x = factor(GT, levels = c(-1, 0, 1)),
    fill = factor(.data[[ornament]], c(0, 1), c('Absent', 'Present')))
  ) +
    geom_bar(position = 'fill') +
    annotate(
      'text', x = 1.5, y = 1, vjust = -0.2, size = 2,
      label = str_glue('R^2 == {round(performance::r2_nagelkerke(m), 2)}'), parse = TRUE
    ) +
    scale_x_discrete(
      limits = c('-1', '0', '1')[c(-1, 0, 1) %in% g$GT],
      labels = labs[c(-1, 0, 1) %in% g$GT]
    ) +
    scale_fill_manual(values = c(present_color, 'white'), guide = 'none',
                      limits = c('Present', 'Absent')) +
    scale_y_continuous(expand = expansion(c(0, 0.2)), name = 'Proportion w/ ornament', breaks = c(0, 0.5, 1)) +
    labs(x = 'Genotype')

  p4 <- gd |>
    pivot_longer(ref_depth:alt_depth, names_to = 'allele', values_to = 'allelic_depth') |>
    ggplot(aes(
      x = factor(GT, levels = c(-1, 0, 1)),
      y = allelic_depth * 2,
      fill = allele
    )) +
    geom_bar(position = 'stack', stat = 'summary', fun = 'median') +
    scale_x_discrete(
      name = 'Genotype',
      limits = c('-1', '0', '1')[c(-1, 0, 1) %in% g$GT],
      labels = labs[c(-1, 0, 1) %in% g$GT]
    ) +
    scale_fill_manual(
      limits = c('ref_depth', 'alt_depth'),
      labels = c('reference', 'alternate'),
      values = c('grey20', 'lightblue'),
      name = 'Allele'
    ) +
    scale_y_continuous(
      name = 'Est. allele count',
      expand = expansion(mult = c(0, 0.05))
    )

  p5 <- make_mini(snp_effects, choice$chr, choice$pos, 4) +
    labs(fill = 'Z-score') +
    theme(text = my_theme$text)
  p1 / p2 / p3 / p4 / p5
}

t6 <- triplet('car_6', '14', 'Orange pattern', orange_snp_effects, car_gwas)
b5 <- triplet('mel_5', '7', 'Black pattern', black_snp_effects, mel_gwas, present_color = 'black')
b7 <- triplet('mel_7', '7', 'Black pattern', black_snp_effects, mel_gwas, present_color = 'black')

l <- theme(axis.title.y = element_blank())
row1 <- (t6[[1]] | (b5[[1]] & l))
row2 <- (t6[[2]] | (b5[[2]] & l))
row3 <- (t6[[3]] | (b5[[3]] + l) | (b7[[3]] + l)) +
  plot_layout(widths = c(1, .5, .5))
row4 <- ((t6[[4]] | (b5[[4]] & l)) & theme(legend.position = 'top')) +
  plot_layout(guide = 'collect')
row5 <- ((t6[[5]] | b5[[5]]) &
  theme(legend.position = 'bottom',
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.2, 'cm'))) +
  plot_layout(guide = 'collect')

complete <- row1 / row2 / row3 / row4 / row5

ggsave('sequencing/gwas_sommer/ornament_peak_figure/ornaments_sexlinked.png', complete,
       width = 18.4/2, height = 12, units = 'cm', dpi = 600)
ggsave('sequencing/gwas_sommer/ornament_peak_figure/ornaments_sexlinked.pdf', complete, device = cairo_pdf,
       width = 18.4/2, height = 12, units = 'cm')

if (FALSE) {
  require(tidyverse)
  require(imager)
  eb <- element_blank()

  bg <- load.image('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
    imager::channel(4) %>%
    as.data.frame()
  pixsets <- c(
    list.files('ornament_analysis/ornament_images', 'car_._new', full.names = TRUE),
    list.files('ornament_analysis/ornament_images', 'mel_.', full.names = TRUE)
  ) |>
    map(\(.x) load.image(.x) |> imager::channel(4) |> as.data.frame()) |>
    setNames(c(paste0('O', 1:7), paste0('B', 1:8))) |>
    bind_rows(.id = 'ornament') |>
    mutate(type = ifelse(str_starts(ornament, 'O'), 'orange', 'black'))

  orn_img <- ggplot(mapping = aes(x, y, alpha = value)) +
    geom_raster(data = bg, fill = 'grey60') +
    geom_raster(
      aes(fill = type),
      filter(pixsets, ornament %in% c('O6', 'B5', 'B7')),
    ) +
    scale_alpha_identity(guide = 'none') +
    scale_fill_identity(guide = 'none') +
    scale_y_reverse() +
    facet_grid(cols = vars(fct_inorder(ornament))) +
    coord_fixed(expand = FALSE) +
    theme(
      axis.text = eb, axis.line = eb, axis.title = eb, axis.ticks = eb,
      strip.background = eb
    )
  ggsave('sequencing/gwas_sommer/ornament_peak_figure/ornament_labels_sexlinked.pdf', orn_img, device = cairo_pdf,
         width = 18.4, height = 2, units = 'cm')
}

