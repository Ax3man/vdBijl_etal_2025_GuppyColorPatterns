library(tidyverse)
library(patchwork)
source('sequencing/genomics_helpers.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')
source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('paper_figures/theme.R')
source('sequencing/gwas_sommer/ornament_peak_figure/make_minis.R')

pheno <- get_phenotypes2()

car_orn_gwas <- prep_sommer_gwas_table(
  str_glue('car_{c(1:3, 5:7)}'), 'orange_ornaments', pval_column = 'p.value'
)
mel_orn_gwas <- prep_sommer_gwas_table(
  str_glue('mel_{c(1, 3, 5:8)}'), 'black_ornaments', pval_column = 'p.value'
)
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

triplet <- function(ornament, chosen_chr, snp_effects, gwas, viz_width = 1e5, skip = 0, present_color = 'orange3') {
  d <- gwas |> filter(trait == ornament, chr2 != 'Un')
  d_sig <- d |> filter(significant)

  #ggplot(d_sig, aes(cum_pos, -log10(p.value), color = source)) + geom_point()
  #d_sig |> filter(chr2 == chosen_chr) |>
  # ggplot(aes(pos / 1e6, -log10(p.value), color = source)) + geom_point()

  if (skip == 0) {
    choice <- slice_min(d_sig |> filter(chr2 == chosen_chr), p.value, with_ties = FALSE)
  } else {
    choice <- slice_min(d_sig |> filter(chr2 == chosen_chr), p.value, n = skip + 1) |>
      slice_max(p.value, with_ties = FALSE)
  }
  genes <- get_genes(choice$chr, choice$pos - 0.5 * viz_width, choice$pos + 0.5 * viz_width)

  d2 <- d |>
    filter(chr2 == chosen_chr, between(pos, choice$pos - viz_width/2, choice$pos + viz_width/2))

  p1 <- ggplot(d2, aes(pos / 1e6, -log10(p.value), color = significant)) +
    geom_segment(
      aes(
        start / 1e6, xend = end / 1e6,
        max(-log10(d2$p.value)) * (1.1 + 0.05 * seq_len(nrow(genes$genes)) %% 3), yend = after_scale(y)
      ),
      genes$genes,
      color = 'black'#,
      #arrow = grid::arrow(length = unit(0.05, 'inch')), lineend = 'round', linejoin = 'mitre'
    ) +
    ggpp::geom_text_s(
      aes(
        (start + end) / 2e6,
        max(-log10(d2$p.value)) * (1.1 + 0.05 * seq_len(nrow(genes$genes)) %% 3),
        label = ifelse(str_starts(gene, 'LOC'), '', gene)
      ),
      genes$genes,
      position = ggpp::position_nudge_to(
        x = (range(d2$pos) + diff(range(d2$pos)) * c(0.1, -0.2)) / 1e6,
        y = max(-log10(d2$p.value)) * 1.25,
        x.action = 'spread'
      ),
      default.color = 'black', color = 'grey80', color.target = 'segment',
      size = 1.5, angle = 30, hjust = 0, segment.linewidth = 0.1
    ) +
    # geom_text(
    # aes(
    #   (start + end) / 2e6, max(-log10(d2$p.value)) * (1.1 + 0.05 * seq_len(nrow(genes$genes))),
    #   label = ifelse(str_starts(gene, 'LOC'), '', gene)
    # ),
    # genes$genes,
    # color = 'black', vjust = -0.5, size = 1.5
    # ) +
    geom_point(size = 0.3) +
    geom_point(data = choice, color = 'firebrick', size = 0.3) +
    scale_color_manual(values = c('grey60', 'black'), guide = 'none') +
    scale_x_continuous(
      expand = c(0, 0),
      name = str_glue('Chr {choice$chr2} (Mb)'),
      breaks = scales::extended_breaks(3)
    ) +
    scale_y_continuous(expand = expansion(c(0, 0.2)), name = expression(-log[10](italic(p)))) +
    coord_cartesian(xlim = range(d2$pos) / 1e6)

  g <- load_vcf(region = str_glue_data(choice, '{chr}:{pos}-{pos}'))$mat |>
    as.data.frame() |> rownames_to_column('sample_name') |>
    left_join(pheno, join_by(sample_name))
  m <- glm(reformulate('V1', ornament), data = g, family = 'binomial')

  r <- choice$ref; a <- choice$alt
  labs <- paste(c(r, r, a), c(r, a, a), sep = '/')

  p2 <- ggplot(g, aes(
    x = factor(V1, levels = c(-1, 0, 1)),
    fill = factor(.data[[ornament]], c(0, 1), c('Absent', 'Present')))
  ) +
    geom_bar(position = 'fill') +
    annotate(
      'text', x = '0', y = 1, vjust = -0.2, size = 2,
      label = str_glue('R^2 == {round(performance::r2_nagelkerke(m), 2)}'), parse = TRUE
    ) +
    scale_x_discrete(limits = c('-1', '0', '1'), labels = labs) +
    scale_fill_manual(values = c(present_color, 'white'), guide = 'none',
                      limits = c('Present', 'Absent')) +
    scale_y_continuous(expand = expansion(c(0, 0.2)), name = 'Proportion w/ ornament', breaks = c(0, 0.5, 1)) +
    labs(x = 'Genotype')
  p3 <- make_mini(snp_effects, choice$chr, choice$pos, 4) +
    labs(fill = 'Z-score') +
    theme(text = my_theme$text)
  p1 / p2 / p3
}

#t1 <- triplet('car_1', '12', orange_snp_effects, car_orn_gwas)
t2 <- triplet('car_2', '5', orange_snp_effects, car_orn_gwas)
t3 <- triplet('car_3', '5', orange_snp_effects, car_orn_gwas)
t5 <- triplet('car_5', '13', orange_snp_effects, car_orn_gwas, skip = 1)
#t6 <- triplet('car_6', '14', orange_snp_effects, car_orn_gwas)
t7 <- triplet('car_7', '9', orange_snp_effects, car_orn_gwas)
#(t1 | t2 | t3 | t5 | t6 | t7) + plot_layout(guide = 'collect')

# car_triplets <- (t2 | t3 | t5 | t7) + plot_layout(guide = 'collect')
# ggsave('sequencing/gwas_sommer/ornament_peak_figure/car_ornaments.png', car_triplets,
#        width = 18.4, height = 6, units = 'cm', dpi = 600)

b1 <- triplet('mel_1', '11', black_snp_effects, mel_orn_gwas, present_color = 'black')
b3 <- triplet('mel_3', '15', black_snp_effects, mel_orn_gwas, present_color = 'black')
#b5 <- triplet('mel_5', '12', black_snp_effects, mel_orn_gwas, present_color = 'black')
b6 <- triplet('mel_6', '10', black_snp_effects, mel_orn_gwas, present_color = 'black')
#b7 <- triplet('mel_7', '7', black_snp_effects, mel_orn_gwas, present_color = 'black')
b8 <- triplet('mel_8', '2', black_snp_effects, mel_orn_gwas, present_color = 'black')
# (b1 | b3 | b5 | b6 | b7 | b8) + plot_layout(guide = 'collect')

# mel_triplets <- (b1  | b3 | b6 | b8) + plot_layout(guide = 'collect')
# ggsave('sequencing/gwas_sommer/ornament_peak_figure/mel_ornaments.png', mel_triplets,
#        width = 18.4, height = 6, units = 'cm', dpi = 600)

# Make one plot for the 8 common ornaments with appreciable autosomal contribution:
l <- theme(axis.title.y = element_blank())
complete <- ((t2 | (t3 & l) | (t5 & l) | (t7 & l) | (b1 & l) | (b3 & l) | (b6 & l) | (b8 & l)) &
               theme(legend.position = 'bottom',
                     legend.key.width = unit(1, 'cm'),
                     legend.key.height = unit(0.2, 'cm'))) +
  plot_layout(guide = 'collect')
ggsave('sequencing/gwas_sommer/ornament_peak_figure/ornaments.png', complete,
       width = 18.4, height = 9, units = 'cm', dpi = 600)
ggsave('sequencing/gwas_sommer/ornament_peak_figure/ornaments.pdf', complete, device = cairo_pdf,
       width = 18.4, height = 9, units = 'cm')

if (FALSE) {
  require(tidyverse)
  require(imager)
  eb <- element_blank()

  bg <- load.image('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
    channel(4) %>%
    as.data.frame()
  pixsets <- c(
    list.files('ornament_analysis/ornament_images', 'car_._new', full.names = TRUE),
    list.files('ornament_analysis/ornament_images', 'mel_.', full.names = TRUE)
  ) |>
    map(\(.x) load.image(.x) |> channel(4) |> as.data.frame()) |>
    setNames(c(paste0('O', 1:7), paste0('B', 1:8))) |>
    bind_rows(.id = 'ornament') |>
    mutate(type = ifelse(str_starts(ornament, 'O'), 'orange', 'black'))

  orn_img <- ggplot(mapping = aes(x, y, alpha = value)) +
    geom_raster(data = bg, fill = 'grey60') +
    geom_raster(
      aes(fill = type),
      filter(pixsets, ornament %in% c('O2', 'O3', 'O5', 'O7', 'B1', 'B3', 'B6', 'B8')),
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
  ggsave('sequencing/gwas_sommer/ornament_peak_figure/ornament_labels.pdf', orn_img, device = cairo_pdf,
         width = 18.4, height = 2, units = 'cm')
}

