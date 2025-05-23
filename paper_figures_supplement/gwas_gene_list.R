library(tidyverse)
library(VariantAnnotation)
library(gggenes)
source('sequencing/gwas/peak_viz/panel_functions/make_gene_panel.R')

ann <- import('sequencing/GCF_000633615.1_Guppy_female_1.0_MT_genomic.gff.gz') %>%
  as.data.frame()

genes <- ann %>% filter(type == 'gene')
exons <- ann %>% filter(type == 'exon') %>%
  dplyr::select(substart = start, subend = end, gene) %>%
  left_join(dplyr::select(genes, seqnames, start, end, gene, strand), 'gene')


gwas <- bind_rows(
  'black pattern' = data.table::fread('sequencing/gwas_sommer/sommer_output/female_mel_PIE.csv.gz'),
  'orange pattern' = data.table::fread('sequencing/gwas_sommer/sommer_output/female_car_PIE.csv.gz'),
  .id = 'color'
)

ROI <- tribble(
  ~label, ~gene, ~chr, ~start, ~end,
  'csf1ra (LG10)', 'csf1r', 'NC_024340.1', 25495046, 25511233,
  'kita (LG4)', 'kit', 'NC_024334.1', 19488737, 19514083,
  'adcy5 (LG2)', 'adcy5', 'NC_024332.1', 13712347, 13798894,
  'aldh1a1 (LG4)', 'aldh1a1', 'NC_024334.1', 14786261, 14791869,
  # extra copy of csf1r?
  'csf1r-2 (Un)', 'csf1r', 'NW_007615014.1', 895156, 910608
) %>%
  mutate(start_wide = start - 10000, end_wide = end + 10000)

map(1:nrow(ROI), \(i) {
  roi <- ROI[i, ]
  d <- semi_join(gwas, roi, join_by(chr, between(pos, start_wide, end_wide)))
  g <- semi_join(genes, roi, join_by(seqnames == chr, overlaps(start, end, start_wide, end_wide))) %>%
    mutate(id = row_number())
  e <- semi_join(exons, roi, join_by(seqnames == chr, overlaps(start, end, start_wide, end_wide))) %>%
    inner_join(g)

  ggplot(d) +
    geom_gene_arrow(
      aes(xmin = start, xmax = end, y = 5 - id * 0.4, forward = strand == '+'),
      g, fill = 'grey90'
    ) +
    geom_subgene_arrow(
      aes(xmin = start, xmax = end, xsubmin = substart, xsubmax = subend, y = 5 - id * 0.4, forward = strand == '+'),
      e, fill = 'grey70', color = NA
    ) +
     geom_text(
       aes(x = (start + end) / 2, y = 5 - id * 0.4, label = gene),
       g, size = 3
    ) +
    geom_point(aes(pos, -log10(p_SHet))) +
    facet_grid(cols = vars(color)) +
    coord_cartesian(xlim = c(roi$start_wide, roi$end_wide)) +
    scale_x_continuous(labels = scales::label_comma(scale = 1e-3)) +
    labs(
      subtitle = roi$label,
      x = 'position (kb)',
      y = expression(-log[10](italic(P)))
    ) +
    theme(strip.background = element_rect(colour = NULL, fill = 'grey90', linewidth = 0))
}) %>%
  patchwork::wrap_plots(ncol = 1)
# 4 warnings that subgenes breaking boundaries are expected because of the duplicated gene name csf1r
ggsave('paper_figures_supplement/gwas_gene_list.png', w = 10, h = 12)


# car_ROI <- tribble(
#   ~label, ~gene, ~chr, ~start, ~end,
#   'ALDH1A1', 'ALDH1A1', 'NC_024334.1', 14786261, 14791869
# ) %>%
#   mutate(start_wide = start - 10000, end_wide = end + 10000) #%>%
# #dplyr::select(-start, -end)
#
# map(1:nrow(car_ROI), \(i) {
#   roi <- car_ROI[i, ]
#   d <- semi_join(gwas, roi, join_by(chr, between(ps, start_wide, end_wide)))
#   g <- semi_join(genes, roi, join_by(seqnames == chr, overlaps(start, end, start_wide, end_wide))) %>%
#     mutate(id = row_number())
#   e <- semi_join(exons, roi, join_by(seqnames == chr, overlaps(start, end, start_wide, end_wide))) %>%
#     inner_join(g)
#
#   ggplot(d) +
#     geom_point(aes(ps, -log10(p_SHet))) +
#     geom_gene_arrow(
#       aes(xmin = start, xmax = end, y = 2 - id * 0.1, forward = strand == '+'),
#       g, fill = 'grey90'
#     ) +
#     geom_subgene_arrow(
#       aes(xmin = start, xmax = end, xsubmin = substart, xsubmax = subend, y = 2 - id * 0.1, forward = strand == '+'),
#       e, fill = 'grey70', color = NA
#     ) +
#     geom_text(
#       aes(x = (start + end) / 2, y = 2 - id * 0.1, label = gene),
#       g, size = 3
#     ) +
#     facet_grid(cols = vars(color)) +
#     coord_cartesian(xlim = c(roi$start_wide, roi$end_wide)) +
#     scale_x_continuous(labels = scales::label_comma(scale = 1e-3)) +
#     labs(
#       subtitle = roi$label,
#       x = 'position (kb)',
#       y = expression(-log[10](italic(P)))
#     )
# }) %>%
#   patchwork::wrap_plots(ncol = 1)
