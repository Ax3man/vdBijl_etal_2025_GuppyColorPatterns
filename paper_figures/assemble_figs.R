library(patchwork)

source('paper_figures/subplots/1_phenotyping_pipeline.R')
source('paper_figures/subplots/1_random_patrilines.R')
source('paper_figures/subplots/1_population_incidence_heatmaps.R')
fig1A <- wrap_elements(full = phenotyping_pipeline + labs(tag = 'A'))
fig1B <- ((example_patrilines + labs(tag = 'B')) | pop_incidence_heatmaps) + plot_layout(widths = c(4, 3))
fig1 <- (fig1A / fig1B) + plot_layout(heights = c(1, 2.2))
ggsave('paper_figures/Fig1.png', fig1, units = 'cm', width = 18, height = 9)
ggsave('paper_figures/Fig1.pdf', fig1, units = 'cm', width = 18, height = 9, device = cairo_pdf)

source('paper_figures/subplots/car_pedigree.R')
source('paper_figures/subplots/car_selection_maps.R')
source('paper_figures/subplots/embedding_selection.R')

Fig2 <- wrap_elements(full = P_car_pedigree) +
  (wrap_elements(full = P_car_selmap) /
     wrap_elements(embedding_selection) +
                  plot_layout(heights = c(1.6, 1))) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = 'A')
ggsave('paper_figures/Fig2.png', Fig2, units = "cm", width = 18, height = 10)
ggsave('paper_figures/Fig2.pdf', Fig2, units = "cm", width = 18, height = 10, device = cairo_pdf)

source('paper_figures/subplots/car_pixel_h2.R')
source('paper_figures/subplots/mel_pixel_h2.R')
Fig3 <- (
  (P_car_pixel_h2) +
  (P_mel_pixel_h2 + theme(strip.text.y.left = element_blank())) &
  theme(legend.position = 'bottom')
) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect')
ggsave('paper_figures/Fig3.png', Fig3, units = "cm", width = 8.8, height = 5)
ggsave('paper_figures/Fig3.pdf', Fig3, units = "cm", width = 8.8, height = 5, device = cairo_pdf)

source('paper_figures/subplots/car_ornament_img.R')
source('paper_figures/subplots/car_ornament_presence_size.R')
source('paper_figures/subplots/car_ornament_h2.R')
#source('paper_figures/subplots/car_ornament_embedding.R')
Fig4 <- (P_car_orn_img) /
  (P_car_orn_h2 + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  (P_car_orn_pres + theme(legend.position = 'right', legend.justification = c(0, 0))) /
  P_car_orn_size /
  #(P_car_orn_emb + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(0.4, 1, 1, 1))
ggsave('paper_figures/Fig4.png', Fig4, units = "cm", width = 18, height = 7.5)
ggsave('paper_figures/Fig4.pdf', Fig4, units = "cm", width = 18, height = 7.5, device = cairo_pdf)
ggsave('paper_figures/Fig4_alt.pdf', Fig4, units = "cm", width = 18, height = 7.5)
ggsave('paper_figures/Fig4.svg', Fig4, units = "cm", width = 18, height = 7.5)

# Note: Figure 5 takes a while to plot, so it's done in subplots/5_manhattan_offtarget.R
# to remake run
# source('paper_figures/subplots/5_manhattan_offtarget.R')

source('paper_figures/subplots/6_texim_prex1.R')

ggsave('paper_figures/Fig6.png', Fig6, units = "cm", width = 8.8, height = 12)
ggsave('paper_figures/Fig6.pdf', Fig6, units = "cm", width = 8.8, height = 12, device = cairo_pdf)
