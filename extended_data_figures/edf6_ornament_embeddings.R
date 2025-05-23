library(tidyverse)

source('paper_figures/subplots/car_ornament_img.R')
source('paper_figures/subplots/car_ornament_embedding.R')

source('paper_figures_supplement/subplots/mel_ornament_img.R')
source('paper_figures_supplement/subplots/mel_ornament_embedding.R')

car <- (P_car_orn_img) /
  (P_car_orn_emb + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  plot_layout(heights = c(0.6, 2))

mel <- (P_mel_orn_img) /
  (P_mel_orn_emb + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  plot_layout(heights = c(0.6, 2))

edf6 <- wrap_plots(car, mel, nrow = 2) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom') &
  scale_color_viridis_c(
    na.value = 'grey80', trans = 'sqrt', breaks = c(.01, .1, .3, .6, 1),
    guide = guide_colorbar(barwidth = unit(8, 'cm'))
  )
ggsave('extended_data_figures/edf6.png', edf6, units = "cm", width = 18, height = 10, dpi = 300)
ggsave('extended_data_figures/edf6.pdf', edf6, units = "cm", width = 18, height = 10, device = quartz, type = 'pdf')
