library(patchwork)
source('paper_figures_supplement/subplots/mel_ornament_img.R')
source('paper_figures_supplement/subplots/mel_ornament_presence_size.R')
source('paper_figures_supplement/subplots/mel_ornament_h2.R')
source('paper_figures/theme.R')

edf4 <- (P_mel_orn_img) /
  (P_mel_orn_h2 + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  (P_mel_orn_pres + theme(legend.position = 'right', legend.justification = c(0, 0))) /
  P_mel_orn_size /
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(0.8, 2, 2, 2))

ggsave('extended_data_figures/edf4.png', edf4, units = "cm", width = 18, height = 8, dpi = 300)
ggsave('extended_data_figures/edf4.pdf', edf4, units = "cm", width = 18, height = 8, device = quartz, type = 'pdf')
