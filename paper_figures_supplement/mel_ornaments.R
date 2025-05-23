source('paper_figures_supplement/subplots/mel_ornament_img.R')
source('paper_figures_supplement/subplots/mel_ornament_presence_size.R')
source('paper_figures_supplement/subplots/mel_ornament_h2.R')
#source('paper_figures_supplement/subplots/mel_ornament_embedding.R')
library(patchwork)

mel_ornaments <- (P_mel_orn_img) /
  (P_mel_orn_h2 + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  (P_mel_orn_pres + theme(legend.position = 'right', legend.justification = c(0, 0))) /
  P_mel_orn_size /
  #(P_mel_orn_emb + theme(legend.position = 'right', legend.justification = c(0, 0.5))) /
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(0.8, 2, 2, 2))
ggsave('paper_figures_supplement/mel_ornaments.png', mel_ornaments, units = "cm", width = 18.4, height = 7.5)
