require(tidyverse)
require(imager)
source('paper_figures/theme.R')
eb <- element_blank()

bg <- load.image('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
  channel(4) %>%
  as.data.frame()
pixsets <- list.files('ornament_analysis/ornament_images', 'mel_', full.names = TRUE) %>%
  map(\(.x) load.image(.x) %>% channel(4) %>% as.data.frame) %>%
  setNames(paste0('mel_', seq_along(.))) %>%
  bind_rows(.id = 'ornament') %>%
  mutate(ornament_lab = fct_relabel(ornament, \(l) str_replace(l, 'mel_', 'B')))

P_mel_orn_img <- ggplot(mapping = aes(x, y, alpha = value)) +
  geom_raster(data = bg, fill = 'grey60') +
  geom_raster(data = pixsets, fill = 'black') +
  scale_alpha_identity(guide = 'none') +
  scale_y_reverse() +
  facet_grid(~ornament_lab) +
  coord_fixed(expand = FALSE) +
  theme(
    axis.text = eb, axis.line = eb, axis.title = eb, axis.ticks = eb,
    strip.background = eb
  )
