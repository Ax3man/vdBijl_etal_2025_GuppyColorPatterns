library(tidyverse)
library(imager)
source('paper_figures/theme.R')

car_images <- list.files('data/carotenoid_coloration_warped', recursive = TRUE, full.names = TRUE) %>%
  str_subset('Icon', negate = TRUE) %>%
  setNames(., .) %>%
  map(\(f) load.image(f) %>% channel(4), .progress = 'loading carotenoid images')

mel_images <- list.files('data/melanic_coloration_warped_v2', recursive = TRUE, full.names = TRUE) %>%
  str_subset('Icon', negate = TRUE) %>%
  setNames(., .) %>%
  map(\(f) load.image(f) %>% channel(4), .progress = 'loading melanin images')

safe_average <- function(unsafe_list) {
  cimgs <- map_lgl(unsafe_list, is.cimg)
  if (any(cimgs)) return(average(unsafe_list[cimgs]))
  return(NA)
}

average_patterns <- data.table::fread('photo_database.csv') %>%
  mutate(
    car_file = glue::glue('data/carotenoid_coloration_warped/{replicate}/{generation}/{unique_id}.png'),
    mel_file = glue::glue('data/melanic_coloration_warped_v2/{replicate}/{generation}/{unique_id}.png')
  ) %>%
  group_by(fish_id, facing_direction) %>%
  summarise(
    average_car = list(safe_average(car_images[car_file])),
    average_mel = list(safe_average(mel_images[mel_file])),
    .groups = 'drop_last'
  ) %>%
  summarise(
    average_car = list(safe_average(average_car)),
    average_mel = list(safe_average(average_mel)),
    .groups = 'drop'
  ) %>%
  summarise(
    average_car = list(safe_average(average_car)),
    average_mel = list(safe_average(average_mel))
  )

car_heatmap_data <- as.data.frame(average_patterns$average_car[[1]]) %>% as_tibble()
mel_heatmap_data <- as.data.frame(average_patterns$average_mel[[1]]) %>% as_tibble()

bg <- load.image('data/extracted_fish_warped/replicate_3/gen_4/20221024_IMG_2089.png') %>%
  channel(4) %>% as.data.frame() %>% filter(value > 0)

p <- ggplot(mapping = aes(x, y, fill = value)) +
  geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
  geom_raster() +
  scale_y_reverse() +
  coord_fixed(expand = FALSE) +
  theme_void(base_size = 7, base_family = 'Arial') +
  theme(
    legend.position = 'bottom', legend.title = element_text(vjust = 0.8),
    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.8, 'cm'),
    legend.justification = 'right'
  )

p_car <- p %+%
  car_heatmap_data +
  #scale_fill_gradient(low = 'grey80', high = 'darkorange', limits = c(0.01, NA), na.value = 'transparent') +
  scale_fill_viridis_c(
    option = 'A',
    limits = c(0.01, NA),
    na.value = 'transparent',
    name = 'Incidence of orange (%)',
    labels = scales::label_percent(),
    breaks = c(0.01, 0.25, 0.5, 0.75)
  ) +
  labs(tag = 'C')

p_mel <- p %+%
  mel_heatmap_data +
  #scale_fill_gradient(low = 'white', high = 'black', limits = c(0.01, NA), na.value = 'transparent') +
  scale_fill_viridis_c(
    option = 'G',
    limits = c(0.01, NA),
    na.value = 'transparent',
    name = 'Incidence of black (%)',
    labels = scales::label_percent(),
    breaks = c(0.01, 0.15, 0.3, 0.45)
  )
pop_incidence_heatmaps <- p_car / p_mel

# artful minis for fig 5
mini <- ggplot(mapping = aes(x, y, fill = value)) +
  geom_raster(aes(alpha = value), data = bg, fill = 'grey70', show.legend = FALSE) +
  geom_raster(show.legend = FALSE) +
  scale_y_reverse() +
  coord_fixed(expand = FALSE) +
  theme_void(base_size = 7, base_family = 'Arial')

mini %+% filter(car_heatmap_data, value > 0.01) +
  scale_fill_gradient(low = 'orange', high = 'white')

mini %+% filter(mel_heatmap_data, value > 0.01) +
  scale_fill_gradient(low = 'black', high = 'white')
