library(tidyverse)
library(posterior)
source('paper_figures/theme.R')

bg_fish <- imager::load.image('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
  #imager::resize_halfXY() %>%
  as.data.frame(wide = 'c') %>% filter(c.4 > 0.1) %>% dplyr::select(x, y, alpha = c.4) %>%
  mutate(x = (x - 1) / 2 + 1, y = (y - 1) / 2 + 1) # correct for half size

posts <- tibble(file = c(
  list.files('quant_gen/cluster_brms_carotenoid/posteriors_1', f = TRUE),
  list.files('quant_gen/cluster_brms_carotenoid/posteriors_2', f = TRUE)
)) %>%
  mutate(file2 = basename(file)) %>%
  separate(file2, c(NA, 'x', 'y', NA)) %>%
  mutate(
    across(c(x, y), parse_number),
    post = map(file, read_rds)
  )

posts_draws <- posts %>%
  mutate(post = map(post, as_draws_df)) %>%
  unnest(post) %>%
  dplyr::rename(
    Ve_id = 'sd_1[1]',
    Va_auto = 'sd_2[1]',
    Va_X = 'sd_3[1]',
    Ve_clutch = 'sd_4[1]',
    Va_Y = 'sd_5[1]'
  ) %>%
  mutate(
    # Note that the model gives us SD, but we want variances (we already named them variances)
    across(Ve_id:Va_Y, \(SD) SD ^ 2),

    Vr = (pi ^ 2) / 3, #logit link variance

    Va = Va_auto + Va_X + Va_Y,
    Ve = Ve_clutch + Ve_id,
    Vtotal = Va + Ve + Vr,

    h2 = Va / (Va + Ve), # residual variance? something with pi
    h2_auto = Va_auto / (Va + Ve),
    h2_X = Va_X / (Va + Ve),
    h2_Y = Va_Y / (Va + Ve),

    rel_id_effect = Ve_id / (Va + Ve),
    rel_clutch_effect = Ve_clutch / (Va + Ve),

    measurement_error = Vr / Vtotal
  )

facet_levels <- c(
  'h2', 'h2_auto', 'rel_clutch_effect', 'h2_X', 'rel_id_effect', 'h2_Y', 'measurement_error'
)
facet_labels <- c(
  expression(italic(h)^2), expression(italic(h)[autosomal] ^ 2), expression(Brood~effect),
  expression(italic(h)[X - linked] ^ 2), expression(Other~environment),
  expression(italic(h)[Y - linked] ^ 2), expression(italic(epsilon))
)

posts_medians <- posts_draws %>%
  dplyr::select(-file) %>%
  group_by(x, y) %>%
  summarise(across(everything(), median), .groups = 'drop') %>%
  dplyr::select(x, y, h2:measurement_error) %>%
  pivot_longer(-c(x, y), names_to = 'source', values_to = 'estimate') %>%
  mutate(source_label = factor(source, facet_levels, facet_labels))

# set up the plot
P_car_pixel_h2 <- posts_medians %>%
  filter(source %in% c('h2', 'h2_auto', 'h2_X', 'h2_Y')) %>%
  ggplot(aes(x, y, fill = estimate)) +
  geom_raster(aes(alpha = alpha), bg_fish, fill = 'grey70', interpolate = TRUE) +
  geom_raster(interpolate = FALSE) +
  # scales
  scale_y_reverse() +
  scale_fill_viridis_c(limits = c(0, 1), name = 'Proportion of variance') +
  scale_alpha_identity() +
  # layout and tweaks
  facet_grid(rows = vars(source_label), labeller = label_parsed, switch = 'y') +
  coord_fixed() +
  labs(title = 'Orange') +
  theme(
    strip.text.y.left = element_text(angle = 0),
    legend.position = 'bottom',
    legend.title = element_text(vjust = 1),
    legend.key.width = grid::unit(25, "points"),
    legend.title.align = 1,
    legend.direction = 'horizontal',
    legend.box.just = 'right',
    axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), strip.background = element_blank(),
    plot.title = element_text(color = 'darkorange', hjust = 0.5),
    plot.margin = margin(0,0,0,0)
  )

