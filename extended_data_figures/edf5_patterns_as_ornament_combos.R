library(tidyverse)
library(patchwork)
library(magick)
library(ggtext)

source('selection_decisions/compile_decisions.R')
source('paper_figures/theme.R')

car_images <- list.files('data/carotenoid_coloration_warped', recursive = TRUE, full.names = TRUE) %>%
  setNames(., basename(.) |> tools::file_path_sans_ext())

mel_images <- list.files('data/melanic_coloration_warped_v2', recursive = TRUE, full.names = TRUE) %>%
  setNames(., basename(.) |> tools::file_path_sans_ext())

make_pattern_image <- function(pattern) {
  fish <- filter(present_10, .data$pattern == .env$pattern)$fish_id
  unique_ids <- filter(ornaments, fish_id %in% fish)$unique_id |> unique()

  orange_avg <- magick::image_read(car_images[unique_ids]) |> image_average()
  orange <- orange_avg |>
    image_channel('alpha') |> image_negate() |> image_threshold('black', '30%') |>
    image_composite(
      image_blank(
        width = image_info(orange_avg)$width, height = image_info(orange_avg)$height, color = 'orange'
      ),
      composite_image = _,
      'Multiply'
    ) |>
    image_transparent('black') |> image_colorize(100, 'orange')

  black <- magick::image_read(mel_images[unique_ids]) |> image_average() |>
    image_channel('alpha') |> image_threshold('black', '30%') |> image_threshold('white', '70%') |>
    image_transparent('white') |> image_colorize(100, 'black')

  bg1 <- image_read('data/extracted_fish_warped/replicate_1/gen_2/20201214_IMG_5421.png') %>%
    image_channel('alpha')
  bg2 <- image_negate(bg1)
  bg3 <- image_colorize(bg1, 100, 'grey80') %>% image_composite(bg2, operator = 'CopyOpacity')

  out <- image_flatten(c(bg3, orange, black))
  out_file <- tempfile(fileext = '.png')
  image_write(out, out_file)
  return(out_file)
}

# read in ornaments from all images
ornaments <- bind_rows(
  read_rds('ornament_analysis/car_ornaments.rds'),
  read_rds('ornament_analysis/mel_ornaments.rds')
) |>
  # join with fish id's, to get individual level info
  left_join(
    data.table::fread('photo_database.csv') |> select(unique_id, fish_id, facing_direction),
    join_by(unique_id)
  )

# reformat, make the pattern for each male, with their ornaments in wide-format
present_10 <- ornaments |>
  summarise(
    present_10 = mean(present_10),
    .by = c(ornament, fish_id, facing_direction)
  ) |>
  summarise(
    present_10 = round(mean(present_10)),
    .by = c(ornament, fish_id)
  ) |>
  pivot_wider(names_from = ornament, values_from = present_10) |>
  unite('pattern', -fish_id, sep = '', remove = FALSE)

# find the incidence of each ornament across the dataset
incidences <- present_10 |> summarise(across(c(-fish_id, -pattern), mean)) |> as_vector()

# calculate the expected frequency of each possible pattern, based on ornament incidence
expected_patterns <- present_10 |>
  select(-fish_id) |>
  distinct() |>
  #complete(across(-pattern)) |> # illegal for some reason?
  complete(
    car_1, car_2, car_3, car_4, car_5, car_6, car_7,
    mel_1, mel_2, mel_3, mel_4, mel_5, mel_6, mel_7, mel_8
  ) |>
  select(-pattern) |> unite('pattern', everything(), sep = '', remove = FALSE) |>
  mutate(
    predicted_frequency = (pick(car_1:mel_8)) |>
      apply(1, \(x) prod(incidences[x != 0]) * prod(1 - incidences[x == 0]))
  ) |>
  select(pattern, predicted_frequency)

# merge the expectations with the observations, fill in 0 values for unobserved patterns
observed_patterns <- present_10 |>
  count(pattern) |>
  mutate(observed_frequency = n / sum(n)) |>
  full_join(expected_patterns, join_by(pattern)) |>
  mutate(across(c(n, observed_frequency), \(x) coalesce(x, 0)))

# Based on these predictions, how many unique patterns would we expect?
expected_unique_patterns <- map_int(
  1:1000,
  \(i) {
    sample.int(
      n = nrow(expected_patterns),
      size = sum(observed_patterns$n),
      replace = TRUE,
      prob = observed_patterns$predicted_frequency
    ) |> n_distinct()
  }
)
obs_pred <- ggplot() +
  geom_histogram(aes(x = expected_unique_patterns, fill = 'Expected'), binwidth = 5) +
  geom_vline(aes(color = 'Observed', xintercept = sum(observed_patterns$n > 0))) +
  geom_pointrange(
    aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = 130),
    brms::posterior_summary(expected_unique_patterns) |> as.data.frame(),
    color = 'darkorange'
  ) +
  scale_fill_manual(values = "darkorange") +
  scale_color_manual(values = "black") +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(c(0, 0.01))) +
  labs(x = 'Unique ornament combinations', y = 'Count', color = NULL, fill = NULL) +
  coord_fixed(ratio = 6) +
  theme(legend.position = 'top')

# What are the most common patterns in our dataset?
p1_data <- slice_max(observed_patterns, observed_frequency, n = 30) |>
  mutate(
    file = map_chr(pattern, make_pattern_image),
    pattern2 = glue::glue("<img src='{file}' width='40'/>")
  )
most_common <- ggplot(p1_data, aes(observed_frequency, fct_reorder(pattern2, observed_frequency))) +
  geom_col(aes(fill = 'Observed'), width = 0.8) +
  #geom_col(aes(fill = 'Predicted'), width = 0.4) +
  geom_point(aes(predicted_frequency, color = 'Predicted'), shape = '|', size = 3) +
  scale_x_continuous(
    labels = scales::label_percent(),
    expand = expansion(c(0.01, 0.05)),
    sec.axis = sec_axis(\(x) x * sum(observed_patterns$n), name = 'Pattern count')
  ) +
  scale_fill_manual(values = c("black", 'grey60'), name = NULL) +
  scale_color_manual(values = "darkorange", name = NULL) +
  labs(x = 'Pattern frequency', y = NULL) +
  theme(axis.text.y = element_markdown(), legend.position = 'top')

edf5 <- most_common + inset_element(obs_pred, 0.5, 0.1, 1, 0.5)

ggsave('extended_data_figures/edf5.png', edf5, units = "cm", width = 18, height = 16, dpi = 300)
ggsave('extended_data_figures/edf5.pdf', edf5, units = "cm", width = 18, height = 16, device = quartz, type = 'pdf')
