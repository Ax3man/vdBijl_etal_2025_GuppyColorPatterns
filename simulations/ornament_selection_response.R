# Understanding the response to selection of ornaments using the breeder's equation

# We need three things:
# 1) the heritability
# 2) the prevalence in the population
# 3) the prevalence among selected males
#
# There are two complications. Firstly, the presence of ornaments is binary, but the breeder's
# equation describes continuous traits. So we use the logit transformation to model the probability
# of the ornament being present; the so-called liability scale.
#
# Secondly, we could only directly select on males. So we need to adjust the effective selection.
# For autosomal loci, the effective effect of selection is only half as strong, but for Y-linked
# effects we do gain the full effect. Since in the breeder's equation the strength of selection and
# the heritability are linearly related, I choose to adjust the heritability instead of the
# effect of selection. Do note that this heritability is biased downwards, since we ignore selection
# effects for X-linked loci through males, as this effect skips one generation. We also selected
# indirectly on females through their breeding values, but because of the large uncertainty
# associated with those values, this was likely much less effective than the selection on males.

library(tidyverse)

predict_selection_responses <- function(ornament) {
  logit <- qlogis; inv.logit <- plogis

  # 1. Heritability of the presence of the ornament (on the liability scale)
  # We get this from the hurdle model presented in fig 4.
  m_df <- glue::glue('quant_gen/ornament_heritability/saved_models/{ornament}.rds') |>
    read_rds() |>
    tidybayes::gather_draws(
    sd_fish_id__hu_Intercept, sd_fish_idA__hu_Intercept, sd_fish_idX__hu_Intercept,
    sd_patriline__hu_Intercept, sd_grow_tank__hu_Intercept
  ) %>%
    ungroup() %>%
    mutate(
      # give clearer names
      .variable = fct_recode(
        factor(.variable),
        env = 'sd_fish_id__hu_Intercept', auto = 'sd_fish_idA__hu_Intercept',
        X = 'sd_fish_idX__hu_Intercept', Y = 'sd_patriline__hu_Intercept',
        tank = 'sd_grow_tank__hu_Intercept'
      ),
      # change to variances (from standard deviations)
      .value = .value ^ 2
    ) %>%
    pivot_wider(names_from = '.variable', values_from = '.value') %>%
    mutate(
      # encode sigma to be the logit link variance
      sigma = (pi ^ 2) / 3,
      Va = auto + X + Y,
      Ve = env + tank,
      Vtotal = Va + Ve + sigma,

      h2 = Va / Vtotal,
      h2_auto = auto / Vtotal,
      h2_X = X / Vtotal,
      h2_Y = Y / Vtotal,
    )
  h2_auto <- posterior::summarise_draws(posterior::rvar(m_df$h2_auto, dim = 1))$mean
  h2_X <- posterior::summarise_draws(posterior::rvar(m_df$h2_X, dim = 1))$mean
  h2_Y <- posterior::summarise_draws(posterior::rvar(m_df$h2_Y, dim = 1))$mean

  # 2 & 3: Find the prevalence of O6 for different groups of the population
  source('selection_decisions/compile_decisions.R')
  orn_per_image <- read_rds('ornament_analysis/car_ornaments.rds') |>
    filter(.data$ornament == .env$ornament) |> dplyr::select(unique_id, present_10)
  orn_per_fish <- data.table::fread('photo_database.csv') |>
    as_tibble() |>
    dplyr::select(fish_id, unique_id, facing_direction) |>
    full_join(orn_per_image, 'unique_id', relationship = "one-to-one") |>
    summarise(present_10 = mean(present_10, na.rm = TRUE), .by = c(fish_id, facing_direction)) |>
    summarise(present_10 = mean(present_10, na.rm = TRUE), .by = fish_id) |>
    mutate(fish_id = toupper(fish_id)) |>
    full_join(
      dplyr::select(selection, fish_id, replicate, generation, selection, decisions),
      by = 'fish_id', relationship = 'one-to-one'
    )

  population_incidences <- orn_per_fish |> summarise(
    pop_incidence = mean(present_10),
    .by = c(replicate, generation, selection)
  ) |>
    mutate(pop_incidence_liability = logit(pop_incidence))

  selected_incidences <- orn_per_fish |>
    filter(decisions %in% c('move to down', 'move to up')) |>
    summarise(
      sel_incidence = mean(present_10),
      .by = c(replicate, generation, decisions)
    ) |>
    mutate(sel_incidence_liability = logit(sel_incidence))

  result <- data.frame(current_gen = c('P', 'P', 'F1', 'F2', 'F3')) |>
    left_join(
      population_incidences |> dplyr::select(-pop_incidence, pop = pop_incidence_liability),
      join_by(current_gen == generation),
      relationship = 'many-to-many'
    ) |>
    mutate(
      selection = str_remove(selection, '_selected'),
      selection = ifelse(row_number() < 7, rep(c('up', 'down'), each = 3), selection)
    ) |>
    pivot_wider(
      id_cols = c(replicate, current_gen),
      names_from = selection, values_from = pop, names_prefix = 'pop_'
    ) |>
    mutate(delta_P = pop_up - pop_down) |>
    left_join(
      selected_incidences |> dplyr::select(-sel_incidence, selected = sel_incidence_liability),
      join_by(replicate, current_gen == generation)
    ) |>
    mutate(decisions = str_remove(decisions, 'move to ')) |>
    pivot_wider(names_from = decisions, values_from = selected, names_prefix = 'selected_') |>
    dplyr::select(-selected_NA) |>
    mutate(ornament = ornament)


  result |> mutate(
    # no x-linked effects for the first generation
    predicted_next_delta_P = ifelse(
      current_gen == 'P',
        (0.5 * h2_auto + h2_Y) * (selected_up - selected_down),
        (0.5 * h2_auto + h2_Y) * (selected_up - selected_down - delta_P) +
          ((1/3) * h2_X) * (lag(selected_up) - lag(selected_down) - lag(delta_P)) + delta_P
    ),
    true_next_delta_P = lead(delta_P),
    .by = replicate,
    .after = 'current_gen'
  )
}

out <- map_dfr(
  paste0('car_', 1:7),
  predict_selection_responses
)

r2 <- with(
  filter(out, !is.na(predicted_next_delta_P), !is.infinite(predicted_next_delta_P)),
  cor(predicted_next_delta_P, true_next_delta_P) ^ 2
)

# r2 = 0.84
cat('r^2 =', r2)

labs <- c('0.1', '1', '10', '100', '1,000'); br <- c(.1, 1, 10, 100, 1000)
out |>
  filter(is.finite(true_next_delta_P), is.finite(predicted_next_delta_P)) |>
  #ggplot(aes(exp(true_next_delta_P * sign(predicted_next_delta_P)), exp(abs(predicted_next_delta_P)))) +
  ggplot(aes(exp(predicted_next_delta_P), exp(true_next_delta_P))) +
  geom_hline(yintercept = 1, col = 'grey60', linewidth = 0.1) +
  geom_vline(xintercept = 1, col = 'grey60', linewidth = 0.1) +
  geom_abline(aes(linetype = '1:1', slope = 1, intercept = 0)) +
  #geom_smooth(method = 'lm', lty = 2, alpha = 0.1, linewidth = 0.5, color = 1) +
  #geom_path(aes(group = ornament), linewidth = 0.1) +
  geom_point(aes(color = ornament)) +
  scale_x_log10(labels = labs, breaks = br) + scale_y_log10(labels = labs, breaks = br) +
  scale_color_discrete(breaks = paste0('car_', 1:7), labels = paste0('O', 1:7)) +
  coord_fixed() +
  facet_grid(cols = vars(str_replace(replicate, '_', ' '))) +
  labs(
    y = 'Observed effect of selection\n(odds ratio)',
    x = 'Predicted effect of selection\n(odds ratio)',
    linetype = NULL,
    color = 'Ornament'
  ) +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() +
  theme(
    legend.position = 'top',
    strip.background = element_rect(color = NA, fill = 'grey90')
  )
ggsave('paper_figures_supplement/breeders_predictions.png', h = 4, w = 8)

