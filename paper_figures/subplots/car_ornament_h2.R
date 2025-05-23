library(tidyverse)
library(ggdist)

source('paper_figures/theme.R')
eb <- element_blank()

if (FALSE) { # we don't need to rerun this part when perfecting figures
  library(tidybayes)
  library(furrr)

  plan(multisession, workers = 7)

  models <- list.files(
    'quant_gen/ornament_heritability/saved_models',
    pattern = 'car',
    full.names = TRUE
  ) %>%
    setNames(., basename(.) %>% tools::file_path_sans_ext())

  all_draws <- future_map_dfr(
    models,
    \(model) {
      model %>%
        read_rds() %>%
        tidybayes::gather_draws(
          sd_fish_id__Intercept, sd_fish_idA__Intercept, sd_fish_idX__Intercept,
          sd_patriline__Intercept, sd_grow_tank__Intercept,

          sd_fish_id__hu_Intercept, sd_fish_idA__hu_Intercept, sd_fish_idX__hu_Intercept,
          sd_patriline__hu_Intercept, sd_grow_tank__hu_Intercept,

          sigma
        ) %>%
        ungroup() %>%
        mutate(
          # give clearer names
          .variable = fct_recode(
            factor(.variable),
            size_env = 'sd_fish_id__Intercept', size_auto = 'sd_fish_idA__Intercept',
            size_X = 'sd_fish_idX__Intercept', size_Y = 'sd_patriline__Intercept',
            size_tank = 'sd_grow_tank__Intercept',

            ap_env = 'sd_fish_id__hu_Intercept', ap_auto = 'sd_fish_idA__hu_Intercept',
            ap_X = 'sd_fish_idX__hu_Intercept', ap_Y = 'sd_patriline__hu_Intercept',
            ap_tank = 'sd_grow_tank__hu_Intercept',

            size_sigma = 'sigma'
          ),
          # change to variances (from standard deviations)
          .value = .value ^ 2
        ) %>%
        separate(.variable, into = c('.response', '.variable'), sep = '_') %>%
        pivot_wider(names_from = '.variable', values_from = '.value') %>%
        mutate(
          .response = fct_recode(factor(.response), 'Absence/Presence' = 'ap', 'Size' = 'size'),

          # encode sigma to be the logit link variance, but only for the presence/absence part
          sigma = coalesce(sigma, (pi ^ 2) / 3),

          Va = auto + X + Y,
          Ve = env + tank,
          Vtotal = Va + Ve + sigma,

          h2 = Va / Vtotal,
          h2_auto = auto / Vtotal,
          h2_X = X / Vtotal,
          h2_Y = Y / Vtotal,
        )
    },
    .id = 'ornament', .options = furrr_options(seed = TRUE)
  )

  plan(sequential)

  variances <- all_draws %>%
    select(ornament, .response, h2:h2_Y) %>%
    pivot_longer(cols = h2:h2_Y)

  write_rds(variances, 'paper_figures/subplots/car_ornament_h2_data.rds')
}

variances <- read_rds('paper_figures/subplots/car_ornament_h2_data.rds')

P_car_orn_h2 <- variances %>%
  ggplot(aes(value, name, fill = fct_rev(.response))) +
  geom_hline(yintercept = 3.5, color = 'grey', linewidth = 0.1) +
  stat_ccdfinterval(
    aes(slab_alpha = after_stat(f)),
    thickness = 1, fill_type = "segments", geom = 'slab', #fill = 'grey20',
    position = position_dodge(width = 0.7), height = 0.75
  ) +
  ggdist::stat_pointinterval(
    linewidth = 0.1, point_size = 0.1, position = position_dodge(width = 0.7), show.legend = FALSE
  ) +
  scale_x_continuous(breaks = c(.25, .5, .75), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_discrete(
    breaks = c('h2', 'h2_auto', 'h2_X', 'h2_Y'), limits = rev(c('h2', 'h2_auto', 'h2_X', 'h2_Y')),
    labels = (expression(italic(h)^2, italic(h)[auto]^2, italic(h)[X]^2, italic(h)[Y]^2))
  ) +
  scale_fill_manual(
    name = NULL,
    limits = c('Absence/Presence', 'Size'),
    values = c(Size = scales::muted('#c4ddc7'), 'Absence/Presence' = scales::muted('#afc8e1')),
    labels = c(Size = 'Size', 'Absence/Presence' = 'Incidence')
  ) +
  scale_slab_alpha_continuous(guide = 'none') +
  labs(y = NULL, x = NULL) +
  facet_grid(~ornament) +
  theme(
    strip.text = eb, strip.background = eb
  )
