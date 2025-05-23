library(tidyverse)

source('quant_gen/prepare_pedigrees.R')

# Get the carotenoid values for each photo
d <- data.table::fread('photo_database.csv') %>%
  as_tibble() %>%
  select(replicate, generation, fish_id, facing_direction, car_perc, mel_perc_v2) %>%
  mutate(fish_idA = fish_id, fish_idX = fish_id) %>%
  # Add the patriline information to the dataset, as well as dam and sire id's
  add_patriline(ped_df) %>%
  left_join(ped_df, by = c('fish_id' = 'animal'))

## ACTUAL FITTING ----------------------------------------------------------------------------------
library(brms)
m <- brm(
  bf(
    # Bivariate outcome
    mvbind(car_perc, mel_perc_v2) ~
      # intercept
      1 +
      # genetic variances (autosomal, x-linked, y-linked) with their correlations
      (1 | Acor | gr(fish_idA, cov = A)) + (1 | Xcor | gr(fish_idX, cov = X)) + (1 | Ycor | patriline) +
      # environmental variances (maternal effects, tank effects, and other) with their correlations
      (1 | Mcor | dam) + (1 | Tcor | grow_tank) + (1 | Ecor | fish_id) +
      # variance due to asymmetry between the sides with the correlation
      (1 | AScor | facing_direction:fish_id)
    # Set the residual correlation, which is the correlation in measurement error
    ) + set_rescor(TRUE),
  data = d, data2 = list(A = pedA, X = pedX),
  backend = 'cmdstanr', chains = 8, cores = 8, warmup = 1500, iter = 2500, threads = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file = 'quant_gen/car_mel_heritability_model.rds', file_refit = 'on_change'
)

# Mean chain execution time: 8.1 days
# Total execution time: 8.7 days

if (FALSE) { # run interactively

  m <- read_rds('quant_gen/car_mel_heritability_model.rds')

  summary(m)
  cowplot::plot_grid(
    pp_check(m, nd = 10, resp = 'carperc'),
    pp_check(m, nd = 10, resp = 'melpercv2'), ncol = 1
  )

  # define function to make the plots
  make_var_plots <- function(var) {
    var_wide <- var %>%
      pivot_wider(names_from = .variable, values_from = .value) %>%
      mutate(
        # Note that we define h2 are the proportion of the *between individual* variance,
        # and so we do not include measurement error and asymmetry (they are small anyway)
        h2 = (auto + X + Y) / (auto + X + Y + maternal + tank + env),
        h2_auto = auto / (auto + X + Y + maternal + tank + env),
        h2_x = X / (auto + X + Y + maternal + tank + env),
        h2_y = Y / (auto + X + Y + maternal + tank + env),
      )

    # heritability in numbers:
    brms::posterior_summary(var_wide$h2)
    brms::posterior_summary(var_wide$h2_auto)
    brms::posterior_summary(var_wide$h2_x)
    brms::posterior_summary(var_wide$h2_y)

    p1 <- ggplot(var, aes(.variable, .value)) +
      stat_ccdfinterval(slab_type = 'ccdf') +
      scale_x_discrete(
        limits = c('auto', 'X', 'Y', 'maternal', 'tank', 'env', 'asymmetry', 'sigma'),
        labels = expression(
          V[A]^italic(auto), V[A]^italic('X-linked'), V[A]^italic('Y-linked'),
          V[E]^italic(maternal), V[E]^italic(tank), V[E]^italic(other),
          V[asymm], atop('', sigma)
        )
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      #coord_cartesian(ylim = c(0, 65)) +
      labs(y = 'Variance contributed', x = NULL) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 11, color = 1))
    h2_2 <- ggplot(var_wide, aes(h2)) +
      stat_halfeye() +
      xlim(0, 1) +
      labs(x = expression(paste('Heritability (', italic(h^2), ')')), y = 'Posterior density') +
      theme_classic()
    X_2 <- ggplot(var_wide, aes(h2_x)) +
      stat_halfeye() +
      xlim(0, 1) +
      labs(x = expression(italic(h)[italic('X-linked')]^2), y = NULL) + #, tag = 'D') +
      theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
    Y_2 <- ggplot(var_wide, aes(h2_y)) +
      stat_halfeye() +
      xlim(0, 1) +
      labs(x = expression(italic(h)[italic('Y-linked')]^2), y = NULL) + #, tag = 'D') +
      theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

    patchwork::wrap_plots(
      p1,
      patchwork::wrap_plots(h2_2, X_2, Y_2, nrow = 1, widths = c(1.2, 1, 1)),
      ncol = 1, heights = c(3, 2)
    )
  }

  ##########################################
  ### Visualize the carotenoid variances ###
  ##########################################

  library(tidybayes)

  var_car <- gather_draws(
    m,
    sd_fish_id__carperc_Intercept,
    sd_fish_idA__carperc_Intercept,
    sd_fish_idX__carperc_Intercept,
    sd_patriline__carperc_Intercept,
    sd_dam__carperc_Intercept,
    sd_grow_tank__carperc_Intercept,
    `sd_facing_direction:fish_id__carperc_Intercept`,
    sigma_carperc
  ) %>%
    ungroup() %>%
    mutate(
      # give clearer names
      .variable = fct_recode(
        factor(.variable),
        env = 'sd_fish_id__carperc_Intercept',
        auto = 'sd_fish_idA__carperc_Intercept',
        X = 'sd_fish_idX__carperc_Intercept',
        Y = 'sd_patriline__carperc_Intercept',
        maternal = 'sd_dam__carperc_Intercept',
        tank = 'sd_grow_tank__carperc_Intercept',
        asymmetry = 'sd_facing_direction:fish_id__carperc_Intercept',
        sigma = 'sigma_carperc'
      ),
      # change to variances (from standard deviations)
      .value = .value ^ 2
    )

  make_var_plots(var_car)
  ggsave(
    'paper_figures_supplement/car_var_components_full.png',
    dpi = 600, w = 7.5, height = 5, bg = 'white'
  )

  ##########################################
  ### Visualize the carotenoid variances ###
  ##########################################

  var_mel <- gather_draws(
    m,
    sd_fish_id__melpercv2_Intercept,
    sd_fish_idA__melpercv2_Intercept,
    sd_fish_idX__melpercv2_Intercept,
    sd_patriline__melpercv2_Intercept,
    sd_dam__melpercv2_Intercept,
    sd_grow_tank__melpercv2_Intercept,
    `sd_facing_direction:fish_id__melpercv2_Intercept`,
    sigma_melpercv2
  ) %>%
    ungroup() %>%
    mutate(
      # give clearer names
      .variable = fct_recode(
        factor(.variable),
        env = 'sd_fish_id__melpercv2_Intercept',
        auto = 'sd_fish_idA__melpercv2_Intercept',
        X = 'sd_fish_idX__melpercv2_Intercept',
        Y = 'sd_patriline__melpercv2_Intercept',
        maternal = 'sd_dam__melpercv2_Intercept',
        tank = 'sd_grow_tank__melpercv2_Intercept',
        asymmetry = 'sd_facing_direction:fish_id__melpercv2_Intercept',
        sigma = 'sigma_melpercv2'
      ),
      # change to variances (from standard deviations)
      .value = .value ^ 2
    )

  make_var_plots(var_mel)
  ggsave(
    'paper_figures_supplement/mel_var_components_full.png',
    dpi = 600, w = 7.5, height = 5, bg = 'white'
  )

  ##################################
  ### Visualize the correlations ###
  ##################################
  cor_draws <- gather_draws(
    m,
    cor_fish_idA__carperc_Intercept__melpercv2_Intercept,
    cor_fish_idX__carperc_Intercept__melpercv2_Intercept,
    cor_patriline__carperc_Intercept__melpercv2_Intercept,
    cor_dam__carperc_Intercept__melpercv2_Intercept,
    cor_grow_tank__carperc_Intercept__melpercv2_Intercept,
    cor_fish_id__carperc_Intercept__melpercv2_Intercept,
    `cor_facing_direction:fish_id__carperc_Intercept__melpercv2_Intercept`,
    rescor__carperc__melpercv2
  )

  split(cor_draws$.value, cor_draws$.variable) %>%
    map_dfr(\(x) brms::posterior_summary(x) %>% as.data.frame(), .id = 'var')

  ggplot(cor_draws, aes(.value, .variable)) +
    geom_vline(xintercept = 0, color = 'grey40', lty = 2) +
    stat_halfeye(normalize = 'all') +
    scale_y_discrete(
      limits = rev(c(
        'cor_fish_idA__carperc_Intercept__melpercv2_Intercept',
        'cor_fish_idX__carperc_Intercept__melpercv2_Intercept',
        'cor_patriline__carperc_Intercept__melpercv2_Intercept',
        'cor_dam__carperc_Intercept__melpercv2_Intercept',
        'cor_grow_tank__carperc_Intercept__melpercv2_Intercept',
        'cor_fish_id__carperc_Intercept__melpercv2_Intercept',
        'cor_facing_direction:fish_id__carperc_Intercept__melpercv2_Intercept',
        'rescor__carperc__melpercv2'
      )),
      labels = rev(c(
        expression(r[A]^auto), expression(r[A]^X), expression(r[A]^Y),
        expression(r[E]^dam), expression(r[E]^tank), expression(r[E]^other), expression(r[asymm]),
        expression(r[sigma])
      )),
      expand = c(0.02, 0)
    ) +
    theme_classic() + theme(axis.text.y = element_text(color = 1, size = 11, hjust = 0)) +
    labs(x = 'Correlation', y = NULL)
  ggsave('paper_figures_supplement/correlation_posteriors.png', w = 4, h = 4, dpi = 600, bg = 'white')

  # d %>%
  #   group_by(fish_id) %>%
  #   summarise(across(c(car_perc, mel_perc_v2), mean)) %>%
  #   ggplot(aes(car_perc, mel_perc_v2)) +
  #   geom_point(alpha = 0.6, shape = 19, stroke = 0, size = 2) +
  #   geom_smooth(method = 'lm', color = 'firebrick') +
  #   theme_minimal() +
  #   labs(x = '% carotenoid body coloration', y = '% melanic body coloration',
  #        subtitle = 'Each point is a male')
  # ggsave('visualization/basic_plots/phenotypic_correlation.png', dpi = 600, w = 5, height = 4, bg = 'white')

}
