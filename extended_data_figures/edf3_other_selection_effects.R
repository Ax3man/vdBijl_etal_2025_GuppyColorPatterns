library(tidyverse)
library(glmmTMB)
library(brms)
library(emmeans)
library(patchwork)
source('paper_figures/theme.R')

## Behavior
dbeh <- data.table::fread('data/behavior_morphology.csv')

dbeh_filtered <- filter(
  dbeh,
  !receptive_female,       # disregard trials with receptive females
  following > 30,          # disregard trials with unreceptive males
  male_activity < 4,       # disregard trials with stressed males
  female_activity < 4      # disregard trials with stressed females
)

# total displays
m_display <- glmmTMB(
  display ~ 1 + selection + standard_length + female_standard_length + (1 | selection:replicate),
  data = dbeh_filtered, family = 'nbinom1'
)
m_short_display <- update(m_display, short_display ~ .)
m_long_display <- update(m_display, long_display ~ .)
m_display_duration <- update(m_display, display_duration ~ .)
m_following <- update(m_display, following ~ ., family = gaussian('log'))
m_sneak <- update(m_display, sneak ~ .)

make_plot <- function(y, lab, model) {
  ggplot(dbeh_filtered, aes(selection, .data[[y]])) +
    geom_violin(color = 1, fill = 'grey80') +
    geom_pointrange(
      aes(y = response, ymin = asymp.LCL, ymax = asymp.UCL),
      emmeans(model, ~ selection, type = 'response') %>% as.data.frame()
    ) +
    scale_x_discrete(limits = c('down-selected', 'up-selected'), labels = c('Down', 'Up')) +
    labs(x = 'Selection', y = lab)
}

behavior <- wrap_plots(
  #make_plot('display', 'Total displays', m_display),
  make_plot('short_display', 'Short displays (count)', m_short_display),
  make_plot('long_display', 'Long displays (count)', m_long_display),
  make_plot('display_duration', 'Display duration (s)', m_display_duration),
  make_plot('sneak', 'Sneak copulation attempts (count)', m_sneak),
  ggplot(dbeh_filtered, aes(selection, following)) +
    geom_violin(color = 1, fill = 'grey80') +
    geom_pointrange(
      aes(y = response, ymin = lower.CL, ymax = upper.CL),
      emmeans(m_following, ~ selection, type = 'response') %>% as.data.frame()
    ) +
    scale_x_discrete(limits = c('down-selected', 'up-selected'), labels = c('Down', 'Up')) +
    labs(x = 'Selection', y = 'Time spent following (s)')
)

## Gross morphology
morph <- data.table::fread('data/behavior_morphology.csv')

m_SL <- brm(
  standard_length ~ 1 + selection + (1 | selection:replicate),
  data = morph, family = 'normal',
  cores = 4, control = list(adapt_delta = 0.99)
)
m_tail_length <- update(m_SL, formula = tail_length ~ ., newdata = morph, cores = 4)
m_tail_area <- update(m_SL, formula = tail_area ~ ., newdata = morph, cores = 4)

make_plot <- function(y, lab, model) {
  ggplot(morph, aes(selection, .data[[y]])) +
    geom_violin(color = 1, fill = 'grey80') +
    geom_pointrange(
      aes(y = emmean, ymin = lower.HPD, ymax = upper.HPD),
      emmeans(model, ~ selection) %>% as.data.frame()
    ) +
    scale_x_discrete(limits = c('down-selected', 'up-selected'), labels = c('Down', 'Up')) +
    labs(x = 'Selection', y = lab)
}

morphology <- (
  make_plot('standard_length', 'Standard length (cm)', m_SL) |
    make_plot('tail_length', 'Tail length (cm)', m_tail_length) |
    make_plot('tail_area', expression(Tail~area~(cm^2)), m_tail_area)
) + plot_annotation(tag_levels = 'A')

## Life history

source('selection_decisions/compile_decisions.R')
breeding <- data.table::fread('data/breeding_data.csv') %>%
  mutate(date = lubridate::as_date(date, format = '%Y/%m/%d'))

dat <- breeding %>%
  inner_join(dplyr::select(selection, fish_id, replicate, generation, selection), c('sire' = 'fish_id')) %>%
  group_by(replicate, generation) %>%
  mutate(
    days_since_first_brood = as.numeric(date - min(date, na.rm = TRUE)),
    generation = factor(generation, c('P', 'F1', 'F2'))
  )

first_broods <- dat %>%
  # drop P generation, since they weren't selected on yet, note that F3 never reproduced.
  filter(generation != 'P') %>%
  group_by(dam) %>%
  slice_min(days_since_first_brood) %>%
  ungroup()

m_time_to_first_brood <- brm(
  days_since_first_brood ~ 1 + selection * generation + (1 | replicate:generation),
  data = first_broods,
  cores = 4, control = list(adapt_delta = 0.999), iter = 4000
)
P_brood_timing <- ggplot(
  first_broods,
  aes(selection, days_since_first_brood, group = interaction(selection, generation))
) +
  geom_violin(color = 1, fill = 'grey80') +
  geom_pointrange(
    aes(y = emmean, ymin = lower.HPD, ymax = upper.HPD),
    as.data.frame(emmeans(m_time_to_first_brood, ~ selection | generation))
  ) +
  facet_grid(
    cols = vars(factor(generation, labels = c(expression(F[1]), expression(F[2])))),
    labeller = label_parsed
  ) +
  scale_x_discrete(limits = c('down_selected', 'up_selected'), labels = c('Down', 'Up')) +
  labs(
    y = 'Time to first brood\n(days from first brood of the generation)',
    x = 'Selection'
  ) +
  theme(legend.position = 'top', strip.background = element_rect(color = NA, fill = 'grey90'))

m_fecundity <- brm(
  nr_offspring ~ selection * generation + (1 | replicate:generation),
  data = first_broods, family = negbinomial(),
  cores = 8, chains = 8, iter = 10000, control = list(adapt_delta = 0.999, max_treedepth = 12)
)

P_fecundity <- ggplot(first_broods, aes(selection, nr_offspring, group = interaction(selection, generation))) +
  geom_violin(color = 1, fill = 'grey80') +
  geom_pointrange(
    aes(y = prob, ymin = lower.HPD, ymax = upper.HPD),
    as.data.frame(emmeans(m_fecundity, ~ selection | generation, type = 'response'))
  ) +
  facet_grid(
    cols = vars(factor(generation, labels = c(expression(F[1]), expression(F[2])))),
    labeller = label_parsed
  ) +
  scale_x_discrete(limits = c('down_selected', 'up_selected'), labels = c('Down', 'Up')) +
  labs(
    y = 'Fecundity\n(nr of offspring in first brood)',
    x = 'Selection'
  ) +
  theme(legend.position = 'top', strip.background = element_rect(color = NA, fill = 'grey90'))

life_history <- (P_fecundity | P_brood_timing)

# Combine and save

edf3 <- (behavior + plot_layout(nrow = 1)) /
  (morphology + plot_layout(nrow = 1)) /
  (life_history) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = margin(3, 3, 3, 3))

ggsave('extended_data_figures/edf3.png', edf3, units = "cm", width = 18, height = 16, dpi = 300)
ggsave('extended_data_figures/edf3.pdf', edf3, units = "cm", width = 18, height = 16, device = quartz, type = 'pdf')

