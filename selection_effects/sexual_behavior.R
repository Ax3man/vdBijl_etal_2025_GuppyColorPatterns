library(tidyverse)
library(emmeans)
library(patchwork)
library(glmmTMB)

d <- data.table::fread('data/behavior_morphology.csv')

d_filtered <- filter(
  d,
  !receptive_female,       # disregard trials with receptive females
  following > 30,          # disregard trials with unreceptive males
  male_activity < 4,       # disregard trials with stressed males
  female_activity < 4      # disregard trials with stressed females
)
count(d_filtered, selection)

# total displays
m_display <- glmmTMB::glmmTMB(
  display ~ 1 + selection + standard_length + female_standard_length + (1 | selection:replicate),
  data = d_filtered, family = 'nbinom1'
)
summary(m_display)
emmeans(m_display, pairwise ~ selection, type = 'response')

m_short_display <- update(m_display, short_display ~ .)
summary(m_short_display)
emmeans(m_short_display, pairwise ~ selection, type = 'response')

m_long_display <- update(m_display, long_display ~ .)
summary(m_long_display)
emmeans(m_long_display, pairwise ~ selection, type = 'response')

m_display_duration <- update(m_display, display_duration ~ .)
summary(m_display_duration)
emmeans(m_display_duration, pairwise ~ selection, type = 'response')

m_following <- update(m_display, following ~ ., family = gaussian('log'))
summary(m_following)
emmeans(m_following, pairwise ~ selection, type = 'response')

m_sneak <- update(m_display, sneak ~ .)
summary(m_sneak)
emmeans(m_sneak, pairwise ~ selection, type = 'response')

make_plot <- function(y, lab, model) {
  ggplot(d_filtered, aes(selection, .data[[y]])) +
    geom_violin(col = 1, linewidth = 0.1, fill = 'grey80') +
    geom_pointrange(
      aes(y = response, ymin = asymp.LCL, ymax = asymp.UCL),
      emmeans(model, ~ selection, type = 'response') %>% as.data.frame()
    ) +
    scale_x_discrete(limits = c('down-selected', 'up-selected')) +
    theme_classic() +
    ylab(lab)
}

wrap_plots(
  #make_plot('display', 'Total displays', m_display),
  make_plot('short_display', 'Short displays', m_short_display),
  make_plot('long_display', 'Long displays', m_long_display),
  make_plot('display_duration', 'Display duration', m_display_duration)
) + plot_annotation(tag_levels = 'A')
ggsave('paper_figures_supplement/display_behavior.png', w = 10, h = 4.5)

wrap_plots(
  make_plot('sneak', 'Sneak copulation attempts', m_sneak),
  ggplot(d_filtered, aes(selection, following)) +
    geom_violin(col = 1, linewidth = 0.1, fill = 'grey80') +
    geom_pointrange(
      aes(y = response, ymin = lower.CL, ymax = upper.CL),
      emmeans(m_following, ~ selection, type = 'response') %>% as.data.frame()
    ) +
    scale_x_discrete(limits = c('down-selected', 'up-selected')) +
    theme_classic() +
    ylab('Time spent following')
) + plot_annotation(tag_levels = 'A')
ggsave('paper_figures_supplement/sneak_behavior.png', w = 6.6, h = 4.5)
