library(tidyverse)
library(brms)
library(rstan)
library(ggdist)
library(patchwork)
source('paper_figures/theme.R')

m <- read_rds('quant_gen/car_mel_heritability_model.rds')

vc <- VarCorr(m, summary = FALSE)

variances <- vc |> map(
  \(comp) bind_rows(
    data.frame(color = 'Orange', value = comp$cov[, 1, 1]),
    data.frame(color = 'Black', value = comp$cov[, 2, 2])
  )
) |>
  list_rbind(names_to = 'component') |>
  mutate(color = factor(color, c('Orange', 'Black')))

panelA <- ggplot(variances, aes(component, value, fill = color)) +
  stat_ccdfinterval(point_size = 1, interval_size_range = c(0.3, 0.6), .width = c(0.66, 0.95)) +
  scale_x_discrete(
    limits = c(
      'fish_idA', 'fish_idX', 'patriline', 'dam', 'grow_tank', 'fish_id',
      'facing_direction:fish_id', 'residual__'
    ),
    labels = expression(
      V[A]^italic(auto), V[A]^italic('X-linked'), V[A]^italic('Y-linked'),
      V[E]^italic(maternal), V[E]^italic(tank), V[E]^italic(other),
      V[asymm], sigma
    )
  ) +
  scale_fill_manual(values = c('Orange' = 'orange3', 'Black' = 'grey40'), guide = 'none') +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(vars(color), scales = 'free_y', nrow = 1) +
  labs(y = 'Variance contributed', x = NULL, tag = 'A') +
  theme(strip.background = element_rect(color = NA, fill = 'grey90'))

var_wide <- variances |>
  mutate(draw = row_number(), .by = c(color, component)) |>
  pivot_wider(names_from = component, values_from = value) |>
  mutate(
    # Note that we define h2 are the proportion of the *between individual* variance,
    # and so we do not include measurement error and asymmetry
    h2 = (fish_idA + fish_idX + patriline) / (fish_idA + fish_idX + patriline + dam + grow_tank + fish_id),
    h2_x = fish_idX / (fish_idA + fish_idX + patriline + dam + grow_tank + fish_id),
    h2_y = patriline / (fish_idA + fish_idX + patriline + dam + grow_tank + fish_id),
  )

h2 <- ggplot(var_wide, aes(h2, color, fill = color)) +
  stat_halfeye(point_size = 1, interval_size_range = c(0.3, 0.6), .width = c(0.66, 0.95)) +
  scale_fill_manual(values = c('Orange' = 'orange3', 'Black' = 'grey40'), guide = 'none') +
  xlim(0, 1) +
  scale_y_discrete(expand = expansion(add = c(0.1, 0)), name = NULL) +
  labs(x = expression(paste('Heritability (', italic(h^2), ')')))

h2x <- h2 + aes(h2_x) + labs(x = expression(italic(h)[italic('X-linked')]^2))
h2y <- h2 + aes(h2_y) + labs(x = expression(italic(h)[italic('Y-linked')]^2))

panelB <- (h2 + labs(tag = 'B') | h2x | h2y) + plot_layout(axes = 'collect')

correlations <- map(vc, \(comp) data.frame(value = comp$cor[, 1, 2])) |>
  list_rbind(names_to = 'component')

panelC <- ggplot(correlations, aes(value, component)) +
  geom_vline(xintercept = 0, color = 'grey40', lty = 2) +
  stat_halfeye(
    normalize = 'groups', point_size = 1, interval_size_range = c(0.3, 0.6), .width = c(0.66, 0.95)
  ) +
  scale_y_discrete(
    limits = rev(c(
      'fish_idA',
      'fish_idX',
      'patriline',
      'dam',
      'grow_tank',
      'fish_id',
      'facing_direction:fish_id',
      'residual__'
    )),
    labels = rev(c(
      expression(r[A]^auto), expression(r[A]^X), expression(r[A]^Y),
      expression(r[E]^dam), expression(r[E]^tank), expression(r[E]^other), expression(r[asymm]),
      expression(r[sigma])
    )),
    expand = c(0.02, 0)
  ) +
  theme(axis.text.y = element_text(hjust = 0)) +
  labs(x = 'Correlation', y = NULL, tag = 'C')

edf1 <- panelA / panelB / (panelC | plot_spacer()) + plot_layout(heights = c(1.5, 1, 2))

ggsave('extended_data_figures/edf1.png', edf1, units = "cm", width = 18, height = 12, dpi = 300)
ggsave('extended_data_figures/edf1.pdf', edf1, units = "cm", width = 18, height = 12, device = quartz, type = 'pdf')
