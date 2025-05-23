library(tidyverse)
source('paper_figures/theme.R')

model_results <- read_rds('visualization/per_pixel_models/mel_selection_model_results.rds')
contrasts <- map_dfr(model_results, 'contrast')
contrasts$generation_label <- with(contrasts, case_when(
  generation == 'F1' ~ 'F[1]', generation == 'F2' ~ 'F[2]', generation == 'F3' ~ 'F[3]',
))

complete_fish <- imager::load.image('data/extracted_fish_warped/replicate_3/gen_2/20210602_IMG_2097.png') %>%
  as.data.frame(wide = 'c') %>% filter(c.4 > 0.1) %>% dplyr::select(x, y, alpha = c.4)
eb <- element_blank()

# note that estimates are one the logg odds ratio (not response) scale
P_mel_selmap <- ggplot(contrasts, aes(x, y, fill = exp(-estimate))) +
  geom_raster(
    aes(x = x, y = y, alpha = alpha),
    fill = 'grey60', data = complete_fish, show.legend = FALSE, inherit.aes = FALSE
  ) +
  geom_raster() +
  scale_alpha_identity(guide = 'none') +
  scale_y_reverse() +
  scico::scale_fill_scico(
    palette = 'vik',
    trans = 'log',
    limits = c(1/31, 31),
    oob = scales::squish,
    breaks = c(1/30, 1/10, 1/3, 1, 3/1, 10, 30),
    labels = c('1/30', '1/10', '1/3', '1', '3', '10', '30'),
    guide = guide_colorbar(title.position = 'top', theme = theme(
      legend.key.height = grid::unit(100, "points"), legend.key.width = unit(0.4, 'lines'),
    ))
  ) +
  facet_grid(generation_label ~ ., switch = 'y', labeller = label_parsed) +
  coord_fixed(expand = FALSE) +
  labs(fill = 'Incidence\nodds ratio') +
  theme(
    legend.position = 'right',
    legend.title.align = 0.5,
    strip.background = eb, strip.text.y.left = element_text(angle = 0),
    axis.text = eb, axis.title = eb, axis.ticks = eb, axis.line = eb,
  )

source('selection_decisions/compile_decisions.R')

ped_df <- data.table::fread('data/pedigree.csv') %>%
  as_tibble() %>%
  mutate(
    sex = str_sub(animal, 1, 1),
    grow_tank = ifelse(is.na(sire), 'source_pop', paste(sire, dam, date_of_birth, sep = '_'))
  ) %>%
  mutate(across(animal:dam, tolower))

mel <- data.table::fread('photo_database.csv') %>%
  group_by(replicate, generation, fish_id, facing_direction) %>%
  summarise(mel_perc_v2 = mean(mel_perc_v2), .groups = 'drop_last') %>%
  summarise(mel_perc_v2 = mean(mel_perc_v2), .groups = 'drop') %>%
  mutate(
    replicate = factor(replicate, paste0('replicate_', 1:3), paste('Replicate', 1:3)),
    generation = factor(
      generation, c('parental_gen_1', 'gen_2', 'gen_3', 'gen_4'), c('P', 'F1', 'F2', 'F3'))
  )

pd <- mel %>%
  left_join(dplyr::select(selection, fish_id, selection) %>% mutate(fish_id = tolower(fish_id))) %>%
  left_join(ped_df, by = c('fish_id' = 'animal')) %>%
  left_join(
    dplyr::select(mel, fish_id, sire_mel_perc = mel_perc_v2, sire_generation = generation),
    by = c('sire' = 'fish_id')
  )

# Expect 300 missing values (P generation has no sires)
P_mel_pedigree <- ggplot(pd, aes(x = mel_perc_v2, y = generation, color = selection)) +
  geom_hline(yintercept = 1:4, col = 'black', linewidth = 0.2) +
  ggridges::geom_density_ridges(
    data = filter(pd, generation == 'P'),
    fill = 1, alpha = 0.5, color = 'black', panel_scaling = FALSE, scale = 1, linewidth = 0.25,
    quantile_lines = TRUE, quantiles = 0.5
  ) +
  ggridges::geom_density_ridges(
    aes(fill = selection, group = interaction(selection, generation)), filter(pd, generation != 'P'),
    alpha = 0.5, color = 'black', panel_scaling = FALSE, scale = .6, linewidth = 0.25,
    quantile_lines = TRUE, quantiles = 0.5
  ) +
  geom_segment(aes(xend = sire_mel_perc, yend = sire_generation), alpha = 0.1, linewidth = 0.15) +
  #geom_point(alpha = .6, size = 1, shape = '|') +
  #facet_wrap( ~ replicate) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(
    limits = rev(c('P', 'F1', 'F2', 'F3')),
    labels = rev(c('P', expression(F[1]), expression(F[2]), expression(F[3]))),
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_color_manual(
    values = c('navy', 'black', 'firebrick'),
    labels = c('Down-selected', 'Up-selected'),
    guide = 'none'
  ) +
  scale_fill_manual(
    values = c('navy', 'firebrick'),
    labels = c('Down-selected', 'Up-selected'),
    guide = guide_legend(override.aes = list(alpha = 1), )
  ) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = 'Black coloration (% body area)', y = 'Generation', fill = NULL) +
  theme(legend.position = 'top')

library(furrr)
library(uwot) # for UMAP
source('selection_decisions/compile_decisions.R')

photo_database <- data.table::fread('photo_database.csv') %>%
  select(replicate, generation, fish_id, facing_direction, unique_id) %>%
  as_tibble() %>%
  mutate(
    image_file = paste0('data/melanic_coloration_warped_v2/', replicate, '/', generation, '/', unique_id, '.png')
  )
##
embedding <- read_rds(
  'dimension_reduction/triplet_loss_encoders/embeddings/mel_model_ped_fullcolor_comparison/embed_dim_5.rds'
) %>%
  magrittr::set_colnames(paste0('V', 1:ncol(.))) %>%
  as.data.frame() %>%
  rownames_to_column('unique_id') %>%
  as_tibble()

umap_model <- umap(embedding[-1], ret_model = TRUE, fast_sgd = TRUE, spread = 2) #, n_neighbors = 50, spread = 5)

embeddings <- bind_cols(
  embedding,
  umap_model$embedding %>% as.data.frame() %>% setNames(paste0('UMAP', 1:2))
) %>%
  left_join(select(photo_database, fish_id, unique_id, facing_direction, image_file), 'unique_id')

emb_df <- embeddings %>%
  group_by(fish_id) %>%
  summarise(
    across(starts_with('UMAP'), mean),
    image_file = sample(image_file, 1)
  ) %>%
  left_join(selection %>% mutate(fish_id = tolower(fish_id)), 'fish_id')

embedding_selection <- ggplot(
  slice_sample(emb_df, prop = 1),
  aes(UMAP1, UMAP2, color = ifelse(generation == 'P', 'Stock', selection))
) +
  geom_point(alpha = 1, size = 0.3) +
  #stat_ellipse(linewidth = 0.2) +
  # geom_text(
  #   aes(x = 2, y = -15, label = paste('n =', n)), count(emb_df, generation),
  #   color = 'grey40', size = 3
  # ) +
  scale_color_manual(
    name = NULL,
    values = c(Stock = 'grey40', down_selected = 'navy', up_selected = 'firebrick'),
    breaks = c('Stock', 'down_selected', 'up_selected'),
    labels = c('Stock', 'Down-selected', 'Up-selected'),
    guide = guide_legend(override.aes = list(size = 0.5))
  ) +
  facet_wrap(vars(generation), nrow = 1) +
  coord_fixed() +
  #theme_classic() +
  labs(x = 'Pattern axis 1', y = 'Pattern axis 2') +
  theme(
    legend.position = 'top', strip.background = element_rect(color = NA, fill = 'grey90'),
    strip.text = element_text(size = 6), legend.box.margin = margin(0,0,0,0),
    legend.margin = margin(0,0,0,0), plot.margin = margin(0,0,0,0)
  )

edf2 <- wrap_elements(full = P_mel_pedigree) +
  (wrap_elements(full = P_mel_selmap) /
     wrap_elements(embedding_selection) +
     plot_layout(heights = c(1.6, 1))) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = 'A')


ggsave('extended_data_figures/edf2.png', edf2, units = "cm", width = 18, height = 10, dpi = 300)
ggsave('extended_data_figures/edf2.pdf', edf2, units = "cm", width = 18, height = 10, device = quartz, type = 'pdf')
