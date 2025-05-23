library(tidyverse)
library(lme4)
library(patchwork)
source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('quant_gen/prepare_pedigrees.R')
source('paper_figures/theme.R')

# 1: The SNP data
pavs <- data.table::fread(
  'grep Texim_26586531 sequencing/gwas_sommer/texim_analyses/Texim.out'
) |>
  set_names(c('chr', 'pos', 'alt', 'presence', 'copy_id', 'file_name')) |>
  mutate(sample_name = file_name |>
           str_remove('/Users/jiachangfu/Desktop/01_short_reads/03_Pret_aln/bam/3.Pret_aln/') |>
           str_remove('_R.q20.sorted.rmdup.bam')
  ) |>
  inner_join(get_sampling_structure(), join_by(sample_name)) |>
  inner_join(get_phenotypes('pa_car_6'), join_by(fish_id, sample_name)) |>
  mutate(phenotype_value_f = factor(phenotype_value, 0:1, c('absent', 'present'))) |>
  add_patriline(ped_df)

# Make one that summarises the support
pavs2 <- pavs |>
  summarise(
    perc_support = mean(presence) * 100,
    copy_present = perc_support > 50,
    .by = c(sample_name, phenotype_value, phenotype_value_f, patriline)
  )

## Compare to phenotype O6 -------------------------------------------------------------------------
m <- pavs2 |>
  glm(phenotype_value ~ copy_present, data = _, family = 'binomial')

g2p_plot <- pavs2 |>
  ggplot(aes(copy_present, fill = phenotype_value_f)) +
  geom_bar(position = 'fill') +
  annotate(
    'text', x = 1.5, y = 0.8, vjust = -0.2, size = 2,
    label = str_glue('R^2 == {round(performance::r2_nagelkerke(m), 2)}'), parse = TRUE
  ) +
  scale_y_continuous(expand = expansion(c(0, 0)), name = 'Proportion w/ ornament', breaks = c(0, 0.5, 1)) +
  scale_x_discrete(limits = c('FALSE', 'TRUE'), labels = c('Absent', 'Present')) +
  scale_fill_manual(values = c(absent = 'white', present = 'orange3'), guide = 'none') +
  labs(x = 'Texim copy 6531', y = 'Proportion w/ ornament')

# 2: The expression data

expression <- data.table::fread('sequencing/gwas_sommer/texim_analyses/Real_Texim_exp.tsv') |>
  filter(copy_name == '6531', tissue == 'T', selection == 'up')

expr_plot <- ggplot(expression, aes(CPM, sex, group = interaction(tank, sex))) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, fill = 'black') +
  geom_text(
    aes(label = sprintf("%.2f", CPM)), position = position_dodge(width = 0.9),
    hjust = c(-0.2, 1.2, -0.2, 1.2, -0.2, 1.2),
    color = rep(c('black', 'white'), times = 3),
    size = 2
  ) +
  scale_x_continuous(expand = c(0, 0), name = 'Expression (CPM)') +
  scale_y_discrete(name = 'Sex', limits = c('F', 'M'), labels = c('Female', 'Male'))

edf10 <- (expr_plot | g2p_plot) + plot_annotation(tag_levels = 'A')

ggsave('extended_data_figures/edf10.png', edf10, units = "cm", width = 12, height = 6, dpi = 300)
ggsave('extended_data_figures/edf10.pdf', edf10, units = "cm", width = 12, height = 6, device = quartz, type = 'pdf')
