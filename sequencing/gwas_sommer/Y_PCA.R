library(irlba)
library(tidyverse)

source('sequencing/genomics_helpers.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')

# Import the VCF file into a GDS file
vcf_file <- "sequencing/gwas/filtered.vcf.gz"

g <- data.table::fread('sequencing/gwas_sommer/snp_inheritance/sommer_snp_inheritance/female_car_PIE.csv.gz') |>
  group_by(chr, start, end, alt) |>
  slice_max(AIC_weight) %>%
  filter(AIC_weight > 0.8, source == 'Y') |>
  ungroup()

if (FALSE) { # only run interactively, otherwise load from disk
  geno <- load_vcf(region = str_glue_data(g, '{chr}:{start}-{end}'))

  pca <- prcomp_irlba(geno$mat, n = 20, center = TRUE, scale. = FALSE)
  rownames(pca$x) <- rownames(geno$mat)
  write_rds(pca, 'sequencing/gwas_sommer/Y_haplotype_PCA.rds', compress = 'gz')
}
pca <- read_rds('sequencing/gwas_sommer/Y_haplotype_PCA.rds')
summary(pca)

# scree plot
scree <- summary(pca)$importance %>%
  as.data.frame() %>%
  dplyr::slice(2) %>%
  pivot_longer(everything()) %>%
  mutate(name = parse_number(name)) %>%
  ggplot(aes(name, value)) + geom_line() + geom_point() + theme_classic() +
  labs(title = 'Scree plot', x = 'PC', y = 'Proportion of variance')

# hclust
dist_matrix <- dist(pca$x) # distance matrix from the first two PCs
hc <- hclust(dist_matrix)
#plot(hc) # to visualize the dendrogram

##
source('quant_gen/prepare_pedigrees.R')
source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('selection_decisions/compile_decisions.R')
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
orn <- read_rds('ornament_analysis/car_ornaments.rds') %>%
  left_join(data.table::fread('photo_database.csv') %>% dplyr::select(unique_id, fish_id)) %>%
  group_by(fish_id, ornament) %>%
  summarise(present = Mode(present)) %>%
  mutate(ornament = str_replace(ornament, 'car_', 'O')) %>%
  pivot_wider(names_from = ornament, values_from = present)

ind_df <- as.data.frame(pca$x) %>% rownames_to_column('sample_name') %>% as_tibble() %>%
  left_join(get_sampling_structure(), join_by(sample_name)) %>%
  left_join(mutate(selection, fish_id = tolower(fish_id)), join_by(fish_id)) %>%
  left_join(orn, join_by(fish_id)) %>%
  add_patriline(ped_df)

ind_df$hclust_groups <- factor(cutree(hc, k = 4)) |> fct_infreq() |> factor(labels = 1:4)

ind_df %>%
  mutate(Y_haplogroup = factor(hclust_groups, levels = 1:4, labels = paste0('Y', 1:4))) %>%
  dplyr::select(sample_name, fish_id, Y_haplogroup) %>%
  data.table::fwrite('sequencing/Y_haplogroups.csv')

ind_df %>%
  pivot_longer(O1:O7, names_to = 'ornament', values_to = 'presence') %>%
  ggplot(aes(hclust_groups, fill = factor(presence, 0:1, c('absent', 'present')))) +
  geom_bar() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c('grey20', 'darkorange'), name = NULL) +
  facet_grid(cols = vars(ornament)) +
  theme_classic() +
  labs(x = 'Assigned Y-haplogroup')

biplot <- ggplot(ind_df, aes(PC1, PC2, color = hclust_groups)) +
  geom_point() + theme_classic() + coord_fixed() +
  scale_color_brewer(palette = "Set1", labels = paste0('Y', 1:4)) +
  labs(
    x = 'PC1 (29%)', y = 'PC2 (23%)',
    color = 'Assigned\nY-haplogroup',
    #subtitle = 'PCA on SNPs associated w/ orange pattern and Y-like inheritance'
  ) +
  theme(legend.position = 'top')

biplot

biplot + geom_point(aes(color = patriline), ind_df %>% group_by(patriline) %>% dplyr::filter(n() > 8))
biplot %+% (ind_df %>% group_by(patriline) %>% filter(n() > 8)) + aes(color = patriline) + stat_ellipse()

biplot + aes(color = selection) + facet_grid(cols = vars(replicate))
biplot + aes(color = factor(O6)) + scale_color_manual(values = c('black', 'darkorange'))

cowplot::plot_grid(biplot, scree, ncol = 2)

ind_df %>%
  pivot_longer(O1:O7, names_to = 'ornament', values_to = 'presence') %>%
  ggplot(aes(hclust_groups, fill = factor(presence, 0:1, c('absent', 'present')))) +
  geom_bar(position = 'fill') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c('grey20', 'darkorange'), name = NULL) +
  facet_grid(cols = vars(ornament)) +
  theme_classic() +
  labs(x = 'Assigned Y-haplotype', y = 'proportion')

ggplot(ind_df, aes(hclust_groups, car_perc)) +
  geom_violin(fill = 'grey80') +
  stat_summary(fun.data = 'mean_cl_boot') +
  labs(x = 'Assigned Y-haplotype', y = 'orange (% body area)') +
  theme_classic()

performance::icc(lme4::lmer(car_perc ~ 1 + (1 | hclust_groups), data = ind_df))
performance::icc(lme4::glmer(O6 ~ 1 + (1 | hclust_groups), data = ind_df, family = 'binomial'))
