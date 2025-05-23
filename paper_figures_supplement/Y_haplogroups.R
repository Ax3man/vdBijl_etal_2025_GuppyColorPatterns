library(irlba)
library(tidyverse)
library(ggtext)

source('sequencing/genomics_helpers.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')

# Import the VCF file into a GDS file
vcf_file <- "sequencing/gwas_sommer/filtered2.vcf.gz"

g <- prep_sommer_gwas_table('car_PIE', pval_column = 'p_SHet') |> filter(significant, source == 'Y')

if (FALSE) { # only run interactively, otherwise load from disk
  geno <- load_vcf(region = str_glue_data(g, '{chr}:{pos}-{pos}'))

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
  labs(x = 'Principal component', y = 'Proportion of variance')

# hclust
dist_matrix <- dist(pca$x) # distance matrix from the first two PCs
hc <- hclust(dist_matrix)
plot(hc, FALSE) # to visualize the dendrogram

##
source('quant_gen/prepare_pedigrees.R')
source('sequencing/gwas/peak_viz/peak_viz_tools.R')
source('selection_decisions/compile_decisions.R')
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
orn <- read_rds('ornament_analysis/car_ornaments.rds') %>%
  left_join(data.table::fread('photo_database.csv') %>% dplyr::select(unique_id, fish_id)) %>%
  group_by(fish_id, ornament) %>%
  summarise(present = Mode(present)) %>%
  pivot_wider(names_from = ornament, values_from = present)

ind_df <- as.data.frame(pca$x) %>% rownames_to_column('sample_name') %>% as_tibble() %>%
  left_join(get_sampling_structure(), join_by(sample_name)) %>%
  left_join(mutate(selection, fish_id = tolower(fish_id)), join_by(fish_id)) %>%
  left_join(orn, join_by(fish_id)) %>%
  add_patriline(ped_df)

ind_df$hclust_groups <- factor(cutree(hc, k = 4))

ind_df %>%
  mutate(Y_haplogroup = factor(hclust_groups, levels = 1:4, labels = paste0('Y', 1:4))) %>%
  dplyr::select(sample_name, fish_id, Y_haplogroup) %>%
  data.table::fwrite('sequencing/Y_haplogroups.csv')

ind_df %>%
  pivot_longer(car_1:car_7, names_to = 'ornament', values_to = 'presence') %>%
  mutate(ornament_label = glue::glue(
    "<img src='ornament_analysis/ornament_figure_labels/{ornament}.png' width='60' />"
  )) %>%
  ggplot(aes(hclust_groups, fill = factor(presence, 0:1, c('absent', 'present')))) +
  geom_bar() +
  scale_y_continuous(expand = c(0, 0), name = 'Count') +
  scale_fill_manual(values = c('grey20', 'darkorange'), name = NULL) +
  facet_grid(cols = vars(ornament_label)) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_markdown()) +
  labs(x = 'Assigned Y-haplogroup')

ggsave('paper_figures_supplement/Y_haplogroup_vs_ornaments.png', w = 8, h = 3.5)

biplot <- ggplot(ind_df, aes(PC1, PC2, color = hclust_groups)) +
  geom_point() + theme_classic() + coord_fixed() +
  scale_color_brewer(palette = "Set1", labels = paste0('Y', 1:4)) +
  labs(
    x = str_glue('PC1 ({scales::percent_format()(summary(pca)$importance[2, 1])})'),
    y = str_glue('PC2 ({scales::percent_format()(summary(pca)$importance[2, 2])})'),
    color = 'Assigned\nY-haplogroup',
    #subtitle = 'PCA on SNPs associated w/ orange pattern and Y-like inheritance'
  ) +
  theme(legend.position = 'top')

patchwork::wrap_plots(biplot, scree, ncol = 2) + patchwork::plot_annotation(tag_levels = 'A')
ggsave('paper_figures_supplement/Y_haplogroups.png', w = 8, h = 3.5)

