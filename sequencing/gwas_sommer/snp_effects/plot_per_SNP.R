# Plot the results in mini-heatmaps ----------------------------------------------------------------

trait <- 'orange_ornaments'

library(tidyverse)
library(imager)

f <- list.files(str_glue('sequencing/gwas_sommer/snp_effects/pixel_results/{trait}'), full.names = TRUE)
# print(length(f))
names(f) <- f |>
  str_remove(fixed(
    str_glue('sequencing/gwas_sommer/snp_effects/pixel_results/{trait}/female_pixel_') |> as.character()
  )) |>
  str_remove(fixed('.csv.gz'))

snp_effects <- map(f, data.table::fread) |>
  list_rbind(names_to = 'pixel') |>
  separate(pixel, c('x', 'y'), sep = '_', convert = TRUE) |>
  mutate(
    Z_score = beta / beta.se,
    SNP = str_glue('{chr}_{pos}_{alt}')
  )

bg <- load.image('data/extracted_fish_warped/replicate_3/gen_4/20221024_IMG_2089.png') |>
  channel(4) |> resize_halfXY() |> as.data.frame() |> filter(value > 0)

# show successful models
# snp_effects |> dplyr::select(x, y) |> distinct() |>
#   ggplot(aes(x, y, fill = 1)) +
#   geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
#   geom_raster() + scale_y_reverse() + coord_fixed(expand = FALSE)

minimal_plot <- function(x, upper_limit) {
  # p <- ggplot(x, aes(x, y, fill = -log10(p.value))) +
  #   geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
  #   geom_raster() +
  #   scale_y_reverse() +
  #   scale_fill_viridis_c(
  #     name = expression(-log[10](italic(P))),
  #     limits = c(3, -log10(min(snp_effects$p.value))), na.value = 'grey80',
  #     option = 'C'
  #   ) +
  #   coord_fixed(expand = FALSE) +
  #   theme_void()# +
  #guides(fill = 'none')

  p <- ggplot(x, aes(x, y, fill = abs(Z_score))) +
    geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_viridis_c(limits = c(0, 5), oob = scales::oob_squish) +
    coord_fixed(expand = FALSE) +
    theme_void()

  filename <- str_glue('sequencing/gwas_sommer/snp_effects/per_snp_plots/{trait}/{x$SNP[1]}.png')
  ggsave(filename, p, w = 4, h = 2)
}

walk(split(snp_effects, snp_effects$SNP), minimal_plot)
