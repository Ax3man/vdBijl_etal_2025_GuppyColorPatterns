# Plot the results in mini-heatmaps ----------------------------------------------------------------
library(tidyverse)
library(imager)

prep_for_minis <- function(trait = 'orange_ornaments') {
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
}

make_mini <- function(snp_effects, chr, pos, upper_limit = 3) {
  bg <- load.image('data/extracted_fish_warped/replicate_3/gen_4/20221024_IMG_2089.png') |>
    imager::channel(4) |> resize_halfXY() |> as.data.frame() |> filter(value > 0)

  # show successful models
  # snp_effects |> dplyr::select(x, y) |> distinct() |>
  #   ggplot(aes(x, y, fill = 1)) +
  #   geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
  #   geom_raster() + scale_y_reverse() + coord_fixed(expand = FALSE)

  minimal_plot <- function(x, upper_limit) {
    ggplot(x, aes(x, y, fill = abs(Z_score))) +
      geom_raster(aes(alpha = value), data = bg, fill = 'grey80', show.legend = FALSE) +
      geom_raster() +
      scale_y_reverse() +
      scale_fill_viridis_c(limits = c(0, upper_limit), oob = scales::oob_squish) +
      coord_fixed(expand = FALSE) +
      theme_void()
  }
  f <- filter(snp_effects, .data$chr == .env$chr, .data$pos == .env$pos)
  if (n_distinct(f$ref) > 1 || n_distinct(f$alt) > 1) {
    f <- f |> filter(ref == last(ref), alt == last(alt))
  }
  minimal_plot(f, upper_limit)
}


