library(tidyverse)
library(patchwork)
library(legendry)
library(sf)
library(furrr)
source('sequencing/genomics_helpers.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')
source('paper_figures/theme.R')
source('sequencing/gwas/peak_viz/peak_viz_tools.R')

manhattan <- function(
    d,
    pval_column = 'p_SHet',
    gwas_colors = c('orange3', 'orange4'),
    buffer_dist = 0.012, #distance parameter passed to st_buffer(); determines point size; finicky AF
    hw = 2,
    chr_gap = 8e6
  ) {

  scaff_labs2 <- scaff_labs |>
    filter(chr2 != '12') |>
    mutate(
      cum_scaff_start = cumsum(lag(scaff_len + chr_gap, default = 1)),
      scaff_mid = cum_scaff_start + scaff_len / 2,
      scaff_color = ifelse(row_number() %% 2 == 0, gwas_colors[1], gwas_colors[2])
    )

  xpos <- 4e7
  ypos <- 5e7

  d2 <- d |>
    mutate(
      significant = p.value.adj < 0.05,
      facet = case_when(
        chr != 'NC_024342.1' & source %in% c('X', 'Y') ~ 'inset',
        TRUE ~ 'standard'
      ) |> factor(levels = c('standard', 'inset'))
    ) |>
    filter(!(chr2 %in% c('Un', 'MT'))) |>
    dplyr::select(-scaff_len, -cum_scaff_start, -cum_pos) |>
    left_join(scaff_labs2, join_by(chr, chr2)) |>
    mutate(
      facet = if_else(chr2 == 12, 'inset', facet),
      cum_pos = pos + cum_scaff_start,

      facet2 = if_else(
        facet == 'standard' & (source != 'Not called' | is.na(source)),
        'Autosomal', 'Sex-linked'
      ),
      cum_pos2 = case_when(
          facet2 == 'Autosomal' ~ cum_pos,
          chr2 == '12' ~ pos,
          source == 'X' ~ xpos,
          source == 'Y' ~ ypos
        ),
      colors = case_when(
        facet2 == 'Autosomal' ~ scaff_color,
        source == 'X' ~ 'firebrick',
        source == 'Y' ~ 'navy',
        chr2 == '12' ~ 'black',
      )
    ) |>
    arrange(desc(.data[[pval_column]]))

  pval_boundary <- filter(d2, significant) |> slice_max(.data[[pval_column]], with_ties = FALSE) |>
    pull(.data[[pval_column]])

  ####
  convert_to_sf <- function(data, x_conv, bdist) {
    data |>
      mutate(y = -log10(.data[[pval_column]]), x = cum_pos2 / x_conv) |>
      group_split(chr2) |>
      future_map(
        \(x, y) x |>
          st_as_sf(coords = c('x', 'y')) |>
          st_buffer(dist = bdist) |>
          st_union(),
        .progress = TRUE, .options = furrr_options(seed = NULL)
      ) |>
      reduce(st_union)
  }

  d3 <- filter(d2, facet2 == 'Autosomal')
  maxy <- max(-log10(d3[[pval_column]]))
  maxx <- max(d3$cum_pos2)
  auto_x_conversion <- maxx / maxy / hw
  auto_buffer <- buffer_dist * maxy

  odd_chrs <- d3 |> filter(chr2 %in% c(1, 3, 5, 7, 9, 11, 14, 16, 18, 20, 22)) |>
    convert_to_sf(auto_x_conversion, auto_buffer)
  even_chrs <- d3 |> filter(chr2 %in% c(2, 4, 6, 8, 10, 13, 15, 17, 19, 21, 23)) |>
    convert_to_sf(auto_x_conversion, auto_buffer)

  d4 <- filter(d2, facet2 == 'Sex-linked')
  maxy <- max(-log10(d4[[pval_column]]))
  maxx <- max(d4$cum_pos2)
  sex_x_conversion <- maxx / maxy / hw * 5
  sex_buffer <- buffer_dist * maxy

  black <- d4 |> filter(colors == 'black') |> convert_to_sf(sex_x_conversion, sex_buffer)
  red <- d4 |> filter(colors == 'firebrick') |> convert_to_sf(sex_x_conversion, sex_buffer)
  blue <- d4 |> filter(colors == 'navy') |> convert_to_sf(sex_x_conversion, sex_buffer)

  auto <- ggplot() +
    geom_sf(
      data = odd_chrs |> fortify() |> mutate(facet = 'Autosomal'),
      fill = gwas_colors[2], color = NA
    ) +
    geom_sf(
      data = even_chrs |> fortify() |> mutate(facet = 'Autosomal'),
      fill = gwas_colors[1], color = NA
    ) +
    coord_sf(expand = FALSE) +
    geom_hline(aes(yintercept = -log10(pval_boundary)), lty = 2, linewidth = 0.5) +
    labs(y = expression(-log[10](italic(p))), x = NULL) +
    scale_x_continuous(
      expand = expansion(c(0.01, 0.01)),
      breaks = scaff_labs2$scaff_mid[c(1:12, 14, 16, 18, 20, 22)] / auto_x_conversion,
      labels = scaff_labs2$chr2[c(1:12, 14, 16, 18, 20, 22)]
    ) +
    facet_grid(vars(facet)) +
    theme(
      strip.background = element_rect(fill = 'grey90', color = NA),
      strip.text = element_text(size = 7),
      axis.ticks.x = element_blank()
    )

  sex <- ggplot() +
    geom_sf(
      data = black |> fortify() |> mutate(facet = 'Sex-linked'),
      fill = 'black', color = NA
    ) +
    geom_sf(
      data = blue |> fortify() |> mutate(facet = 'Sex-linked'),
      fill = 'mediumblue', color = NA
    ) +
    geom_sf(
      data = red |> fortify() |> mutate(facet = 'Sex-linked'),
      fill = 'firebrick', color = NA
    ) +
    coord_sf(expand = FALSE) +
    geom_hline(aes(yintercept = -log10(pval_boundary)), lty = 2, linewidth = 0.5) +
    labs(y = NULL, x = NULL) +
    scale_x_continuous(
      limits = c(1, ypos + 5e6) / sex_x_conversion,
      breaks = c(13219565, xpos, ypos) / sex_x_conversion, labels = c('12', 'X', 'Y'),
      guide = guide_axis_nested(
        key = key_range_manual((xpos - 5e6) / sex_x_conversion, (ypos + 5e6) / sex_x_conversion, 'Cross-mapping')
      ),
      expand = expansion(mult = c(0.01, 0.2))
    ) +
    facet_grid(vars(facet)) +
    theme(
      strip.background = element_rect(fill = 'grey90', color = NA),
      strip.text = element_text(size = 7),
      axis.ticks.x = element_blank()
    )

  both <- auto | sex

  #ggsave('test.pdf', both, units = "cm", width = 8.8, height = 5, device = cairo_pdf)
  ####

  return(both)
}

car_PIE_gwas <- prep_sommer_gwas_table('car_PIE', name = NA, pval_column = 'p_SHet')
mel_PIE_gwas <- prep_sommer_gwas_table('mel_PIE', name = NA, pval_column = 'p_SHet')
car_orn_gwas <- prep_sommer_gwas_table(
  str_glue('car_{c(1:3, 5:7)}'), name = NA, pval_column = 'p.value'
)
mel_orn_gwas <- prep_sommer_gwas_table(
  str_glue('mel_{c(1, 3, 5:8)}'), name = NA, pval_column = 'p.value'
)

plan(multisession, workers = 16)

A <- manhattan(filter(car_orn_gwas, trait == 'car_2') |> slice_sample(n = 100000), pval_column = 'p.value')
B <- manhattan(filter(mel_orn_gwas, trait == 'mel_1') |> slice_sample(n = 100000), pval_column = 'p.value', gwas_colors = c('grey40', 'grey10'))
C <- manhattan(car_PIE_gwas |> slice_sample(n = 100000))
D <- manhattan(mel_PIE_gwas |> slice_sample(n = 100000), gwas_colors = c('grey40', 'grey10'))

A <- manhattan(filter(car_orn_gwas, trait == 'car_2'), pval_column = 'p.value')
B <- manhattan(filter(mel_orn_gwas, trait == 'mel_1'), pval_column = 'p.value', gwas_colors = c('grey40', 'grey10'))
C <- manhattan(car_PIE_gwas)
D <- manhattan(mel_PIE_gwas, gwas_colors = c('grey40', 'grey10'))

fig5 <- wrap_plots(
  A, B, C, D,
  ncol = 2
) + plot_annotation(tag_levels = 'A')

plan(sequential)

ggsave('paper_figures/Fig5.pdf', fig5, units = "cm", width = 20, height = 10, device = cairo_pdf)
