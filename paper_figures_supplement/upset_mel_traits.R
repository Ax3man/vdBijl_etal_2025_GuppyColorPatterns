library(ComplexUpset)
library(tidyverse)

source('paper_figures/theme.R')
source('sequencing/gwas_sommer/gwas_sommer_tools.R')

# first, load GWAS results, and calculate q-values:
main <- prep_sommer_gwas_table('mel_PIE', pval_column = 'p_SHet') |>
  filter(significant)

components <- prep_sommer_gwas_table(str_glue('mel_{c(1, 3, 5:8)}'), 'black_ornaments') |>
  mutate(trait = paste0('B', parse_number(trait))) |>
  filter(significant)

table(components$trait)
nrow(components) / nrow(main)

# pivot and combine to get the 0/1 data for upset
to_plot <- bind_rows(
  main |> mutate(trait = 'pattern'),
  components
) |> distinct() |>
  dplyr::select(trait, chr, pos, alt, source) |>
  as_tibble() |>
  mutate(tmp = 1L) |>
  pivot_wider(names_from = trait, values_from = tmp, values_fill = 0L) |>
  drop_na(source) |>
  mutate(source = ifelse(source == 'Not called', 'auto', source))

# Summarise the amount of cross-mapping variants for the pattern
filter(to_plot, chr != 'NC_024342.1', str_starts(chr, 'NC'), pattern == 1) |>
  dplyr::count(source) |> mutate(frac = n / sum(n))
# source     n   frac
#  X         46 0.0817
#  Y        252 0.448
#  auto     265 0.471
# out of    563

make_upset <- function(d, min_size, k_labels = FALSE) {
  if (k_labels) {
    labels <- scales::label_comma(scale = 1/1000, suffix = 'k')
  } else {
    labels <- scales::label_comma()
  }

  u <- upset(
    as.data.frame(d),
    intersect = rev(c('pattern', paste0('B', c(1, 3 , 5:8)))),
    min_size = min_size, name = NULL, sort_sets = FALSE,

    base_annotations = list(
      Inheritance = (
        ggplot(mapping = aes(fill = factor(source, rev(c('Y', 'X', 'auto'))))) +
          geom_bar() +
          scale_fill_manual(
            values = c(auto = 'black', 'X' = 'firebrick', 'Y' = 'blue3'),
            limits = c('auto', 'X', 'Y'),
            labels = c('Autosomal', 'X-linked', 'Y-linked'),
          ) +
          scale_y_continuous(labels = labels) +
          labs(y = 'Intersection size', fill = 'Inferred pattern\nof inheritance\n(>80% support)')
      )
    )
  ) & theme_get()

  eb <- element_blank()
  u[[2]] <- u[[2]] + theme(
    panel.grid.major.y = eb, axis.text.x = eb, axis.ticks.x = eb, axis.title.x = eb, axis.line.x = eb
  )
  u[[3]]$layers[[1]] <- NULL

  u[[3]] <- u[[3]] +
    aes(fill = factor(source, rev(c('Y', 'X', 'auto')))) +
    scale_fill_manual(
      values = c(auto = 'black', 'X' = 'firebrick', 'Y' = 'blue3'),
      limits = c('auto', 'X', 'Y'),
      labels = c('Autosomal', 'X-linked', 'Y-linked'),
      guide = 'none'
    ) +
    theme(
      axis.text.y = eb, axis.ticks.y = eb, axis.title.y = eb, axis.line.y = eb
    ) +
    scale_y_reverse(labels = labels)
  u[[4]]$layers[[1]] <- NULL
  u[[4]] <- u[[4]] + theme(
    axis.text.x = eb, axis.ticks.x = eb, axis.title.x = eb, axis.line.x = eb
  ) + ylab(NULL)

  return(u)
}

upset_all <- make_upset(to_plot, min_size = 10, k_labels = TRUE)
upset_auto <- make_upset(
  dplyr::filter(to_plot, chr != 'NC_024342.1', str_starts(chr, 'NW', TRUE)),
  min_size = 1
)

# what proportion of ornament loci were captured by the pattern?
to_plot |> dplyr::filter(B1 | B3 | B5 | B6 | B7 | B8) |>
  dplyr::count(pattern = as.logical(pattern)) |> mutate(frac = n / sum(n))
# what proportion of pattern loci were captured by the ornaments?
to_plot |> dplyr::filter(pattern == 1) |>
  dplyr::count(ornament = B1 | B3 | B5 | B6 | B7 | B8) |> mutate(frac = n / sum(n))


ggsave('paper_figures_supplement/upset_mel_traits.png', upset_all, width = 6, height = 4)
ggsave('paper_figures_supplement/upset_mel_traits_autosomal.png', upset_auto, width = 6, height = 4)
