library(tidyverse)
source('sequencing/genomics_helpers.R')

files <- list.files('sequencing/gwas_sommer/snp_inheritance/sommer_snp_inheritance/female', f = T)
files <- str_subset(files, '342', negate = TRUE)

link <- map(files, data.table::fread) |> list_rbind()

count(link, source) |> mutate(prop = n / sum(n), perc = prop * 100) |> print()

included_loci <- load_vcf(reference = 'female')$df

link2 <- semi_join(link, included_loci)

count(link2, source) |> mutate(prop = n / sum(n), perc = prop * 100) |> print()
