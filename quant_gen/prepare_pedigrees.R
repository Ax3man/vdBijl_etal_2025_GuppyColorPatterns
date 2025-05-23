suppressPackageStartupMessages({
  library(nadiv)
  library(tidyverse)
})

add_patriline <- function(data, pedigree) {
  # Create new column with the patriline, empty for now
  data <- mutate(data, patriline = NA)
  # Loop through the pedigree.
  for (i in seq_len(nrow(pedigree))) {
    id <- pedigree[i, 'animal', drop = TRUE]
    # If the fish is female, we use NA
    if (pedigree[i, 'sex', drop = TRUE] == 'F') {
      current_id <- NA
    } else {
      # Initialize the pedigree chain, with the current individuals sire
      current_id <- id
      while (TRUE) {
        # Find the sire
        sire <- pedigree[pedigree$animal == current_id, 'sire', drop = TRUE]
        # if the new sire is NA, break the while loop and keep the previous id
        if (is.na(sire)) break
        # Otherwise change the id to be one level up, so from son to father and look for more sires
        current_id <- sire
      }
    }
    # assign current_id as the patriline for all photos of the current iteration
    data$patriline[data$fish_id == id] <- current_id
  }
  return(data)
}

# Load the pedigree data
ped_df <- data.table::fread('data/pedigree.csv') %>%
  as_tibble() %>%
  mutate(
    sex = str_sub(animal, 1, 1),
    across(c(animal, sire, dam), tolower),
    grow_tank = ifelse(is.na(sire), 'source_pop', paste(sire, dam, date_of_birth, sep = '_'))
  )

pedA <- makeA(as.data.frame(dplyr::select(ped_df, ID = animal, Dam = dam, Sire = sire)))
A <- as.matrix(pedA)
pedA_inv <- makeAinv(as.data.frame(dplyr::select(ped_df, animal, dam, sire)))$Ainv
colnames(pedA_inv) <- rownames(pedA_inv)
A_inv <- as.matrix(pedA_inv)

# Make the additive relationship matrix for the X-chromosome
suppressWarnings({
  S <- makeS(as.data.frame(dplyr::select(ped_df, animal, dam, sire, sex)), 'M', 'ngdc', TRUE)
})
pedX <- S$S
colnames(pedX) <- rownames(pedX)
X <- as.matrix(pedX)
pedX_inv <- S$Sinv
colnames(pedX_inv) <- rownames(pedX_inv)
X_inv <- as.matrix(pedX_inv)

# Make the additive relationship matrix for the Y.
makeY <- function(ped_df) {
  df <- ped_df |>
    filter(sex != 'F') |>
    rename(fish_id = animal) |>
    add_patriline(ped_df)
  Y <- outer(df$patriline, df$patriline, FUN = function(x, y) as.numeric(x == y))
  rownames(Y) <- colnames(Y) <- df$fish_id

  return(Y)
}
Y <- makeY(ped_df)

