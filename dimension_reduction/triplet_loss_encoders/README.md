# Phenotype embeddin using triplet learning

The R scripts in this folder were used to generate the patternspace phenotype embeddings.

`fit_triplet_loss_model.R`: This file contains the function that trains the network. It takes the images of patterns, combines it with pedigree information, generates triplets, and trains a network on those triplets. It will save the training history, and the model weights to disk.

`car_model_ped_fullcolor_comparison.R`: This will train 5 networks on the orange patterns, from 1 dimensional patternspace up to 5 dimensional patternspace.

`mel_model_ped_fullcolor_comparison.R`: This will train 5 networks on the black patterns, from 1 dimensional patternspace up to 5 dimensional patternspace.

`triplet_loss_tools.R` This file contains helper functions that are used in the other files.

`embeddings`: The actual embeddings are saved in this folder. Only the 5 dimensional embeddings are used in the paper.
