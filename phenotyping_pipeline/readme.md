# The phenotyping pipeline

This folder contains the scripts used to ingest, annotate and transform the guppy photos. The code used to actually train the deep learning models used in the pipeline can be found within the `phenotyping_models/` folder. The rest of the files here contain the code to build and run the pipeline itself. It contains
- `keras_bg_removal`: This Unet takes the original guppy image and predicts a mask of the guppy, allowing for background removal.
- `keras_carotenoid_detection`: This Unet takes extracted guppy images (i.e. without background), and predicts a mask of all carotenoid (orange and yellow) ornamentation on the fish.
- `keras_melanin_detection_v2`: Does the same as previous, but for the melanic (black) ornamentation.
- `keras_landmark_detection`: This convolutional net predicts heatmaps for likely locations of four fixed landmarks on the fish. Scripts in this folder convert landmark coordinates to heatmaps for this purpose. Taking the highest prediction in the output heatmap gives the most likely location for the landmark.

`new_photo_intake.R` was ran to detect newly added photos, run all the pipeline components on those new photos, add them to the database, and generate data files. This file therefore contains all the code to apply the pipeline to guppy images, and save the results. It uses the `phenotyping_models` and applies all the computational steps in between, such as image rotation, placemement of sliding landmarks etc.

`generate_consensus_shape.R`, `generate_warped_images.R` and `process_full_landmarks.R` were run to programmatically generate the consensus shape, the images warped to the consensus shape, and the sets of "full landmarks" (along the whole outline of the fish), respectively.

These files are used by the scripts above and contain various functions and tools: 

- `tools.R`
- `index_photos.R`
- `fish_extraction.R`
- `auto_place_landmarks.R`
- `carotenoid_extraction.R`
- `melanin_extraction_v2.R`
- `place_manual_landmarks.R`
