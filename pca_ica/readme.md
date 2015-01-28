`run_pca` will perform PCA on a movie, and save the result:
- `compute_pca`: Carries out the actual PCA computation

`run_ica` will perform ICA on PCA filter-trace pairs, and save the result:
- `compute_spatiotemporal_ica_input`: Mixes the PCA filter-trace pairs according to the "spatiotemporal parameter" mu
- `compute_ica_weights`: Generates the ICA transformation matrix via FastICA
- `compute_ica_pairs`: Applies the ICA transformation matrix to the PCA pairs, to produce ICA pairs
