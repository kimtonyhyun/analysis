## Standard data formats

This document describes the expected contents of stored data files.

#### PCA mat files:

Stores the results of PCA decomposition of the movie, used by a number of different cell extraction methods (e.g. ICA, SSS). File name should be `pca_*.mat`.

Top level variables are: `info`, `filters`, `traces`, `S`.
```
info.movie_height
info.movie_width
info.movie_frames
info.num_PCs

info.trim.enabled: "Trim" removes certain (presumably non-cell) pixels from the PCA computation
info.trim.idx_kept: Indices to the "kept" pixels after the trim operation

info.medfilt.enabled
info.medfilt.halfwidth

filters: [num_PCs x num_pixels]
traces: [num_PCs x num_frames]
S: 1-D vector of singular values from the PCA decomposition
```

#### "Reconstructed" mat files:

Contains cell candidates (i.e. trace and filter pairs) to be classified. For historical reasons, referred to as "reconstructions." File name should be `rec_*.mat`.

Top level variables are: `info`, `filters`, `traces`
```
info.type: Name of the cell extraction method (e.g. `ica`, `reconstruction`, ...)
info.num_pairs: Number of filter-trace pairs in file

filters: [height x width x num_pairs]
traces: [num_frames x num_pairs]
```
