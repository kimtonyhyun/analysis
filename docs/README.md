## Standard data formats

This document describes the expected contents of stored data files.

#### "Reconstructed" mat files:

Contains cell candidates (i.e. trace and filter pairs) to be classified. For historical reasons, referred to as "reconstructions." File name should be `rec_*.mat`.

Top level variables are: `info`, `filters`, `traces`
```
info.type: Name of the cell extraction method (e.g. `ica`, `cellmax`, ...)
info.num_pairs: Number of filter-trace pairs in file

filters: [height x width x num_pairs]
traces: [num_frames x num_pairs]
```
