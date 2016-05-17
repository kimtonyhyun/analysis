## Use of `DaySummary` objects

The Matlab-based analysis workflow for the strategy shifting experiment at Herrin is organized around a `DaySummary` object. The basic idea behind the DaySummary class is to encapsulate all relevant information for a single session of the strategy shifting protocol.

1. [`DaySummary` QuickStart](ds_quickstart.md)
2. Classifying cells using `DaySummary`
3. Cross-dataset alignment
4. Use of `MultiDay`

## Standard data formats

This section describes the expected contents of stored data files.

#### "Reconstructed" mat files:

Contains cell candidates (i.e. trace and filter pairs) to be classified. For historical reasons, referred to as "reconstructions." File name should be `rec_*.mat`.

Top level variables are: `info`, `filters`, `traces`, which must have the following structure:
```
info.type: Name of the cell extraction method (e.g. `ica`, `cellmax`, ...)
info.num_pairs: Number of filter-trace pairs in file

filters: [height x width x num_pairs]
traces: [num_frames x num_pairs]
```
