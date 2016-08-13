## Herrin analysis workflow

The Matlab-based analysis workflow for the strategy shifting experiment at Herrin is organized around a `DaySummary` object. The basic idea behind the `DaySummary` class is to encapsulate all relevant information for a single session of the strategy shifting protocol.

1. [`DaySummary` QuickStart](ds_quickstart.md)
  * Using `DaySummary` with non-PlusMaze datasets
  * Using `browse_rasters`
2. Classifying cells using `DaySummary`
3. [Cross-dataset alignment](alignment.md)
4. [Use of `MultiDay`](multiday.md)
  * Use of `browse_multiday`

## Standard data formats

This section describes the expected contents of stored data files.

#### "Rec" (mat) files:

Contains cell candidates (i.e. trace and filter pairs) to be classified. Stored as native Matlab "mat" file. (In the past, "rec" stood for "reconstruction", but this is now vestigial.) File name should be `rec_*.mat`.

Top level variables are: `info`, `filters`, `traces`, which must have the following contents:
```
info.type: Name of the cell extraction method (e.g. `ica`, `cellmax`, ...)
info.num_pairs: Number of filter-trace pairs in file

filters: [height x width x num_pairs]
traces: [num_frames x num_pairs]
```
