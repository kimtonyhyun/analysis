## Using DaySummary with non-PlusMaze datasets

The `DaySummary` object has been retrofitted so that it can encapsulate a dataset (i.e. set of filters and traces) with no associated "PlusMaze text file" (e.g. "c11m1d12.txt", as described in the [DaySummary Quickstart](ds_quickstart.md)).

To instantiate the DaySummary instance this way, call the constructor as follows:
```
ds = DaySummary([], 'path_to_rec_dir')
```
