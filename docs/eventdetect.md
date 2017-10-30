## Detecting events from fluorescence traces

Begin by instantiating a `DaySummary` object that contains the traces for event detection (see the [`DaySummary` Quickstart](docs/ds_quickstart.md) for background info):
```
>> sources = data_sources
sources = 
  struct with fields:
         maze: '_data/c14m4d15_ti2.txt'
     behavior: '_data/c14m4d15_ti2.mp4'
     tracking: '_data/c14m4d15_ti2.xy'
    miniscope: '_data/c14m4d15_gfix_rm_cr_mc_cr_norm40_dff_ti2.hdf5'

>> ds = DaySummary(sources.maze, 'cm01-fix');
30-Oct-2017 11:59:38: Loaded 364 filters and traces (export_rec) from cm01-fix\rec_161217-182928.mat
30-Oct-2017 11:59:38: Loaded trial metadata from _data/c14m4d15_ti2.txt
  Computing trace correlations between all sources... Done (0.3 sec)
  Computing auxiliary spatial parameters... Done (5.3 sec)
  Computing distances between all sources... Done (0.3 sec)
30-Oct-2017 11:59:47: Loaded classification from cm01-fix\class_161217-183059.txt
>> 
```

Remarks:
- Do not apply probe trial elimination (e.g. the `'noprobe'` flag) on `DaySummary` instantiation. Global trace properties (e.g. the baseline standard deviation) will be affected by omission of trials. Let's pre-emptively avoid future confusion by always applying event detection on the full trace from each session.
- Event detection will be performed for _all_ sources in the `DaySummary`, whether or not the source has been classified to be a cell. For this reason, it's advised to perform event detection on `DaySummary` instances containing only classified cells.

(In principle, event detection may be performed on the raw fluorescence traces without using the `DaySummary` wrapper. However, the `DaySummary` organization provides a convenient interface for accessing relevant behavioral parameters (e.g. frame indices associated with trial) which we'll make use of.)
