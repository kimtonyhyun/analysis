### Quickstart for DaySummary usage

The Matlab-based analysis workflow for the strategy shifting experiment at Herrin is based on the use of a `DaySummary` object. The basic idea behind the `DaySummary` class is to encapsulate all relevant information for a single session of the strategy shifting protocol. This tutorial will explain the basic instantiation and usage of the `DaySummary` object.

#### Preferred folder structure for session data

The following screenshot illustrates my (Tony's) preferred folder structure for session data, using the "c11m1d12" dataset as an example:

![Tony's preferred folder structure](https://raw.githubusercontent.com/schnitzer-lab/analysis/kimth/ds-docs/docs/ds_folder_structure.PNG?token=AB_C3xtbUy5yEYCYtQyMMzBJkIygOKY_ks5XQzrcwA%3D%3D)

The experimental data is stored in the `_data` subdirectory:

- __"c11m1d12_ti2.txt"__: The "PlusMaze text file" containing trial metadata (e.g. for each trial mouse start, goal, choice, frame indices). _Necessary_ for instantiating a `DaySummary` object.
- __"c11m1d12_ti2.mp4"__: The overhead behavioral video associated with the session. Optional input for instantiating `DaySummary`.
- __"c11m1d12_ti2.xy"__: Mouse spatial tracking data. Each line is the (x,y) coordinate of the mouse position for a given frame. Optional input for instantiating `DaySummary`.
- __"c11m1d12_gfix_rm_mc_cr_norm40_dff_ti2.hdf5"__: The DFF movie associated with the session. Not an input to `DaySummary` (but is used during manual cell classification).
 
Cell extraction data is stored in a separate `cm01` subdirectory:

- __"rec_151125-102235.mat"__: Matlab data file (we'll call it a "rec" file), containing the cell filters and traces as produced by the cell extraction algorithm (here, CELLMax). _Necessary_ for instantiating a `DaySummary` object.
- __"class_151125-172203.txt"__: Text file that contains the classification label for each of the candidate sources (filter-trace pair).
 
Note: The bare minimum for instantiating a `DaySummary` object is the PlusMaze text file ("c11m1d12_ti2.txt" in the above example) and the filter-trace Matlab data file ("rec_151125-102235.mat").

#### Instantiating the `DaySummary` object

For convenience, I like to hard-code a "sources" struct for each session dataset, so that I don't have to spell out the individual filenames each time:
```
>> sources = data_sources
sources = 
         maze: '_data/c11m1d12_ti2.txt'
     behavior: '_data/c11m1d12_ti2.mp4'
     tracking: '_data/c11m1d12_ti2.xy'
    miniscope: '_data/c11m1d12_gfix_rm_mc_cr_norm40_dff_ti2.hdf5'
```

The "sources" struct above can be directly used for the instantiation of a `DaySummary` object. The struct _must_ have a `maze` subfield with the name of the PlusMaze text file. The `behavior` and `tracking` subfields are optional, and the `DaySummary` constructor will automatically load them as well when available. The `miniscope` subfield is not used by `DaySummary`.

The `DaySummary` object is instantiated as follows:
```
>> m1d12 = DaySummary(sources, 'cm01');
16-May-2016 10:39:59: Loaded trial metadata from _data/c11m1d12_ti2.txt
16-May-2016 10:40:06: Loaded 1223 filters and traces (cellmax) from cm01\rec_151125-102235.mat
  Computing auxiliary spatial parameters... Done (25.2 sec)
  Computing distances between all sources... Done (1.9 sec)
16-May-2016 10:40:43: Loaded classification from cm01\class_151125-172203.txt
16-May-2016 10:40:45: Loaded behavior video from "_data/c11m1d12_ti2.mp4"
16-May-2016 10:40:45: Loaded tracking data from "_data/c11m1d12_ti2.xy"
```

The second argument to the `DaySummary` constructor (`'cm01'` in the above example) indicates the subdirectory containing the "rec" file (i.e. cell filter and trace data). The classification file will be automatically loaded if available in the same subdirectory.

#### Minimal instantiation of the `DaySummary` object

The above example used the "sources" struct to automatically load the behavior and tracking data. Below, I demonstrate the "minimal" instantiation option, where only the PlusMaze text file is provided:
```
>> m1d12 = DaySummary('_data/c11m1d12_ti2.txt', 'cm01');
16-May-2016 10:46:28: Loaded trial metadata from _data/c11m1d12_ti2.txt
16-May-2016 10:46:35: Loaded 1223 filters and traces (cellmax) from cm01\rec_151125-102235.mat
  Computing auxiliary spatial parameters... Done (25.3 sec)
  Computing distances between all sources... Done (2.0 sec)
16-May-2016 10:47:12: Loaded classification from cm01\class_151125-172203.txt
```
