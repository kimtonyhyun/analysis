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

- __"rec_151125-102235.mat"__: Matlab data file, containing the cell filters and traces as produced by the cell extraction algorithm (here, CELLMax). _Necessary_ for instantiating a `DaySummary` object.
- __"class_151125-172203.txt"__: Text file that contains the classification label for each of the candidate sources (filter-trace pair).
 
Note: The bare minimum for instantiating a `DaySummary` object is the PlusMaze text file ("c11m1d12_ti2.txt" in the above example) and the filter-trace Matlab data file ("rec_151125-102235.mat").

