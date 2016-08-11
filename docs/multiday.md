## Quickstart to `MultiDay` usage

A single `DaySummary` instance encapsulates the information associated with a single session of the strategy shifting experiment. A `MultiDay` object is built on top of `DaySummary` and is intended to handle the low-level mechanics of multi-session cell alignment.

#### Instantiating a `MultiDay` instance

In this tutorial, we will create a `MultiDay` object from three `DaySummary` instances: `m1d12`, `m1d13`, and `m1d14` (i.e. Days 12, 13, 14 of Cohort 11, Mouse 1). We assume the following `DaySummary` variables exist in the Matlab workspace:
```
>> m1d12

m1d12 = 
  DaySummary with properties:
              cells: [1223x1 struct]
             trials: [110x1 struct]
          num_cells: 1223
         num_trials: 110
      trial_indices: [110x4 int32]
    full_num_frames: 13851

>> m1d13

m1d13 = 
  DaySummary with properties:
              cells: [1244x1 struct]
             trials: [155x1 struct]
          num_cells: 1244
         num_trials: 155
      trial_indices: [155x4 int32]
    full_num_frames: 19414

>> m1d14

m1d14 = 
  DaySummary with properties:
              cells: [1138x1 struct]
             trials: [100x1 struct]
          num_cells: 1138
         num_trials: 100
      trial_indices: [100x4 int32]
    full_num_frames: 11878
```

Additionally, in order to instantiate a `MultiDay` object based on the above `DaySummary` instances, we also need the "match matrices" between the days (see [here](alignment.md) for more information). We assume that the following match matrices exist in the workspace: `m_12to13`, `m_13to12`, `m_12to14`, `m_14to12`, `m_13to14`, `m_14to13`.

Instntiation of a `MultiDay` object requires a list of `DaySummary`s (`ds_list`) and a list of match matrices (`match_list`). They are formatted as follows:
```
>> ds_list = {12, m1d12; 13, m1d13; 14, m1d14}

ds_list = 
    [12]    [1x1 DaySummary]
    [13]    [1x1 DaySummary]
    [14]    [1x1 DaySummary]

>> match_list = {12, 13, m_12to13, m_13to12; 12, 14, m_12to14, m_14to12; 13, 14, m_13to14, m_14to13}

match_list = 
    [12]    [13]    {1223x1 cell}    {1244x1 cell}
    [12]    [14]    {1223x1 cell}    {1138x1 cell}
    [13]    [14]    {1244x1 cell}    {1138x1 cell}
```

With `ds_list` and `match_list`, the `MultiDay` object can be created as:
```
>> md = MultiDay(ds_list, match_list);
11-Aug-2016 12:54:27: Day 12 has 468 classified cells (out of 1223)
11-Aug-2016 12:54:27: Day 13 has 507 classified cells (out of 1244)
11-Aug-2016 12:54:27: Day 14 has 489 classified cells (out of 1138)
  Removed 4 inconsistent matches!
11-Aug-2016 12:54:27: Found 341 matching classified cells across all days

>> md

md = 
  MultiDay with properties:
         valid_days: [12 13 14]
           num_days: 3
          num_cells: 341
    matched_indices: [341x3 double]
           sort_day: 12
```

#### Important side note: Definition of cross-day match and "match conflicts"

#### Basic usage of the `MultiDay` object
