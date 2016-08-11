## Quickstart to `MultiDay` usage

A single `DaySummary` instance encapsulates the information associated with a single session of the strategy shifting experiment. A `MultiDay` object is built on top of `DaySummary` and is intended to handle the low-level mechanics of multi-session cell alignment.

#### Instantiating a `MultiDay` instance

In this tutorial, we will create a `MultiDay` object from three `DaySummary` instances. We assume the following `DaySummary` variables exist in the Matlab workspace: `m1d12`, `m1d13`, `m1d14` (i.e. Days 12, 13, 14 of Cohort 11 / Mouse 1).

Additionally, in order to instantiate a `MultiDay` object based on the above `DaySummary` instances, we also need the "match matrices" between the days. We assume that the following match matrices exist in the workspace: `m_12to13`, `m_13to12`, `m_12to14`, `m_14to12`, `m_13to14`, `m_14to13`. (See [here](alignment.md) for more information on match matrices.)

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

#### Definition of cross-day match

The following three-day example illustrates the behavior of the cross-day alignment in `MultiDay`. Suppose we have a single cell that shows up on three days as: Cell A from Day 1, Cell B from Day 2, and Cell C from Day 3. An edge (shown in red, below) between two cells indicates that there is a bidirectional match.

The simplest case is when there is a match between all day pairs, as follows:

![Full match](md_simple-case.png)

In this case, it clearly makes sense to call A / B / C to be an aligned cell.

However, the `MultiDay` matching algorithm would also consider the following relationships to indicate that A / B / C are aligned: 

![All match cases](md_all-cases.png)

Thus, it is possible to provide fewer than the full set of matches (in `match_list`) when instantiating `MultiDay`.

#### Match conflicts



#### Basic usage of the `MultiDay` object
