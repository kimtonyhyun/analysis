## Quickstart to `MultiDay` usage

A single `DaySummary` instance encapsulates the information associated with a single session of the strategy shifting experiment. A `MultiDay` object is built on top of `DaySummary` and is intended to handle the low-level mechanics of multi-session cell alignment.

#### Instantiating a `MultiDay` instance

In this tutorial, we will create a `MultiDay` object from three `DaySummary` instances: `m1d12`, `m1d13`, and `m1d14`; i.e. Days 12, 13, 14 of Cohort 11, Mouse 1. We assume the following `DaySummary` variables exist in the Matlab workspace:
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
