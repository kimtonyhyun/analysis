## Cross-dataset alignment using DaySummary

In this tutorial, we assume two `DaySummary` objects (`m1d12` and `m1d13`) that represent the two datasets to be aligned:
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
```

The alignment procedure can be invoked as:
```
>> [match_12to13, match_13to12] = run_alignment(m1d12, m1d13);
```

At this point, you are prompted to select 4 pairs of cells from each dataset. I find it easiest to alternate clicks between the datasets:
