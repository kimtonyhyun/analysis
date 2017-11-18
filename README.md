analysis
========

Repository for code that is specific to the strategy-shifting experiment at Herrin.

Basic formatting of data objects:
- Movies: `[height width time]`
- Traces: `[time cell-idx]`
- Cell images: `[height width cell-idx]`

Some required toolboxes:
- Bioinformatics toolbox: for `MultiDay` instantiation (uses `graphconncomp`)
- Signal Processing toolbox: for event detection (uses `butter`, `filtfilt`)
