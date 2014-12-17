analysis
========

Repository for code that is _specific_ to the strategy-shifting experiment at Herrin. (General data preprocessing code should go into the lab-wide repository!)

We should strictly observe the lab-wide formatting for data objects:
- Movies: `[height width time]`
- Traces: `[cell-idx time]`
- Events (binarized): `[cell-idx time]`
- Events (binarized, sparse): Struct with fields: `nFrames: num_frames, eventTimes: {cell-idx, [events]}`
- Cell images: `[height width cell-idx]`
