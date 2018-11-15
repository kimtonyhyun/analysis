## Explanation of PlusMaze metadata (txt) file format

An example PlusMaze metadata (txt) file is here:

The first five lines of the example file shows:
```
east north north 14.224 1 105 183 285
east north north 13.010 286 394 443 545
west south south 13.122 546 654 706 808
west south south 12.905 809 916 964 1066
east north north 12.675 1067 1175 1217 1319
...
```

Each line of the txt file represents a trial, in the format:
```
<start-arm> <goal-arm> <choice-arm> <time> <frame-idx1> <frame-idx2> <frame-idx3> <frame-idx4>
```
where:
- `<*-arm>` is one of `{'east', 'west', 'north', 'south'}`. Choice arm is the one selected by the animal. The trial is considered to be correct if the goal and choice arms are equal.
- `<time>` is the time (in seconds) of the trial. This parameter is not actually used by `DaySummary` for computing any derivative quantities.
- `<frame-idx1>` is the frame index corresponding to the _first frame_ of the trial.
- `<frame-idx2>` is the frame index corresponding to the _opening_ of the _start_ gate.
- `<frame-idx3>` is the frame index corresponding to the _closing_ of the _end_ gate.
- `<frame-idx4>` is the frame index corresponding to the _last frame_ of the trial.
