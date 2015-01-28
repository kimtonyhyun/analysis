### Manual classification

The top level script for manual classification is `classify_ics`, which:

1. First, examine the temporal IC trace via `view_trace`:
  * `view_trials_in_trace_by_color`: Color-code an IC trace by trial
  * `view_superimposed_trials_in_trace`: Look for trial-phase sensitivity of an IC trace
2. Second, examine the spatial IC filter via `view_ic_over_movie_interactively`:
  * `threshold_ic_filter`: Generates the outline of the IC spatial filter
  * `parse_active_frames`: Highlights portions of the trace that likely contain activity
 
In addition, support functions `save_classification` and `load_classification` allow the user to continue manual classification over multiple sittings.
