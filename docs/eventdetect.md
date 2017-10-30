## Detecting events from fluorescence traces

Begin by instantiating a `DaySummary` object that contains the traces for event detection (see the [`DaySummary` Quickstart](docs/ds_quickstart.md) for background info):

(In principle, event detection may be performed on the raw fluorescence traces without using the `DaySummary` wrapper. However, the `DaySummary` organization provides a convenient interface for accessing relevant behavioral parameters (e.g. frame indices associated with trial) which we'll make use of.)
