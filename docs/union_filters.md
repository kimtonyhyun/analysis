## Generating cross-day "union" filters

At a high-level, the generation of "union" filters across multiple days consists of the following steps:

1. Run alignment (`run_alignment.m`) across all days.
2. For each day, import (potentially) missed filters from the other sessions. Each import will generate a `rec_*.mat` file.
3. Sort all imported filter sets.
4. "Resolve" possible duplicates.
