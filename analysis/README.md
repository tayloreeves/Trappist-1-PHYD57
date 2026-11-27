# analysis

Stability check, edge-case runs and plots.

## Using the stability check

```cpp
#include "analysis/stability.hpp"

StabilityLimits L = make_limits(x, y, z);              // once per copy, at the start

StabilityResult r = check_stability(x, y, z, m, L);    // every few steps
if (r.code != STABLE) {
    // save copy_id, time, stability_name(r.code), r.i, r.j
}
```

Index 0 is the star, 1-7 are planets b-h. Positions as separate x, y, z arrays.
Works in code units and also compiles with nvcc.

## Results format

Please save results as `results.csv` with these columns so the plots work:

```
copy_id,unstable,t_instability,reason,planet_i,planet_j
```

Time in years. Reason is stable / close_encounter / escape / infall.

## Files

- `stability.hpp` - the stability check
- `test_stability.cpp` - tests for it
- `sandbox.cpp` - small integrator for edge cases
- `trappist_params.py` - planet data, resonances, Hill spacing
- `plot_results.py` - plots for the main runs
- `plot_edgecases.py` - resonance + mass scaling figure
- `results/` - edge-case outputs and figure
