# analysis

Stability check, edge-case runs and plots.

## Using the stability check (stability.hpp)

```cpp
#include "analysis/stability.hpp"

// once per copy, from the starting state (same layout as the integrator: x[nd*i+k])
StabLimits L = stab_setup(x, v, m, R, nd, G);     // R = physical radii, G = 4 pi^2 in AU/yr/Msun

// every few steps
StabResult r = stab_check(x, m, L);
if (r.code != STABLE) {
    // save copy_id, time, stab_name(r.code), r.i, r.j and stop this copy
}
if (r.hill_flags) {
    // a planet pair is inside its mutual Hill radius: log it, keep going
}
```

Body 0 is the star, 1-7 are planets b-h. Works in 2D or 3D, and compiles with nvcc for the GPU.

## What it checks

Stops the run:
- collision: two bodies touch (d < R_i + R_j)
- tidal_disruption: inside the Roche limit of the star or another planet, d = 2.44 R_primary (rho_primary/rho_secondary)^(1/3)
- close_encounter: two planets within 2 Roche-lobe radii (assignment; `L.use_close_encounter = false` to turn off)
- escape: further than 2x planet h's starting distance (assignment)
- infall: inside planet b's starting periapsis minus 2 of its Hill radii (assignment; `L.use_infall = false` to turn off)

Flags only (run continues):
- hill_flags: bit `pair_index(i,j)` set while planets i and j are inside their mutual Hill radius

## Results format

```
copy_id,unstable,t_instability,reason,planet_i,planet_j
```

Time in years. For star-planet events planet_i = 0.

## Files

- `stability.hpp` - the stability check
- `test_stability.cpp` - tests (`g++ -O2 -std=c++17 test_stability.cpp -o test_stability`)
- `sandbox.cpp` - small Forest-Ruth integrator for edge cases
- `trappist_params.py` - planet data, resonances, Hill spacing
- `plot_results.py` - plots for the main runs
- `plot_edgecases.py` - resonance + mass scaling figure
- `results/` - edge-case outputs and figure
