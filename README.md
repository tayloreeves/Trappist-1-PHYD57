# Trappist-1-PHYD57
PHYD57 A4 Group Assignment


This can serve as the main README doc for the whole repo


## integrator (rk4)

rk4 n-body integrator, newtonian gravity. units: au, yr, msun (g = 4 pi^2).

state layout: x[nd*i+k] = coordinate k of body i (same for v, a).
nd = 2 by default, compile with -Dnd=3 for 3d.

functions (src/integrator.h):
- accel(n,x,m,a)           acceleration of every body
- energy(n,x,v,m)          total energy, for accuracy checks
- rk4_step(n,x,v,m,dt,w)   one step in place. w = scratch of wsize(n) doubles,
                           one per thread, caller-owned (no allocation inside)

build and test: gcc -O3 -Wall src/integrator.c tests/test_integrator.c -o test_integrator -lm

integrator.h is the shared interface: changes go through pull requests.


# analysis (Praditi)

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

Analysis edited (Praditi)
analysis/stability.hpp, same layout as the integrator (x[nd*i+k]), 2D or 3D.

- stab_setup(x,v,m,R,nd,G)    once per copy, from the starting state
- stab_check(x,m,L)           every few steps

stops the run on: collision, tidal disruption (roche limit), close encounter
(2 roche lobes, assignment), escape (2x planet h), infall (b's periapsis - 2 hill radii).
hill sphere entries are only flagged (hill_flags), the run continues.

results format: copy_id,unstable,t_instability,reason,planet_i,planet_j
more details: analysis/README.md
