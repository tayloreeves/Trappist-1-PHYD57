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