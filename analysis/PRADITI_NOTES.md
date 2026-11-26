# Notes - Praditi

My part: intro and conclusion of the talk, stability check code, a few edge-case runs, plots.( this is what I make of my role open to changing)

## Data
From Agol et al. 2021. Star mass 0.0898 Msun. Masses, periods and a are in trappist_params.py.
Eccentricities are all below 0.01 and the orbits are almost flat.

## For the intro
- Neighbouring planets are near resonances: 8:5, 5:3, 3:2, 3:2, 4:3, 3:2 (0.2-1.3% off)
- Laplace relations hold for every set of three planets, so the resonances are linked
- Planets are only 6.6-12.8 mutual Hill radii apart (f-g closest)
- Systems this packed usually go unstable, so the resonances are probably what keep it stable
- To read: Tamayo et al. 2017

## Sandbox results
- dt = P_b/100 gives energy error ~2e-8 over 10 yr. Halving dt makes it ~16x smaller on average
- Mass scaling, 8 copies, 1000 yr:
  - x1 and x3: all stable
  - x10: all unstable, median 25 yr
  - x30: all unstable, median 0.06 yr
  - x100: all unstable, median 0.005 yr
  - so it breaks somewhere around 3-4 mutual Hill radii
- Moving planet h off resonance: nothing happens in 1000 yr, need much longer runs
- Speed ~144 yr/s on 2 threads, so long runs need the GPU

## Edge case ideas
- mass scaling (done)
- break a resonance
- equally spaced toy system
- increase one planet's eccentricity
- Lyapunov time / frequency drift

## Questions
- Which parameters do we perturb, by how much, how many copies, how long?
- Random starting phases mess up the resonances. Keep them or also try the observed ones?
- Infall rule: I use 95% of b's starting distance since b starts on the line. OK?
- Roche radius: I used the Hill radius. Does that match L12?
- Results format: copy_id,unstable,t_instability,reason,planet_i,planet_j
- How often to run the check? I use every 4 steps
