# TRAPPIST-1 data (Agol et al. 2021) + resonance and Hill spacing numbers
# python trappist_params.py

import numpy as np
from fractions import Fraction

M_EARTH_SUN = 3.0034896e-6
M_STAR, M_STAR_ERR = 0.0898, 0.0023   # solar masses

names = ["b", "c", "d", "e", "f", "g", "h"]
P     = np.array([1.510826, 2.421937, 4.049219, 6.101013, 9.207540, 12.352446, 18.772866])  # days
P_err = np.array([0.000006, 0.000018, 0.000026, 0.000035, 0.000032, 0.000054, 0.000214])
m     = np.array([1.374, 1.308, 0.388, 0.692, 1.039, 1.321, 0.326])   # Earth masses
m_err = np.array([0.069, 0.056, 0.012, 0.022, 0.031, 0.038, 0.020])
a     = np.array([0.01154, 0.01580, 0.02227, 0.02925, 0.03849, 0.04683, 0.06189])  # AU
a_err = np.array([0.00010, 0.00013, 0.00019, 0.00025, 0.00033, 0.00040, 0.00053])
ecosw = np.array([-0.00215, 0.00055, -0.00496, 0.00433, -0.00840, 0.00380, -0.00365])
esinw = np.array([0.00217, 0.00001, 0.00267, -0.00461, -0.00051, 0.00128, -0.00002])
ecc = np.hypot(ecosw, esinw)

mu = m * M_EARTH_SUN / M_STAR   # planet/star mass ratio

if __name__ == "__main__":
    print("planet  P[d]      m[Me]  a[AU]    e       m/M*")
    for i in range(7):
        print(f"  {names[i]}    {P[i]:8.4f}  {m[i]:5.3f}  {a[i]:.5f}  {ecc[i]:.4f}  {mu[i]:.3e}")

    # nearest simple ratio p:q for each neighbouring pair
    print("\nperiod ratios:")
    for i in range(6):
        r = P[i+1] / P[i]
        frac = Fraction(r).limit_denominator(5)
        off = (r / float(frac) - 1) * 100
        print(f"  {names[i+1]}/{names[i]} = {r:.4f}  ~ {frac.numerator}:{frac.denominator}  ({off:+.2f}%)")

    # three-body Laplace relations: p/P1 - (p+q)/P2 + q/P3 ~ 0
    print("\nLaplace relations:")
    pq = [(2, 3), (1, 2), (2, 3), (1, 2), (1, 1)]
    for i, (p_, q_) in enumerate(pq):
        val = p_/P[i] - (p_+q_)/P[i+1] + q_/P[i+2]
        print(f"  {names[i]}{names[i+1]}{names[i+2]} (p={p_}, q={q_}): {val:+.2e} /day  (relative {val*P[i]:+.1e})")

    # spacing in mutual Hill radii
    print("\nspacing in mutual Hill radii:")
    for i in range(6):
        RH = ((mu[i] + mu[i+1]) / 3) ** (1/3) * (a[i] + a[i+1]) / 2
        print(f"  {names[i]}-{names[i+1]}: {(a[i+1] - a[i]) / RH:5.1f}")

    # single-planet Hill radius (our close encounter limit is 2x this)
    print("\nHill radius:")
    for i in range(7):
        rH = a[i] * (mu[i] / 3) ** (1/3)
        print(f"  {names[i]}: {rH:.2e} AU ({100*rH/a[i]:.2f}% of a)")

    T0 = P[0] / (2*np.pi)
    print(f"\ncode units: 1 time unit = {T0:.5f} d, 1 yr = {365.25/T0:.1f} time units")
