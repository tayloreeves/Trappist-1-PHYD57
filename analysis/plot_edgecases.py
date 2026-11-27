# resonance chain + mass scaling figure
# run ./sandbox massscan 1000 8 first, then: python plot_edgecases.py

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from trappist_params import P, names, mu, a

fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))

# left: period ratios vs exact resonances
res = [8/5, 5/3, 3/2, 3/2, 4/3, 3/2]
lab = ["8:5", "5:3", "3:2", "3:2", "4:3", "3:2"]
ratios = P[1:] / P[:-1]
x = np.arange(6)
ax[0].scatter(x, ratios, s=60, color="C0", zorder=3, label="observed")
ax[0].scatter(x, res, marker="_", s=600, color="C3", zorder=2, label="exact resonance")
for i in range(6):
    ax[0].annotate(f"{lab[i]}\n+{100*(ratios[i]/res[i]-1):.2f}%", (x[i], ratios[i]),
                   textcoords="offset points", xytext=(0, 10), ha="center", fontsize=8)
ax[0].set_xticks(x, [f"{names[i+1]}/{names[i]}" for i in range(6)])
ax[0].set_ylabel("period ratio")
ax[0].set_ylim(1.25, 1.8)
ax[0].set_title("TRAPPIST-1 is a resonant chain")
ax[0].legend(loc="lower left", fontsize=8)

# right: instability time vs closest spacing (heavier planets = closer in Hill radii)
d = np.genfromtxt("results/massscan.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")

def delta_min(f):
    mm = mu * f
    RH = ((mm[:-1] + mm[1:]) / 3) ** (1/3) * (a[:-1] + a[1:]) / 2
    return ((a[1:] - a[:-1]) / RH).min()

for f in np.unique(d["mass_factor"]):
    sel = d["mass_factor"] == f
    D = delta_min(f)
    uns = d["unstable"][sel] == 1
    t = d["t_instability"][sel]
    ax[1].scatter(np.full(uns.sum(), D), t[uns], color="C3", s=25)
    if (~uns).any():
        ax[1].scatter(np.full((~uns).sum(), D), t[~uns], color="C2", marker="^", s=40)
    ax[1].annotate(f"mass x{f:g}", (D, t.max()), textcoords="offset points",
                   xytext=(0, 8), ha="center", fontsize=8)

ax[1].scatter([], [], color="C3", label="went unstable")
ax[1].scatter([], [], color="C2", marker="^", label="stable to end of run")
ax[1].set_yscale("log")
ax[1].set_ylim(1e-3, 1e4)
ax[1].set_xlim(1, 7.3)
ax[1].set_xlabel("closest planet spacing [mutual Hill radii]")
ax[1].set_ylabel("instability time [yr]")
ax[1].set_title("Heavier planets -> tighter spacing -> faster chaos")
ax[1].legend(fontsize=8, loc="upper left")

fig.tight_layout()
fig.savefig("results/preface_edgecases.png", dpi=150)
print("saved results/preface_edgecases.png")
