# plots for the main runs: survival curve, instability times, which planets
# python plot_results.py results.csv    (or --fake to test with made-up data)
# columns: copy_id,unstable,t_instability,reason,planet_i,planet_j

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLANETS = ["", "b", "c", "d", "e", "f", "g", "h"]


def make_fake(path="results_fake.csv", K=1000, t_end=1e6, seed=42):
    # fake data just to test the plots
    rng = np.random.default_rng(seed)
    with open(path, "w") as f:
        f.write("copy_id,unstable,t_instability,reason,planet_i,planet_j\n")
        for k in range(K):
            if rng.random() < 0.35:
                t = 10 ** rng.uniform(2, np.log10(t_end))
                reason = rng.choice(["close_encounter", "escape", "infall"], p=[0.8, 0.12, 0.08])
                i = int(rng.integers(1, 7))
                j = i + 1 if reason == "close_encounter" else -1
                f.write(f"{k},1,{t:.6e},{reason},{i},{j}\n")
            else:
                f.write(f"{k},0,{t_end:.6e},stable,-1,-1\n")
    return path


def main(path):
    d = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    K = len(d)
    uns = d["unstable"] == 1
    t = d["t_instability"]
    t_end = t.max()
    n_uns = uns.sum()

    # binomial error on the unstable fraction
    p = n_uns / K
    err = np.sqrt(p * (1 - p) / K)
    print(f"copies: {K}")
    print(f"unstable by {t_end:.3g} yr: {n_uns} = {100*p:.1f} +/- {100*err:.1f} %")
    if n_uns:
        print(f"median instability time: {np.median(t[uns]):.3g} yr")
        for r in ["close_encounter", "escape", "infall"]:
            print(f"  {r}: {(d['reason'][uns] == r).sum()}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    # survival curve with error band
    times = np.logspace(np.log10(max(t[uns].min() if n_uns else 1, 1)), np.log10(t_end), 200)
    surv = np.array([1 - (uns & (t <= tt)).sum() / K for tt in times])
    s_err = np.sqrt(surv * (1 - surv) / K)
    ax[0].plot(times, surv, color="C0")
    ax[0].fill_between(times, surv - s_err, surv + s_err, color="C0", alpha=0.25)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("time [yr]")
    ax[0].set_ylabel("fraction still stable")
    ax[0].set_title("Survival of TRAPPIST-1 copies")

    # when copies go unstable
    if n_uns:
        bins = np.logspace(np.log10(t[uns].min()), np.log10(t_end), 25)
        ax[1].hist(t[uns], bins=bins, color="C1")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("instability time [yr]")
    ax[1].set_ylabel("number of copies")
    ax[1].set_title("When copies go unstable")

    # which planets are involved
    counts = np.zeros(8, dtype=int)
    for i, j in zip(d["planet_i"][uns], d["planet_j"][uns]):
        counts[i] += 1
        if j > 0:
            counts[j] += 1
    ax[2].bar(PLANETS[1:], counts[1:], color="C2")
    ax[2].set_xlabel("planet")
    ax[2].set_ylabel("times involved in instability")
    ax[2].set_title("Which planets cause trouble")

    fig.tight_layout()
    out = path.rsplit(".", 1)[0] + "_plots.png"
    fig.savefig(out, dpi=150)
    print(f"saved {out}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--fake":
        main(make_fake())
    elif len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        print("usage: python plot_results.py results.csv   (or --fake)")
