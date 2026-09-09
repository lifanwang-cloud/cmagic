#!/usr/bin/env python
"""CMAGIC Hubble residuals vs SALT3 c and x1, labeled by numeric snid
(matching the per-object panel filenames in figures/sdss_panels/)."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "sdss_salt3_comparison.csv"))))
D = []
for r in rows:
    if not (r["cm_mu"] and not r["cm_failed"] and r["x1"] and r["c"]):
        continue
    D.append((float(r["cm_mu"]) - float(r["mu_lcdm"]), float(r["c"]), float(r["x1"]),
              r["cm_mode"], r["snid"]))
res = np.array([d[0] for d in D]); res -= np.median(res)
c = np.array([d[1] for d in D]); x1 = np.array([d[2] for d in D])
modes = [d[3] for d in D]; snids = [d[4] for d in D]

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
for ax, x, xl in [(axes[0], c, "SALT3 color $c$"), (axes[1], x1, "SALT3 stretch $x_1$")]:
    for mk, mode, col in [("s", "L", "#c9631a"), ("D", "S", "#2a7a2a")]:
        sel = [i for i, m in enumerate(modes) if m == mode]
        ax.scatter(x[sel], res[sel], marker=mk, s=45, facecolors="none",
                   edgecolors=col, label=f"mode {mode} (n={len(sel)})")
    for k, i in enumerate(np.argsort(x)):
        dx, dy = [(4, 5), (4, -10), (-4, 5), (-4, -10)][k % 4]
        ax.annotate(snids[i], (x[i], res[i]), textcoords="offset points",
                    xytext=(dx, dy), ha="left" if dx > 0 else "right", fontsize=6.5,
                    color={"L": "#c9631a", "S": "#2a7a2a"}[modes[i]])
    r_p = np.corrcoef(x, res)[0, 1]
    slope = np.polyfit(x, res, 1)[0]
    xx = np.linspace(x.min(), x.max(), 5)
    ax.plot(xx, np.polyval(np.polyfit(x, res, 1), xx), "k--", lw=1,
            label=f"r = {r_p:+.2f}; slope {slope:+.2f}")
    ax.axhline(0, color="grey", lw=0.7)
    ax.set_xlabel(xl); ax.legend(fontsize=9, frameon=False); ax.grid(alpha=0.25)
axes[0].set_ylabel(r"CMAGIC $\mu - \mu_{\Lambda \rm CDM}$ (grey removed) [mag]")
fig.suptitle("CMAGIC Hubble residuals vs SALT3 light-curve parameters (labels = snid; panels in figures/sdss_panels/)", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "sdss_cmagic_residual_vs_x1c.png"), dpi=170)
print("wrote figures/sdss_cmagic_residual_vs_x1c.png; n =", len(D))
