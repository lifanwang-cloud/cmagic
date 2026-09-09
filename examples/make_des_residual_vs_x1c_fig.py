#!/usr/bin/env python
"""DES-SN5YR: standardized CMAGIC Hubble residuals vs SALT3 c and x1,
every point labeled by its DES SNID (cookbook 6.1 acceptance criterion:
post-fit correlations consistent with zero). Uses the clean n=52 set from
des_standardize_and_fig.py (des_standardized.csv) and covariate errors from
des_salt3_comparison.csv.
Writes figures/des_cmagic_residual_vs_x1c.png."""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
std = list(csv.DictReader(open(os.path.join(HERE, "des_standardized.csv"))))
salt = {r["snid"]: r for r in csv.DictReader(
    open(os.path.join(HERE, "des_salt3_comparison.csv")))}

snids = [r["snid"] for r in std]
modes = [r["mode"] for r in std]
res = np.array([float(r["mu_corr"]) - float(r["mu_lcdm"]) for r in std])
emu = np.array([float(r["emu_corr"]) for r in std])
x1 = np.array([float(r["x1_ext"]) for r in std])
c = np.array([float(r["c_ext"]) for r in std])
ex1 = np.array([float(salt[s].get("ex1") or 0.) for s in snids])
ec = np.array([float(salt[s].get("ec") or 0.) for s in snids])
w = 1 / emu**2
res = res - np.sum(w * res) / np.sum(w)

fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.0), sharey=True)
COL = {"L": "#c9631a", "S": "#2a7a2a"}
for ax, x, xe, xl in [(axes[0], c, ec, "SALT3 color $c$"),
                      (axes[1], x1, ex1, "SALT3 stretch $x_1$")]:
    for mk, mode in [("s", "L"), ("D", "S")]:
        sel = [i for i, m in enumerate(modes) if m == mode]
        ax.errorbar(x[sel], res[sel], yerr=emu[sel],
                    xerr=xe[sel] if np.any(xe[sel] > 0) else None,
                    fmt=mk, ms=7, mfc="none", mec=COL[mode], ecolor=COL[mode],
                    elinewidth=0.6, ls="none",
                    label=f"mode {mode} (n={len(sel)})")
    for k, i in enumerate(np.argsort(x)):
        dx, dy = [(4, 6), (4, -11), (-4, 6), (-4, -11)][k % 4]
        ax.annotate(snids[i], (x[i], res[i]), textcoords="offset points",
                    xytext=(dx, dy), ha="left" if dx > 0 else "right",
                    fontsize=5.5, color=COL.get(modes[i], "k"))
    r_w = float(np.corrcoef(x, res)[0, 1])
    # weighted correlation (the fitted quantity)
    xm = np.sum(w * x) / np.sum(w); rm = np.sum(w * res) / np.sum(w)
    r_wt = float(np.sum(w * (x - xm) * (res - rm)) /
                 np.sqrt(np.sum(w * (x - xm)**2) * np.sum(w * (res - rm)**2)))
    xx = np.linspace(x.min(), x.max(), 5)
    ax.plot(xx, np.polyval(np.polyfit(x, res, 1, w=w), xx), "k--", lw=1,
            label=f"weighted r = {r_wt:+.2f} (unwtd {r_w:+.2f})")
    ax.axhline(0, color="grey", lw=0.7)
    ax.set_xlabel(xl)
    ax.legend(fontsize=9, frameon=False)
    ax.grid(alpha=0.25)
axes[0].set_ylabel(r"standardized $\mu - \mu_{\Lambda \rm CDM}$ [mag]")
fig.suptitle("DES-SN5YR: CMAGIC residuals vs SALT3 parameters AFTER the "
             "fitted corrections (cookbook 6.1); clean n=52, labels = DES SNID",
             fontsize=11)
fig.tight_layout()
out = os.path.join(HERE, "figures", "des_cmagic_residual_vs_x1c.png")
fig.savefig(out, dpi=170)
print("wrote", out)
