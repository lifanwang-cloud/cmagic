#!/usr/bin/env python
"""Hubble-residual diagram, v0.3 standardized: SALT3 + CMAGIC (corrected mu,
inflated errors), CMAGIC points labeled by SDSS numeric snid. Quotes BOTH the
unweighted rms and the error-weighted rms (w = 1/sigma^2) for each method."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
comp = list(csv.DictReader(open(os.path.join(HERE, "sdss_salt3_comparison.csv"))))
stan = list(csv.DictReader(open(os.path.join(HERE, "sdss_standardized.csv"))))

S = [(float(r["z"]), float(r["mu"]) - float(r["mu_lcdm"]), float(r["emu"]))
     for r in comp if r["mu"] and not r["salt3_cut"]]
S = np.array(S); eS = np.sqrt(S[:, 2] ** 2 + 0.10 ** 2)
gS = np.average(S[:, 1], weights=1 / eS ** 2)

C = [(float(r["z"]), float(r["mu_corr"]) - float(r["mu_lcdm"]), float(r["emu_corr"]),
      r["mode"], r["snid"]) for r in stan]
Carr = np.array([(a, b, c) for a, b, c, *_ in C])
modes = [m for *_, m, _ in C]; snids = [s for *_, s in C]
gC = np.average(Carr[:, 1], weights=1 / Carr[:, 2] ** 2)

def stats(res, err):
    rms_u = np.std(res)
    w = 1 / err ** 2
    rms_w = np.sqrt(np.sum(w * res ** 2) / np.sum(w))
    return rms_u, rms_w

sU, sW = stats(S[:, 1] - gS, eS)
cU, cW = stats(Carr[:, 1] - gC, Carr[:, 2])

fig, ax = plt.subplots(figsize=(11, 5.8))
ax.errorbar(S[:, 0], S[:, 1] - gS, yerr=eS, fmt="o", ms=5, color="#1f6fb4",
            ecolor="#9ec4e0", elinewidth=1, capsize=2, zorder=3,
            label=f"SALT3 (n={len(S)}): rms {sU:.3f} unwtd / {sW:.3f} wtd")
for mk, mode, col in [("s", "L", "#c9631a"), ("D", "S", "#2a7a2a")]:
    sel = [i for i, m in enumerate(modes) if m == mode]
    if sel:
        a = Carr[sel]
        u, w = stats(a[:, 1] - gC, a[:, 2])
        ax.errorbar(a[:, 0], a[:, 1] - gC, yerr=a[:, 2], fmt=mk, ms=6, mfc="none",
                    color=col, elinewidth=1, capsize=2, zorder=4,
                    label=f"CMAGIC v0.3 mode {mode} (n={len(sel)}): rms {u:.3f} / {w:.3f} wtd")
for k, i in enumerate(np.argsort(Carr[:, 0])):
    dx, dy = [(4, 6), (4, -11), (-4, 6), (-4, -11)][k % 4]
    ax.annotate(snids[i], (Carr[i, 0], Carr[i, 1] - gC), textcoords="offset points",
                xytext=(dx, dy), ha="left" if dx > 0 else "right", fontsize=6.5,
                color={"L": "#c9631a", "S": "#2a7a2a"}[modes[i]])
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("redshift")
ax.set_ylabel(r"$\mu - \mu_{\Lambda \rm CDM}(H_0{=}72,\ \Omega_m{=}0.3)$  [mag]")
ax.set_title("SDSS-II Hubble residuals, v0.3 standardized chain (weighted grey offsets removed; CMAGIC labeled by SDSS snid)", fontsize=10.5)
ax.legend(fontsize=9, loc="lower left", frameon=False)
ax.set_ylim(-1.1, 1.1); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "sdss_hubble_residuals_v03.png"), dpi=170)
print(f"SALT3 rms unwtd/wtd = {sU:.3f}/{sW:.3f}; CMAGIC all = {cU:.3f}/{cW:.3f}")
