#!/usr/bin/env python
"""Combined SDSS-II Hubble-residual figure: SALT3 and CMAGIC (modes distinguished),
error bars shown, each method's grey offset removed and quoted, same-object pairs
connected. Reads examples/sdss_salt3_comparison.csv; writes
examples/figures/sdss_hubble_residuals_both.png."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "sdss_salt3_comparison.csv"))))
S, C = [], []
for r in rows:
    try:
        z = float(r["z"]); mul = float(r["mu_lcdm"])
    except ValueError:
        continue
    if r["mu"] and not r["salt3_cut"]:
        S.append((z, float(r["mu"]) - mul, float(r["emu"])))
    if r["cm_mu"] and not r["cm_failed"]:
        C.append((z, float(r["cm_mu"]) - mul, float(r["cm_emu"]), r["cm_mode"]))
S = np.array(S, dtype=float)
Carr = np.array([(a, b, c) for a, b, c, _ in C], dtype=float)
modes = [m for *_, m in C]
gS = np.median(S[:, 1]); gC = np.median(Carr[:, 1])
eS = np.sqrt(S[:, 2] ** 2 + 0.10 ** 2)

fig, ax = plt.subplots(figsize=(10, 5.5))
zc = {round(z, 6): dm - gC for (z, dm, e), m in zip(Carr, modes)}
for z, dm, e in S:
    if round(z, 6) in zc:
        ax.plot([z, z], [dm - gS, zc[round(z, 6)]], color="lightgrey", lw=0.8, zorder=1)
ax.errorbar(S[:, 0], S[:, 1] - gS, yerr=eS, fmt="o", ms=5, color="#1f6fb4",
            ecolor="#9ec4e0", elinewidth=1, capsize=2, zorder=3,
            label=f"SALT3 (n={len(S)}; grey {gS:+.2f} removed; err incl. $\\sigma_{{int}}$=0.10)")
for mk, mode, col in [("s", "L", "#c9631a"), ("D", "S", "#2a7a2a")]:
    sel = [i for i, m in enumerate(modes) if m == mode]
    if sel:
        a = Carr[sel]
        ax.errorbar(a[:, 0], a[:, 1] - gC, yerr=a[:, 2], fmt=mk, ms=6, mfc="none",
                    color=col, elinewidth=1, capsize=2, zorder=4,
                    label=f"CMAGIC mode {mode} (n={len(sel)})")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("redshift")
ax.set_ylabel(r"$\mu - \mu_{\Lambda \rm CDM}(H_0{=}72,\ \Omega_m{=}0.3)$  [mag]")
ax.set_title(f"SDSS-II Hubble residuals: SALT3 vs CMAGIC (same photometry; CMAGIC grey {gC:+.2f} removed)", fontsize=11)
ax.legend(fontsize=9, loc="lower left", frameon=False)
ax.set_ylim(-1.6, 1.6); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "sdss_hubble_residuals_both.png"), dpi=170)
print("wrote figures/sdss_hubble_residuals_both.png")
