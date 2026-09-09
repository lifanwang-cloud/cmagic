#!/usr/bin/env python
"""Hubble-residual diagram under the SYMMETRIC treatment (cookbook 10.3, 5.1):
SALT3 with its Tripp coefficients REFIT on the same 28-object sample that
CMAGIC's standardization is fit on -- the apples-to-apples version of the
v0.3 figure. Top: per-object residuals (CMAGIC labeled by SDSS snid).
Bottom: tercile-binned weighted means for the three treatments -- SALT3 with
fixed fiducial coefficients, SALT3 refit, CMAGIC standardized -- the redshift
drift decomposition of cookbook 10.3 made visible.

Run: python examples/make_hubble_residual_fig_symmetric.py
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA_FID, BETA_FID, M_B = 0.14, 3.1, -19.36

salt = {r["snid"]: r for r in csv.DictReader(open(os.path.join(HERE, "sdss_salt3_comparison.csv")))}
stan = list(csv.DictReader(open(os.path.join(HERE, "sdss_standardized.csv"))))

snids = [r["snid"] for r in stan]
z = np.array([float(r["z"]) for r in stan])
mu_l = np.array([float(r["mu_lcdm"]) for r in stan])
modes = [r["mode"] for r in stan]
cm_res = np.array([float(r["mu_corr"]) for r in stan]) - mu_l
cm_e = np.array([float(r["emu_corr"]) for r in stan])
s_mu = np.array([float(salt[s]["mu"]) for s in snids])
s_emu = np.array([float(salt[s]["emu"]) for s in snids])
x1 = np.array([float(salt[s]["x1"]) for s in snids])
c = np.array([float(salt[s]["c"]) for s in snids])
n = len(snids)


def fit_sigint(resid, e, A=None):
    """Robust WLS with intrinsic term iterated to chi2/dof = 1 (WS06 convention)."""
    if A is None:
        A = np.ones((len(resid), 1))
    sig = 0.10
    for _ in range(80):
        w = 1 / (e**2 + sig**2)
        coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], resid * np.sqrt(w), rcond=None)
        r = resid - A @ coef
        chi2 = np.sum(w * r**2) / max(len(resid) - A.shape[1], 1)
        if abs(chi2 - 1) < 1e-4:
            break
        sig = max(1e-4, sig * chi2**0.4)
    return coef, r, sig, 1 / (e**2 + sig**2)


# SALT3, fixed fiducial coefficients (grey-offset removed)
_, res_fix, sig_fix, w_fix = fit_sigint(s_mu - mu_l, s_emu)

# SALT3, Tripp coefficients refit on this sample (the CMAGIC treatment)
mB = s_mu - ALPHA_FID * x1 + BETA_FID * c + M_B
A = np.column_stack([np.ones(n), -x1, c])
coef, res_fit, sig_fit, w_fit = fit_sigint(mB - mu_l, s_emu, A)
M0f, alpf, betf = coef
Cov = np.linalg.inv(A.T @ np.diag(w_fit) @ A)
e_fit = np.sqrt(s_emu**2 + sig_fit**2)

# CMAGIC standardized (errors already inflated by the v0.3.1 chain)
w_cm = 1 / cm_e**2
res_cm = cm_res - np.sum(w_cm * cm_res) / np.sum(w_cm)


def rms_pair(res, err):
    w = 1 / err**2
    return float(np.std(res)), float(np.sqrt(np.sum(w * res**2) / np.sum(w)))


sU, sW = rms_pair(res_fit, e_fit)
cU, cW = rms_pair(res_cm, cm_e)

fig, (ax, axb) = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True,
                              gridspec_kw=dict(height_ratios=[2.4, 1]))
ax.errorbar(z, res_fit, yerr=e_fit, fmt="o", ms=5, color="#1f6fb4",
            ecolor="#9ec4e0", elinewidth=1, capsize=2, zorder=3,
            label=(f"SALT3, coefficients refit on sample "
                   f"($\\alpha$={alpf:+.2f}, $\\beta$={betf:+.2f}; n={n}): "
                   f"rms {sU:.3f} / {sW:.3f} wtd"))
for mk, mode, col in [("s", "L", "#c9631a"), ("D", "S", "#2a7a2a")]:
    sel = [i for i, m in enumerate(modes) if m == mode]
    if sel:
        u, w = rms_pair(res_cm[sel], cm_e[sel])
        ax.errorbar(z[sel], res_cm[sel], yerr=cm_e[sel], fmt=mk, ms=6, mfc="none",
                    color=col, elinewidth=1, capsize=2, zorder=4,
                    label=f"CMAGIC v0.3.1 mode {mode} (n={len(sel)}): rms {u:.3f} / {w:.3f} wtd")
for k, i in enumerate(np.argsort(z)):
    dx, dy = [(4, 6), (4, -11), (-4, 6), (-4, -11)][k % 4]
    ax.annotate(snids[i], (z[i], res_cm[i]), textcoords="offset points",
                xytext=(dx, dy), ha="left" if dx > 0 else "right", fontsize=6.5,
                color={"L": "#c9631a", "S": "#2a7a2a"}[modes[i]])
ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel(r"$\mu - \mu_{\Lambda \rm CDM}(H_0{=}72,\ \Omega_m{=}0.3)$  [mag]")
ax.set_title("SDSS-II Hubble residuals, symmetric treatment: both methods' "
             "standardization fit on the same 28 objects", fontsize=10.5)
ax.legend(fontsize=8.5, loc="lower left", frameon=False)
ax.set_ylim(-0.85, 0.85)
ax.grid(alpha=0.25)

# bottom: tercile-binned weighted means, three treatments
edges = np.quantile(z, [0, 1 / 3, 2 / 3, 1])
edges[-1] += 1e-9
print(f"z tercile edges: {edges.round(3)}   "
      f"SALT3 refit: alpha={alpf:+.3f}+/-{np.sqrt(Cov[1,1]):.3f}, "
      f"beta={betf:+.3f}+/-{np.sqrt(Cov[2,2]):.3f}, sig_int={sig_fit:.3f}")
print(f'{"bin":5s} {"n":>2s} {"SALT3 fixed":>14s} {"SALT3 refit":>14s} {"CMAGIC":>14s}')
drift = {}
for name, res, w, col, mk, off in [
        ("SALT3 fixed fiducial", res_fix, w_fix, "#8a8a8a", "o", -0.004),
        ("SALT3 refit", res_fit, w_fit, "#1f6fb4", "o", 0.0),
        ("CMAGIC standardized", res_cm, w_cm, "#c9631a", "s", +0.004)]:
    ms, es, zc = [], [], []
    for i in range(3):
        m = (z >= edges[i]) & (z < edges[i + 1])
        ms.append(np.sum(w[m] * res[m]) / np.sum(w[m]))
        es.append(1 / np.sqrt(np.sum(w[m])))
        zc.append(z[m].mean())
    drift[name] = (ms[2] - ms[0], float(np.hypot(es[2], es[0])))
    axb.errorbar(np.array(zc) + off, ms, yerr=es, fmt=mk + "-", ms=7, lw=1.2,
                 color=col, mfc="none" if name.startswith("CMAGIC") else col,
                 capsize=3, label=f"{name}: drift {drift[name][0]:+.3f} $\\pm$ {drift[name][1]:.3f}")
lab = ["low", "mid", "high"]
for i in range(3):
    m = (z >= edges[i]) & (z < edges[i + 1])
    row = [f"{np.sum(wv[m]*rv[m])/np.sum(wv[m]):+.3f}+/-{1/np.sqrt(np.sum(wv[m])):.3f}"
           for rv, wv in [(res_fix, w_fix), (res_fit, w_fit), (res_cm, w_cm)]]
    print(f"{lab[i]:5s} {m.sum():2d} {row[0]:>14s} {row[1]:>14s} {row[2]:>14s}")
for k, v in drift.items():
    print(f"{k:22s} low->high drift {v[0]:+.3f} +/- {v[1]:.3f} ({abs(v[0])/v[1]:.1f} sigma)")
axb.axhline(0, color="k", lw=0.8)
axb.set_xlabel("redshift")
axb.set_ylabel("binned weighted\nmean residual [mag]")
axb.set_title("tercile means: the redshift drift under each treatment (cookbook 10.3)",
              fontsize=9.5)
axb.legend(fontsize=8, frameon=False, loc="upper left")
axb.grid(alpha=0.25)
fig.tight_layout()
out = os.path.join(HERE, "figures", "sdss_hubble_residuals_symmetric.png")
fig.savefig(out, dpi=170)
print("wrote", out)
print(f"weighted rms: SALT3 refit {sW:.3f}, CMAGIC {cW:.3f} "
      f"(unweighted {sU:.3f} / {cU:.3f})")
