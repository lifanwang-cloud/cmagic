"""Decompose the SALT3 redshift drift on the SDSS-II subset (cookbook 10.3).

The v0.3.1 residual table showed SALT3's weighted Hubble residual drifting
+0.161 +/- 0.049 mag (3.3 sigma) from the low to the high redshift bin while
CMAGIC stayed flat. This script tests whether the drift is SALT3 or the
comparison's construction: CMAGIC's standardization coefficients were FITTED
on this sample (section 6.1) while SALT3 used FIXED fiducial Tripp
coefficients (alpha = 0.14, beta = 3.1). In a magnitude-limited survey the
covariates x1, c drift with z, so any coefficient mismatch maps that drift
into a residual-vs-z slope.

Test: refit (M0, alpha, beta) for SALT3 on the same 28 objects with the same
robust-WLS treatment, and recompute the z-binned weighted residuals.

Inputs: sdss_salt3_comparison.csv, sdss_standardized.csv (both produced by
sdss_salt3_comparison.py and sdss_standardize.py).
Run:  python examples/sdss_salt3_drift_diagnostic.py
"""
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA_FID, BETA_FID, M_B = 0.14, 3.1, -19.36


def main():
    salt = {r['snid']: r for r in csv.DictReader(
        open(os.path.join(HERE, 'sdss_salt3_comparison.csv')))}
    std = list(csv.DictReader(
        open(os.path.join(HERE, 'sdss_standardized.csv'))))

    snids = [r['snid'] for r in std]
    z = np.array([float(r['z']) for r in std])
    mu_l = np.array([float(r['mu_lcdm']) for r in std])
    cm_mu = np.array([float(r['mu_corr']) for r in std])
    cm_emu = np.array([float(r['emu_corr']) for r in std])
    s_mu = np.array([float(salt[s]['mu']) for s in snids])
    s_emu = np.array([float(salt[s]['emu']) for s in snids])
    x1 = np.array([float(salt[s]['x1']) for s in snids])
    c = np.array([float(salt[s]['c']) for s in snids])
    n = len(snids)
    print(f'common fit set n={n}')

    def werr(resid, e, sig_int):
        w = 1 / (e**2 + sig_int**2)
        m = np.sum(w * resid) / np.sum(w)
        return m, 1 / np.sqrt(np.sum(w)), w

    def solve_sigint(resid, e):
        # intrinsic term iterated until chi2/dof = 1 (WS06 convention)
        sig = 0.10
        for _ in range(60):
            m, _, w = werr(resid, e, sig)
            chi2 = np.sum(w * (resid - m)**2) / (len(resid) - 1)
            if abs(chi2 - 1) < 1e-4:
                break
            sig = max(1e-4, sig * chi2**0.4)
        return sig

    # 1. SALT3 with the FIXED fiducial coefficients (as shipped)
    res_fix = s_mu - mu_l
    sig_fix = solve_sigint(res_fix, s_emu)
    m_fix, _, w_fix = werr(res_fix, s_emu, sig_fix)
    res_fix = res_fix - m_fix

    # 2. SALT3 with coefficients FITTED on this sample (the CMAGIC treatment)
    #    mu = mB + alpha*x1 - beta*c - M_B  =>  reconstruct mB from shipped mu
    mB = s_mu - ALPHA_FID * x1 + BETA_FID * c + M_B
    y = mB - mu_l
    A = np.column_stack([np.ones(n), -x1, c])   # y = M0 - alpha*x1 + beta*c
    sig = 0.10
    for _ in range(80):
        w = 1 / (s_emu**2 + sig**2)
        coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w),
                                   rcond=None)
        r = y - A @ coef
        chi2 = np.sum(w * r**2) / (n - 3)
        if abs(chi2 - 1) < 1e-4:
            break
        sig = max(1e-4, sig * chi2**0.4)
    M0f, alpf, betf = coef
    Cov = np.linalg.inv(A.T @ np.diag(1 / (s_emu**2 + sig**2)) @ A)
    ealp, ebet = np.sqrt(Cov[1, 1]), np.sqrt(Cov[2, 2])
    res_fit = y - A @ coef
    w_sfit = 1 / (s_emu**2 + sig**2)
    print(f'SALT3 fitted on sample: alpha={alpf:+.3f}+/-{ealp:.3f} '
          f'(fid {ALPHA_FID}), beta={betf:+.3f}+/-{ebet:.3f} (fid {BETA_FID}), '
          f'sig_int={sig:.3f} (fixed-coef sig_int={sig_fix:.3f})')

    # 3. CMAGIC standardized as shipped (errors already inflated)
    res_cm = cm_mu - mu_l
    m_cm, _, w_cm = werr(res_cm, cm_emu, 0.0)
    res_cm = res_cm - m_cm

    # z terciles of the common set
    edges = np.quantile(z, [0, 1 / 3, 2 / 3, 1])
    edges[-1] += 1e-9
    lab = ['low', 'mid', 'high']
    print(f'\nz tercile edges: {edges.round(3)}')
    print(f'{"bin":5s} {"n":>2s} {"<z>":>6s} {"<x1>":>6s} {"<c>":>7s} '
          f'{"SALT3 fixed":>12s} {"SALT3 fitted":>13s} {"CMAGIC":>14s}')
    rows = []
    for i in range(3):
        m = (z >= edges[i]) & (z < edges[i + 1])

        def bw(res, w):
            mm = np.sum(w[m] * res[m]) / np.sum(w[m])
            return mm, 1 / np.sqrt(np.sum(w[m]))
        f0, ef0 = bw(res_fix, w_fix)
        f1, ef1 = bw(res_fit, w_sfit)
        f2, ef2 = bw(res_cm, w_cm)
        rows.append((f0, ef0, f1, ef1, f2, ef2))
        print(f'{lab[i]:5s} {m.sum():2d} {z[m].mean():6.3f} '
              f'{x1[m].mean():+6.2f} {c[m].mean():+7.3f} '
              f'{f0:+.3f}+/-{ef0:.3f} {f1:+.3f}+/-{ef1:.3f} '
              f'  {f2:+.3f}+/-{ef2:.3f}')

    for name, lo, hi in [('SALT3 fixed ', rows[0][:2], rows[2][:2]),
                         ('SALT3 fitted', rows[0][2:4], rows[2][2:4]),
                         ('CMAGIC      ', rows[0][4:6], rows[2][4:6])]:
        d, ed = hi[0] - lo[0], float(np.hypot(hi[1], lo[1]))
        print(f'{name} low->high drift: {d:+.3f} +/- {ed:.3f}  '
              f'({abs(d) / ed:.1f} sigma)')

    # decomposition: coefficient mismatch x covariate drift
    mlo = (z >= edges[0]) & (z < edges[1])
    mhi = (z >= edges[2]) & (z < edges[3])
    dx1 = x1[mhi].mean() - x1[mlo].mean()
    dc = c[mhi].mean() - c[mlo].mean()
    pred = (ALPHA_FID - alpf) * dx1 - (BETA_FID - betf) * dc
    print(f'\ncovariate drift low->high: d<x1>={dx1:+.2f}, d<c>={dc:+.3f}')
    print(f'predicted drift from (fiducial - fitted) coefficient mismatch: '
          f'{pred:+.3f} mag')

    for name, res, w in [('SALT3 fixed ', res_fix, w_fix),
                         ('SALT3 fitted', res_fit, w_sfit),
                         ('CMAGIC      ', res_cm, w_cm)]:
        Az = np.column_stack([np.ones(n), z - z.mean()])
        cf, *_ = np.linalg.lstsq(Az * np.sqrt(w)[:, None], res * np.sqrt(w),
                                 rcond=None)
        Cv = np.linalg.inv(Az.T @ np.diag(w) @ Az)
        print(f'{name} slope vs z: {cf[1]:+.2f} +/- {np.sqrt(Cv[1, 1]):.2f} '
              f'mag per unit z')

    for nm, v in [('x1', x1), ('c', c)]:
        print(f'r(SALT3 fixed resid, {nm}) = '
              f'{np.corrcoef(res_fix, v)[0, 1]:+.2f}   '
              f'r({nm}, z) = {np.corrcoef(v, z)[0, 1]:+.2f}')


if __name__ == '__main__':
    main()
