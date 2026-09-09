"""DES-SN5YR: sample standardization + symmetric SALT3 comparison
(cookbook 10.4). Consumes des_validation[_assisted].csv and
des_salt3_comparison.csv — no refitting.

1. CMAGIC standardization (section 6.1): M0 + delta_S + alpha_C x1 + beta_C c
   robust WLS with OUR SALT3 covariates; host_term precedence honored.
2. SALT3 both ways: fixed fiducial Tripp, and (M0, alpha, beta) refit on the
   same objects — the symmetric treatment (cookbook 10.3).
3. Weighted residual table, z-tercile drift, and the two-panel figure.

Run:  python examples/des_standardize_and_fig.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import cmagic                                                # noqa: E402
from des_validation import mu_lcdm, H0, OM                   # noqa: E402

ALPHA_FID, BETA_FID, M_B = 0.14, 3.1, -19.36


def main():
    for name in ('des_validation_assisted.csv', 'des_validation.csv'):
        if os.path.exists(os.path.join(HERE, name)):
            cmrows = list(csv.DictReader(open(os.path.join(HERE, name))))
            print(f'CMAGIC results: {name}')
            break
    salt = {r['snid']: r for r in csv.DictReader(
        open(os.path.join(HERE, 'des_salt3_comparison.csv')))}

    # ---- common fit set: CMAGIC ok + SALT3 ok ----
    rows = []
    for r in cmrows:
        s = salt.get(r['snid'])
        if r['failed'] or not r.get('mu') or s is None or s['salt3_cut'] \
                or not s.get('mu'):
            continue
        rows.append(dict(
            snid=r['snid'], z=float(r['z']), mu=float(r['mu']),
            emu=float(r['emu']) if r['emu'] else 0.15,
            mode=r['mode'], host_term=float(r['host_term'] or 0.),
            x1=float(s['x1']), c=float(s['c']),
            s_mu=float(s['mu']), s_emu=float(s['emu']),
            mu_l=float(r['mu_lcdm'])))
    n = len(rows)
    print(f'common fit set n={n} '
          f'(modes {dict((m, sum(1 for r in rows if r["mode"] == m)) for m in "LQSR")})')

    # ---- 1. CMAGIC standardization ----
    S = cmagic.standardize_sample(rows, covariates=('x1', 'c'),
                                  fit_mode_offsets=True, h0=H0, om=OM)
    print('coefficients:', {k: (round(v[0], 3), round(v[1], 3))
                            for k, v in S.coefficients.items()})
    print('inflation:', {k: round(v, 3) for k, v in S.inflation.items()})
    print('post-fit correlations:', {k: round(v, 2)
                                     for k, v in S.post_correlations.items()})
    z = np.array([r['z'] for r in rows])
    mu_l = np.array([r['mu_l'] for r in rows])
    keep = [i for i in range(n) if S.mu_corr[i] is not None]
    z, mu_l = z[keep], mu_l[keep]
    rows = [rows[i] for i in keep]
    cm_res = np.array([S.mu_corr[i] for i in keep]) - mu_l
    cm_e = np.array([S.emu_corr[i] for i in keep])
    n = len(rows)
    w_cm = 1 / cm_e**2
    cm_res = cm_res - np.sum(w_cm * cm_res) / np.sum(w_cm)
    modes = [r['mode'] for r in rows]
    snids = [r['snid'] for r in rows]

    # ---- 2. SALT3 fixed + refit ----
    s_mu = np.array([r['s_mu'] for r in rows])
    s_emu = np.array([r['s_emu'] for r in rows])
    x1 = np.array([r['x1'] for r in rows])
    c = np.array([r['c'] for r in rows])

    def fit_sigint(resid, e, A=None):
        if A is None:
            A = np.ones((len(resid), 1))
        sig = 0.10
        for _ in range(80):
            w = 1 / (e**2 + sig**2)
            coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None],
                                       resid * np.sqrt(w), rcond=None)
            r_ = resid - A @ coef
            chi2 = np.sum(w * r_**2) / max(len(resid) - A.shape[1], 1)
            if abs(chi2 - 1) < 1e-4:
                break
            sig = max(1e-4, sig * chi2**0.4)
        return coef, r_, sig, 1 / (e**2 + sig**2)

    _, res_fix, sig_fix, w_fix = fit_sigint(s_mu - mu_l, s_emu)
    mB = s_mu - ALPHA_FID * x1 + BETA_FID * c + M_B
    A = np.column_stack([np.ones(n), -x1, c])
    coef, res_fit, sig_fit, w_fit = fit_sigint(mB - mu_l, s_emu, A)
    M0f, alpf, betf = coef
    Cov = np.linalg.inv(A.T @ np.diag(w_fit) @ A)
    e_fit = np.sqrt(s_emu**2 + sig_fit**2)
    print(f'SALT3 refit: alpha={alpf:+.3f}+/-{np.sqrt(Cov[1,1]):.3f}, '
          f'beta={betf:+.3f}+/-{np.sqrt(Cov[2,2]):.3f}, sig_int={sig_fit:.3f} '
          f'(fixed-coef sig_int={sig_fix:.3f})')

    def rms_pair(res, err):
        w = 1 / err**2
        return float(np.std(res)), float(np.sqrt(np.sum(w * res**2) / np.sum(w)))

    sU, sW = rms_pair(res_fit, e_fit)
    cU, cW = rms_pair(cm_res, cm_e)
    print(f'weighted rms: SALT3 refit {sW:.3f}, CMAGIC {cW:.3f} '
          f'(unweighted {sU:.3f} / {cU:.3f})')
    for mode in ('L', 'S'):
        sel = [i for i, m in enumerate(modes) if m == mode]
        if sel:
            u, w = rms_pair(cm_res[sel], cm_e[sel])
            print(f'  CMAGIC mode {mode}: n={len(sel)} rms {u:.3f} / {w:.3f} wtd')
    r_p = float(np.corrcoef(res_fit, cm_res)[0, 1])
    print(f'residual correlation r = {r_p:+.2f} +/- '
          f'{(1 - r_p**2) / np.sqrt(max(n - 3, 1)):.2f}')

    # ---- 3. figure ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (ax, axb) = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2.4, 1]))
    ax.errorbar(z, res_fit, yerr=e_fit, fmt='o', ms=5, color='#1f6fb4',
                ecolor='#9ec4e0', elinewidth=1, capsize=2, zorder=3,
                label=(f'SALT3, coefficients refit on sample '
                       f'($\\alpha$={alpf:+.2f}, $\\beta$={betf:+.2f}; n={n}): '
                       f'rms {sU:.3f} / {sW:.3f} wtd'))
    for mk, mode, col in [('s', 'L', '#c9631a'), ('D', 'S', '#2a7a2a')]:
        sel = [i for i, m in enumerate(modes) if m == mode]
        if sel:
            u, w = rms_pair(cm_res[sel], cm_e[sel])
            ax.errorbar(z[sel], cm_res[sel], yerr=cm_e[sel], fmt=mk, ms=6,
                        mfc='none', color=col, elinewidth=1, capsize=2,
                        zorder=4,
                        label=f'CMAGIC mode {mode} (n={len(sel)}): '
                              f'rms {u:.3f} / {w:.3f} wtd')
    for k, i in enumerate(np.argsort(z)):
        dx, dy = [(4, 6), (4, -11), (-4, 6), (-4, -11)][k % 4]
        ax.annotate(snids[i], (z[i], cm_res[i]), textcoords='offset points',
                    xytext=(dx, dy), ha='left' if dx > 0 else 'right',
                    fontsize=5.5,
                    color={'L': '#c9631a', 'S': '#2a7a2a'}.get(modes[i], 'k'))
    ax.axhline(0, color='k', lw=0.8)
    ax.set_ylabel(r'$\mu - \mu_{\Lambda \rm CDM}(H_0{=}72,\ \Omega_m{=}0.3)$  [mag]')
    ax.set_title('DES-SN5YR Hubble residuals, symmetric treatment '
                 '(both standardizations fit on the same objects)',
                 fontsize=10.5)
    ax.legend(fontsize=8.5, loc='lower left', frameon=False)
    ax.set_ylim(-1.0, 1.0)
    ax.grid(alpha=0.25)

    edges = np.quantile(z, [0, 1 / 3, 2 / 3, 1])
    edges[-1] += 1e-9
    lab = ['low', 'mid', 'high']
    print(f'z tercile edges: {edges.round(3)}')
    print(f'{"bin":5s} {"n":>2s} {"<z>":>6s} {"SALT3 fixed":>14s} '
          f'{"SALT3 refit":>14s} {"CMAGIC":>14s}')
    drift = {}
    for name, res, w, col, mk, off in [
            ('SALT3 fixed fiducial', res_fix - np.sum(w_fix * res_fix) / np.sum(w_fix),
             w_fix, '#8a8a8a', 'o', -0.004),
            ('SALT3 refit', res_fit, w_fit, '#1f6fb4', 'o', 0.0),
            ('CMAGIC standardized', cm_res, w_cm, '#c9631a', 's', +0.004)]:
        ms, es, zc = [], [], []
        for i in range(3):
            m = (z >= edges[i]) & (z < edges[i + 1])
            ms.append(np.sum(w[m] * res[m]) / np.sum(w[m]))
            es.append(1 / np.sqrt(np.sum(w[m])))
            zc.append(z[m].mean())
        drift[name] = (ms[2] - ms[0], float(np.hypot(es[2], es[0])))
        axb.errorbar(np.array(zc) + off, ms, yerr=es, fmt=mk + '-', ms=7,
                     lw=1.2, color=col,
                     mfc='none' if name.startswith('CMAGIC') else col,
                     capsize=3,
                     label=f'{name}: drift {drift[name][0]:+.3f} '
                           f'$\\pm$ {drift[name][1]:.3f}')
    for i in range(3):
        m = (z >= edges[i]) & (z < edges[i + 1])
        cols = []
        for rv, wv in [(res_fix - np.sum(w_fix * res_fix) / np.sum(w_fix), w_fix),
                       (res_fit, w_fit), (cm_res, w_cm)]:
            cols.append(f'{np.sum(wv[m]*rv[m])/np.sum(wv[m]):+.3f}'
                        f'+/-{1/np.sqrt(np.sum(wv[m])):.3f}')
        print(f'{lab[i]:5s} {m.sum():2d} {z[m].mean():6.3f} {cols[0]:>14s} '
              f'{cols[1]:>14s} {cols[2]:>14s}')
    for k, v in drift.items():
        print(f'{k:22s} low->high drift {v[0]:+.3f} +/- {v[1]:.3f} '
              f'({abs(v[0])/v[1]:.1f} sigma)')
    axb.axhline(0, color='k', lw=0.8)
    axb.set_xlabel('redshift')
    axb.set_ylabel('binned weighted\nmean residual [mag]')
    axb.set_title('tercile means: the redshift drift under each treatment',
                  fontsize=9.5)
    axb.legend(fontsize=8, frameon=False, loc='upper left')
    axb.grid(alpha=0.25)
    fig.tight_layout()
    out = os.path.join(HERE, 'figures', 'des_hubble_residuals_symmetric.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=170)
    print('wrote', out)

    # ---- CSV of the standardized set ----
    with open(os.path.join(HERE, 'des_standardized.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['snid', 'z', 'mode', 'x1_ext', 'c_ext', 'mu_corr',
                    'emu_corr', 'mu_lcdm'])
        for i, r in enumerate(rows):
            w.writerow([r['snid'], r['z'], r['mode'], r['x1'], r['c'],
                        round(cm_res[i] + mu_l[i], 4), round(cm_e[i], 4),
                        r['mu_l']])


if __name__ == '__main__':
    main()
