"""v0.3 sample-level standardization of the SDSS-II validation set
(cookbook sections 6.1 and 10.1): rerun the blind chain (with the v0.3 boundary
gates), fit the sample-level M0 + delta_S + alpha_C x1 + beta_C c model with
external SALT3 covariates (from sdss_salt3_comparison.csv), and report/plot the
pre/post comparison. The PI's acceptance criterion: post-fit residual
correlations vs x1 and c consistent with zero.

Outputs: sdss_standardized.csv, figures/sdss_cmagic_residual_vs_x1c_post.png,
figures/sdss_hubble_pre_post.png.
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from sdss_validation import load_cache, mu_lcdm, H0, OM      # noqa: E402
import cmagic                                                # noqa: E402


def main():
    meta, phot = load_cache()
    covs, coverr = {}, {}
    for r in csv.DictReader(open(os.path.join(HERE, 'sdss_salt3_comparison.csv'))):
        if r.get('x1') and r.get('c') and not r.get('salt3_cut'):
            covs[r['snid']] = (float(r['x1']), float(r['c']))
            coverr[r['snid']] = (float(r.get('ex1') or 0.),
                                 float(r.get('ec') or 0.))
    results, ids, newly_gated = [], [], []
    FILTERS = {'g': 'sdssg', 'r': 'sdssr', 'i': 'sdssi', 'z': 'sdssz'}
    prev = {r['snid']: r for r in csv.DictReader(
        open(os.path.join(HERE, 'sdss_validation.csv')))}
    for m in meta:
        rows = phot[m['snid']]
        table = dict(mjd=np.array([float(r['mjd']) for r in rows]),
                     band=np.array([r['band'] for r in rows]),
                     mag=np.array([float(r['mag']) for r in rows]),
                     emag=np.array([float(r['emag']) for r in rows]))
        try:
            res = cmagic.distance(photometry=table, z=m['z'], filters=FILTERS,
                                  engine='auto', ebv_mw=m['mwebv'], n_mc=120,
                                  seed=42, h0=H0)
        except Exception as e:
            print(f"  {m['snid']}: exception {str(e)[:50]}")
            continue
        res.provenance['snid'] = m['snid']
        results.append(res); ids.append(m['snid'])
        was = prev.get(m['snid'], {}).get('failed', '')
        if res.failed and not was:
            newly_gated.append((m['snid'], res.failed,
                                ';'.join(res.flags)[:60]))
    okmask = [r.failed is None for r in results]
    print(f'rerun: {sum(okmask)} of {len(results)} fit '
          f'(v0.2 had 36); newly gated: {newly_gated}')
    ext = [covs.get(s) for s in ids]
    fitset = [r for r, s in zip(results, ids)
              if r.failed is None and covs.get(s)]
    fitext = [covs[s] for r, s in zip(results, ids)
              if r.failed is None and covs.get(s)]
    fitids = [s for r, s in zip(results, ids)
              if r.failed is None and covs.get(s)]
    S = cmagic.standardize_sample(fitset, covariates=('x1', 'c'),
                                  fit_mode_offsets=True, h0=H0, om=OM,
                                  external_covariates=fitext)
    print('coefficients:', S.coefficients)
    print('post-fit correlations:', S.post_correlations)
    print('inflation factors:', S.inflation)
    print('pre :', {k: v for k, v in S.stats_pre.items()})
    print('post:', {k: v for k, v in S.stats_post.items()})
    # ---- CSV ----
    with open(os.path.join(HERE, 'sdss_standardized.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['snid', 'z', 'mode', 'x1_ext', 'c_ext', 'mu_raw',
                    'mu_corr', 'emu_corr', 'mu_lcdm'])
        for i, (r, s) in enumerate(zip(fitset, fitids)):
            zz = r.provenance.get('true_z')
            w.writerow([s, zz, r.mode, fitext[i][0], fitext[i][1], r.mu,
                        S.mu_corr[i], S.emu_corr[i],
                        round(float(mu_lcdm(zz)), 4)])
    # ---- figures ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    zs = np.array([r.provenance['true_z'] for r in fitset])
    modes = [r.mode for r in fitset]
    x1 = np.array([e[0] for e in fitext]); c = np.array([e[1] for e in fitext])
    res_post = np.array([S.mu_corr[i] - mu_lcdm(zs[i])
                         for i in range(len(fitset))])
    res_pre = np.array([r.mu - mu_lcdm(z) for r, z in zip(fitset, zs)])
    res_pre -= np.median(res_pre)
    # residual vs x1/c, POST-correction, same style/labels as the pre figure
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
    emu_c = np.array([S.emu_corr[i] or 0.15 for i in range(len(fitset))])
    ex1 = np.array([coverr.get(s_, (0., 0.))[0] for s_ in fitids])
    ecc = np.array([coverr.get(s_, (0., 0.))[1] for s_ in fitids])
    for ax, x, xe, xl in [(axes[0], c, ecc, 'SALT3 color $c$'),
                          (axes[1], x1, ex1, 'SALT3 stretch $x_1$')]:
        for mk, mode, col in [('s', 'L', '#c9631a'), ('D', 'S', '#2a7a2a')]:
            sel = [i for i, m in enumerate(modes) if m == mode]
            ax.errorbar(x[sel], res_post[sel], yerr=emu_c[sel],
                        xerr=xe[sel] if np.any(xe[sel] > 0) else None,
                        fmt=mk, ms=7, mfc='none', mec=col, ecolor=col,
                        elinewidth=0.6, ls='none',
                        label=f'mode {mode} (n={len(sel)})')
        for k, i in enumerate(np.argsort(x)):
            dx, dy = [(4, 5), (4, -10), (-4, 5), (-4, -10)][k % 4]
            ax.annotate(fitids[i], (x[i], res_post[i]),
                        textcoords='offset points', xytext=(dx, dy),
                        ha='left' if dx > 0 else 'right', fontsize=6.5,
                        color={'L': '#c9631a', 'S': '#2a7a2a'}[modes[i]])
        r_p = np.corrcoef(x, res_post)[0, 1]
        xx = np.linspace(x.min(), x.max(), 5)
        ax.plot(xx, np.polyval(np.polyfit(x, res_post, 1), xx), 'k--', lw=1,
                label=f'r = {r_p:+.2f} (post-fit)')
        ax.axhline(0, color='grey', lw=0.7)
        ax.set_xlabel(xl); ax.legend(fontsize=9, frameon=False)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r'corrected $\mu - \mu_{\Lambda \rm CDM}$ [mag]')
    fig.suptitle('CMAGIC residuals vs SALT3 parameters AFTER the v0.3 fitted '
                 'corrections (cookbook 6.1); labels = snid', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'figures',
                             'sdss_cmagic_residual_vs_x1c_post.png'), dpi=170)
    # Hubble residuals pre/post
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)
    emu_pre = np.array([(5 / np.log(10)) * r.eD_mpc / r.D_mpc
                        if (r.eD_mpc and r.D_mpc) else 0.15 for r in fitset])
    for a, rr, ee, ttl in ((a1, res_pre, emu_pre, 'raw chain (grey removed)'),
                           (a2, res_post, emu_c, 'standardized '
                            '(M0, delta_S, alpha_C, beta_C removed)')):
        for mk, mode, col in [('s', 'L', '#c9631a'), ('D', 'S', '#2a7a2a')]:
            sel = [i for i, m in enumerate(modes) if m == mode]
            a.errorbar(zs[sel], rr[sel], yerr=ee[sel], fmt=mk, ms=6, mfc='none',
                       mec=col, ecolor=col, elinewidth=0.6, ls='none',
                       label=f'mode {mode}')
        a.axhline(0, color='k', lw=0.7)
        a.set_ylabel('residual [mag]')
        w_ = 1 / ee ** 2
        rms_w = float(np.sqrt(np.average((rr - np.average(rr, weights=w_)) ** 2,
                                         weights=w_)))
        a.set_title(f'{ttl}: rms {np.std(rr):.3f} (weighted {rms_w:.3f})',
                    fontsize=10)
        a.grid(alpha=0.3); a.legend(fontsize=8)
        a.set_ylim(-0.9, 0.9)
    a2.set_xlabel('z')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'figures', 'sdss_hubble_pre_post.png'),
                dpi=170)
    print('wrote sdss_standardized.csv, figures/sdss_cmagic_residual_vs_x1c_post'
          '.png, figures/sdss_hubble_pre_post.png')
    return S


if __name__ == '__main__':
    main()
