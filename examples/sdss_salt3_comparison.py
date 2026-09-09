"""SALT3 head-to-head comparison on the SDSS-II validation subset (PI request).

Runs SALT3 (sncosmo) on the identical cached photometry used by
sdss_validation.py, applies the standard Tripp standardization
(alpha = 0.14, beta = 3.1, M_B = -19.36; mB = -2.5 log10(x0) + 10.635), and
compares with the CMAGIC results in sdss_validation.csv against
LCDM (H0 = 72, Om = 0.3). SALT3 validity cuts: fit converged, |x1| < 3,
|c| < 0.3 (standard cosmology-sample cuts).

Outputs: examples/figures/sdss_salt3_vs_cmagic.png,
examples/sdss_salt3_comparison.csv, and a printed comparison table.
Run:  python examples/sdss_salt3_comparison.py    (needs cmagic[highz])
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from sdss_validation import load_cache, mu_lcdm, H0, OM      # noqa: E402

import sncosmo                                               # noqa: E402
from astropy.table import Table                              # noqa: E402

ALPHA, BETA, M_B = 0.14, 3.1, -19.36
X1_CUT, C_CUT = 3.0, 0.3


def salt3_fit(rows, z, mwebv):
    data = dict(time=[], band=[], flux=[], fluxerr=[], zp=[], zpsys=[])
    for r in rows:
        f = 10 ** (-0.4 * (float(r['mag']) - 25.))
        data['time'].append(float(r['mjd']))
        data['band'].append('sdss' + r['band'])
        data['flux'].append(f)
        data['fluxerr'].append(f * float(r['emag']) * np.log(10) / 2.5)
        data['zp'].append(25.); data['zpsys'].append('ab')
    tab = Table(data)
    model = sncosmo.Model(source='salt3', effects=[sncosmo.F99Dust(r_v=3.1)],
                          effect_names=['mw'], effect_frames=['obs'])
    model.set(z=z, mwebv=mwebv)
    t0g = float(tab['time'][np.argmax(tab['flux'])])
    res, fitted = sncosmo.fit_lc(tab, model, ['t0', 'x0', 'x1', 'c'],
                                 bounds={'t0': (t0g - 20, t0g + 20),
                                         'x1': (-4, 4), 'c': (-0.4, 0.6)})
    p = dict(zip(res.param_names, res.parameters))
    mB = -2.5 * np.log10(p['x0']) + 10.635
    mu = mB + ALPHA * p['x1'] - BETA * p['c'] - M_B
    # mu error from the fit covariance (x0, x1, c terms)
    emu = 0.1
    try:
        vn = res.vparam_names
        C = np.array(res.covariance)
        J = np.zeros(len(vn))
        for k, name in enumerate(vn):
            J[k] = {'x0': -2.5 / (np.log(10) * p['x0']), 'x1': ALPHA,
                    'c': -BETA}.get(name, 0.)
        emu = float(np.sqrt(J @ C @ J))
    except Exception:
        pass
    ex1 = ec = None
    try:
        vn = res.vparam_names
        C = np.array(res.covariance)
        if 'x1' in vn:
            ex1 = float(np.sqrt(C[vn.index('x1'), vn.index('x1')]))
        if 'c' in vn:
            ec = float(np.sqrt(C[vn.index('c'), vn.index('c')]))
    except Exception:
        pass
    return dict(t0=p['t0'], x1=p['x1'], c=p['c'], ex1=ex1, ec=ec, mu=float(mu),
                emu=max(emu, 0.02), chisq=float(res.chisq), ndof=int(res.ndof))


def main():
    meta, phot = load_cache()
    cm = {r['snid']: r for r in csv.DictReader(open(
        os.path.join(HERE, 'sdss_validation.csv')))}
    out = []
    for m in meta:
        row = dict(snid=m['snid'], z=m['z'])
        try:
            s = salt3_fit(phot[m['snid']], m['z'], m['mwebv'])
            row.update(s)
            if abs(s['x1']) > X1_CUT:
                row['salt3_cut'] = 'x1'
            elif abs(s['c']) > C_CUT:
                row['salt3_cut'] = 'c'
            else:
                row['salt3_cut'] = ''
        except Exception as e:
            row['salt3_cut'] = 'fit_failed'
            row['note'] = str(e)[:50]
        c = cm.get(m['snid'], {})
        row['cm_mu'] = float(c['mu']) if c.get('mu') else None
        row['cm_emu'] = float(c['emu']) if c.get('emu') else None
        row['cm_mode'] = c.get('mode', '')
        row['cm_failed'] = c.get('failed', '')
        row['K_syst'] = float(c['K_syst']) if c.get('K_syst') else None
        row['mu_lcdm'] = float(mu_lcdm(m['z']))
        out.append(row)
    s_ok = [r for r in out if r.get('mu') is not None and not r['salt3_cut']]
    c_ok = [r for r in out if r['cm_mu'] is not None and not r['cm_failed']]
    both = [r for r in out if r in s_ok and r['cm_mu'] is not None
            and not r['cm_failed']]
    print(f'SALT3: {len(s_ok)}/{len(out)} pass fit + cuts '
          f'(|x1|<{X1_CUT}, |c|<{C_CUT}); CMAGIC: {len(c_ok)}; '
          f'overlap: {len(both)}')

    def stats(rs, mukey, emukey):
        res = np.array([r[mukey] - r['mu_lcdm'] for r in rs])
        emu = np.array([r[emukey] or 0.15 for r in rs])
        off = float(np.median(res))
        rms = float(np.std(res - off))
        chi2 = float(np.mean(((res - off) / emu) ** 2))
        return off, rms, chi2, res - off
    off_s, rms_s, chi_s, rs_s = stats(both, 'mu', 'emu')
    off_c, rms_c, chi_c, rs_c = stats(both, 'cm_mu', 'cm_emu')
    bothL = [r for r in both if r['cm_mode'] == 'L']
    _, rms_sL, _, _ = stats(bothL, 'mu', 'emu')
    _, rms_cL, _, _ = stats(bothL, 'cm_mu', 'cm_emu')
    print('method | N  | grey offset | rms  | chi2/dof | rms(mode L only)')
    print(f'SALT3  | {len(both):2d} | {off_s:+.3f}      | {rms_s:.3f} | '
          f'{chi_s:5.2f}    | {rms_sL:.3f}')
    print(f'CMAGIC | {len(both):2d} | {off_c:+.3f}      | {rms_c:.3f} | '
          f'{chi_c:5.2f}    | {rms_cL:.3f}')
    # correlation of residuals (independence test; W03/WS06: ~0.15 at low z)
    r_p = float(np.corrcoef(rs_s, rs_c)[0, 1])
    er = (1 - r_p ** 2) / np.sqrt(max(len(both) - 3, 1))
    print(f'residual correlation (overlap): Pearson r = {r_p:.2f} +/- {er:.2f}')
    # Delta-mu outliers
    print('Delta-mu = mu_cmagic - mu_salt3 (offset-corrected); >3 sigma flagged:')
    nflag = 0
    for r, ds, dc in zip(both, rs_s, rs_c):
        dmu = dc - ds
        sig = np.sqrt((r['cm_emu'] or 0.15) ** 2 + (r['emu'] or 0.1) ** 2)
        if abs(dmu) > 3 * sig:
            nflag += 1
            print(f"  {r['snid']:6s} z={r['z']:.3f} mode={r['cm_mode']} "
                  f"dmu={dmu:+.2f} ({abs(dmu) / sig:.1f} sigma) "
                  f"K_syst={r['K_syst']}")
    print(f'{nflag} objects flagged of {len(both)}')
    # ---- CSV ----
    with open(os.path.join(HERE, 'sdss_salt3_comparison.csv'), 'w',
              newline='') as f:
        fields = ['snid', 'z', 'mu', 'emu', 'x1', 'ex1', 'c', 'ec',
                  'salt3_cut', 'cm_mu',
                  'cm_emu', 'cm_mode', 'cm_failed', 'K_syst', 'mu_lcdm']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, '') for k in fields})
    # ---- figure ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    a1, a2, a3 = axes
    zg = np.linspace(0.04, 0.37, 80)
    a1.plot(zg, [mu_lcdm(x) for x in zg], 'k-', lw=1, label='LCDM(72, 0.3)')
    for r in s_ok:
        a1.errorbar(r['z'], r['mu'] - off_s, yerr=r['emu'], fmt='s', color='grey',
                    ms=4, alpha=0.7, elinewidth=0.6)
    for r in c_ok:
        col = 'tab:blue' if r['cm_mode'] == 'L' else 'tab:red'
        a1.errorbar(r['z'], r['cm_mu'] - off_c, yerr=r['cm_emu'] or 0.15, fmt='o',
                    color=col, ms=4, alpha=0.8, elinewidth=0.6)
    a1.plot([], [], 's', color='grey', label=f'SALT3 ({len(s_ok)})')
    a1.plot([], [], 'o', color='tab:blue', label='CMAGIC L')
    a1.plot([], [], 'o', color='tab:red', label='CMAGIC S')
    a1.set_xlabel('z'); a1.set_ylabel('mu (own grey offset removed)')
    a1.legend(fontsize=7); a1.grid(alpha=0.3)
    a1.set_title('SDSS-II: SALT3 vs CMAGIC')
    for r, ds, dc in zip(both, rs_s, rs_c):
        col = 'tab:blue' if r['cm_mode'] == 'L' else 'tab:red'
        a2.errorbar(ds, dc, xerr=r['emu'], yerr=r['cm_emu'] or 0.15, fmt='o',
                    color=col, ms=5, alpha=0.8, elinewidth=0.6)
    lim = max(np.max(np.abs(rs_s)), np.max(np.abs(rs_c))) * 1.1
    a2.plot([-lim, lim], [-lim, lim], 'k:', lw=0.8)
    a2.axhline(0, color='k', lw=0.5); a2.axvline(0, color='k', lw=0.5)
    a2.set_xlabel('SALT3 residual [mag]')
    a2.set_ylabel('CMAGIC residual [mag]')
    a2.set_title(f'overlap n={len(both)}: Pearson r = {r_p:.2f} +/- {er:.2f}')
    a2.grid(alpha=0.3)
    for r, ds, dc in zip(both, rs_s, rs_c):
        col = 'tab:blue' if r['cm_mode'] == 'L' else 'tab:red'
        a3.errorbar(r['z'], dc - ds,
                    yerr=np.hypot(r['emu'], r['cm_emu'] or 0.15), fmt='o',
                    color=col, ms=5, alpha=0.8, elinewidth=0.6)
    a3.axhline(0, color='k', lw=0.6)
    a3.set_xlabel('z'); a3.set_ylabel('mu_CMAGIC - mu_SALT3 [mag]')
    a3.set_title('per-object method difference')
    a3.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'figures', 'sdss_salt3_vs_cmagic.png'),
                dpi=170)
    print('wrote figures/sdss_salt3_vs_cmagic.png, sdss_salt3_comparison.csv')


if __name__ == '__main__':
    main()
