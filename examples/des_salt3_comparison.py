"""SALT3 head-to-head on the DES-SN5YR validation subset (cookbook 10.4).

Runs SALT3 (sncosmo) on the identical cached photometry used by
des_validation.py, with the fixed fiducial Tripp standardization
(alpha = 0.14, beta = 3.1, M_B = -19.36) for continuity with the SDSS-II
section; the symmetric refit lives in make_hubble_residual_fig_symmetric_des.
Cross-check: our (x1, c, mu) against the release's own SALT3 fit
(x1_des, c_des, mu_des from DES-Dovekie_Metadata.csv) on the same objects.

Outputs: examples/des_salt3_comparison.csv and a printed table.
Run:  python examples/des_salt3_comparison.py    (needs cmagic[highz])
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from des_validation import load_cache, mu_lcdm, H0, OM       # noqa: E402

import sncosmo                                               # noqa: E402
from astropy.table import Table                              # noqa: E402

ALPHA, BETA, M_B = 0.14, 3.1, -19.36
X1_CUT, C_CUT = 3.0, 0.3


def salt3_fit(rows, z, mwebv):
    data = dict(time=[], band=[], flux=[], fluxerr=[], zp=[], zpsys=[])
    for r in rows:
        f = 10 ** (-0.4 * (float(r['mag']) - 25.))
        data['time'].append(float(r['mjd']))
        data['band'].append('des' + r['band'])
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
    emu = 0.1
    ex1 = ec = None
    try:
        vn = res.vparam_names
        C = np.array(res.covariance)
        J = np.zeros(len(vn))
        for k, name in enumerate(vn):
            J[k] = {'x0': -2.5 / (np.log(10) * p['x0']), 'x1': ALPHA,
                    'c': -BETA}.get(name, 0.)
        emu = float(np.sqrt(J @ C @ J))
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
    cm = {}
    for name in ('des_validation_assisted.csv', 'des_validation.csv'):
        cmf = os.path.join(HERE, name)
        if os.path.exists(cmf):
            cm = {r['snid']: r for r in csv.DictReader(open(cmf))}
            print(f'CMAGIC results from {name}')
            break
    out = []
    for k, m in enumerate(meta):
        row = dict(snid=m['snid'], z=m['z'])
        try:
            s = salt3_fit(phot[m['snid']], m['z'], m['mwebv'])
            row.update(s)
            row['salt3_cut'] = ('x1' if abs(s['x1']) > X1_CUT else
                                'c' if abs(s['c']) > C_CUT else '')
        except Exception as e:
            row['salt3_cut'] = 'fit_failed'
            row['note'] = str(e)[:50]
        for key in ('x1_des', 'c_des', 'mu_des', 'muerr_des'):
            row[key] = m[key]
        c = cm.get(m['snid'], {})
        row['cm_mu'] = float(c['mu']) if c.get('mu') else None
        row['cm_emu'] = float(c['emu']) if c.get('emu') else None
        row['cm_mode'] = c.get('mode', '')
        row['cm_failed'] = c.get('failed', '')
        row['K_syst'] = float(c['K_syst']) if c.get('K_syst') else None
        row['mu_lcdm'] = round(float(mu_lcdm(m['z'])), 4)
        out.append(row)
        if (k + 1) % 20 == 0:
            print(f'[{k+1}/{len(meta)}]', flush=True)
    s_ok = [r for r in out if r.get('mu') is not None and not r['salt3_cut']]
    both = [r for r in s_ok if r['cm_mu'] is not None and not r['cm_failed']]
    print(f'SALT3 ok: {len(s_ok)}/{len(out)}; overlap with CMAGIC: {len(both)}')
    # cross-check vs the release fits
    dx1 = np.array([r['x1'] - r['x1_des'] for r in s_ok])
    dc = np.array([r['c'] - r['c_des'] for r in s_ok])
    dmu = np.array([r['mu'] - r['mu_des'] for r in s_ok])
    print(f'vs DES release SALT3: d_x1 median {np.median(dx1):+.3f} '
          f'(rms {np.std(dx1):.3f}); d_c median {np.median(dc):+.3f} '
          f'(rms {np.std(dc):.3f}); d_mu median {np.median(dmu):+.3f} '
          f'(rms {np.std(dmu):.3f}) [release MU is bias-corrected]')

    def stats(rs, mukey, emukey):
        res = np.array([r[mukey] - r['mu_lcdm'] for r in rs])
        emu = np.array([r[emukey] or 0.15 for r in rs])
        off = float(np.median(res))
        return off, float(np.std(res - off)), res - off, emu
    off_s, rms_s, rs_s, es = stats(both, 'mu', 'emu')
    off_c, rms_c, rs_c, ec = stats(both, 'cm_mu', 'cm_emu')
    print(f'overlap n={len(both)}: SALT3 rms {rms_s:.3f}, CMAGIC rms {rms_c:.3f}')
    r_p = float(np.corrcoef(rs_s, rs_c)[0, 1]) if len(both) > 3 else np.nan
    er = (1 - r_p ** 2) / np.sqrt(max(len(both) - 3, 1))
    print(f'residual correlation: r = {r_p:.2f} +/- {er:.2f}')
    with open(os.path.join(HERE, 'des_salt3_comparison.csv'), 'w',
              newline='') as f:
        fields = ['snid', 'z', 'mu', 'emu', 'x1', 'ex1', 'c', 'ec', 'salt3_cut',
                  'x1_des', 'c_des', 'mu_des', 'muerr_des', 'cm_mu', 'cm_emu',
                  'cm_mode', 'cm_failed', 'K_syst', 'mu_lcdm']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, '') for k in fields})
    print('wrote des_salt3_comparison.csv')


if __name__ == '__main__':
    main()
