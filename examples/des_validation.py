"""DES-SN5YR validation of the high-z CMAGIC chain (COOKBOOK section 10.4).

Data: the public DES-SN5YR data release (DES Collaboration 2024; light curves
Sanchez et al. 2024), https://github.com/des-science/DES-SN5YR ,
doi:10.5281/zenodo.12720777. The repository ships a cached subset
(examples/data/des/): 122 cosmology-sample SNe Ia at zHEL 0.05-0.65 (all
z < 0.2 plus 25 per bin above, seed 42), SNANA SMP FLUXCAL converted to AB
magnitudes at SNR >= 3 — the same convention as the SDSS-II cache, so both
validation runs feed the identical pipeline. The release's own SALT3 x1/c/mB
and MU are carried in the meta table for cross-checks.

Precedent: Conley et al. 2006, ApJ 644, 1 — the blind CMAGIC cosmology from
21 high-z SNe ("data sets not observed in a manner optimized for CMAGIC").

Outputs (next to this script): des_hubble.png, des_validation.csv, printed
summary (mode census, rms by z bin, K_syst, failure census, and the
photon-noise share of eB_BV0.6 — the low-z-vs-high-z 'mystery' diagnostic).
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cmagic

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'data', 'des')
H0, OM = 72.0, 0.3
C_KMS = 299792.458
FILTERS = {'g': 'desg', 'r': 'desr', 'i': 'desi', 'z': 'desz'}


def mu_lcdm(z):
    zz = np.linspace(0., z, 300)
    ez = np.sqrt(OM * (1 + zz) ** 3 + (1 - OM))
    dc = C_KMS / H0 * np.trapezoid(1. / ez, zz)
    return 5 * np.log10((1 + z) * dc) + 25.


def load_cache():
    meta, phot = [], {}
    with open(os.path.join(CACHE, 'des_ia_subset_meta.csv')) as fh:
        for r in csv.DictReader(ln for ln in fh if not ln.startswith('#')):
            for k in ('z', 'zhd', 'mwebv', 'pkmjd', 'x1_des', 'ex1_des',
                      'c_des', 'ec_des', 'mB_des', 'mu_des', 'muerr_des'):
                r[k] = float(r[k])
            meta.append(r)
    with open(os.path.join(CACHE, 'des_ia_subset_phot.csv')) as fh:
        for r in csv.DictReader(ln for ln in fh if not ln.startswith('#')):
            phot.setdefault(r['snid'], []).append(r)
    return meta, phot


def dm15_from_release_salt3(x1, c):
    """External dm15: the release SALT3 (x1, c) rest-frame B decline —
    the same convention as the pipeline's internal dm15_rest (highz.py)."""
    import sncosmo
    mr = sncosmo.Model(source='salt3')
    mr.set(z=0., t0=0., x0=1e-10, x1=x1, c=c)
    return float(mr.bandmag('bessellb', 'vega', 15.)
                 - mr.bandmag('bessellb', 'vega', 0.))


def main():
    assisted = '--assisted' in sys.argv
    only = [a for a in sys.argv[1:] if not a.startswith('-')] or None
    meta, phot = load_cache()
    if only:
        meta = [m for m in meta if m['snid'] in only]
    tag = '_assisted' if assisted else ''
    if '--replot' in sys.argv:
        out = []
        for r in csv.DictReader(open(os.path.join(HERE, f'des_validation{tag}.csv'))):
            for k in ('z', 'mu', 'emu', 'mu_lcdm', 'K_syst', 'eB_fit'):
                r[k] = float(r[k]) if r.get(k) else None
            r['failed'] = r['failed'] or None
            out.append(r)
        return finish(out, tag)
    print(f'{len(meta)} DES-SN5YR cosmology-sample Ia in the cache'
          + (' [ASSISTED: t_bmax=release PKMJD, dm15 from release SALT3 x1/c]'
             if assisted else ' [BLIND]'))
    out = []
    for k, m in enumerate(meta):
        rows = phot[m['snid']]
        table = dict(mjd=np.array([float(r['mjd']) for r in rows]),
                     band=np.array([r['band'] for r in rows]),
                     mag=np.array([float(r['mag']) for r in rows]),
                     emag=np.array([float(r['emag']) for r in rows]))
        kw = {}
        if assisted:
            kw = dict(t_bmax=m['pkmjd'],
                      dm15=dm15_from_release_salt3(m['x1_des'], m['c_des']))
        try:
            res = cmagic.distance(photometry=table, z=m['z'], filters=FILTERS,
                                  engine='auto', ebv_mw=m['mwebv'], n_mc=100,
                                  seed=42, h0=H0, **kw)
        except Exception as e:
            out.append(dict(snid=m['snid'], z=m['z'], failed='exception',
                            note=str(e)[:60]))
            print(f'[{k+1}/{len(meta)}] {m["snid"]} EXC {str(e)[:40]}',
                  flush=True)
            continue
        prov = res.provenance or {}
        beta_used = prov.get('beta_used', prov.get('beta_fixed', 1.94))
        rb = prov.get('rb_host', 3.1)
        host_term = (rb - beta_used) * (res.E_host or 0.)
        out.append(dict(
            snid=m['snid'], z=m['z'], failed=res.failed, mode=res.mode,
            mu=res.mu,
            emu=None if res.eD_mpc is None
            else round(5 / np.log(10) * res.eD_mpc / res.D_mpc, 4),
            eB_fit=getattr(res, 'eB_BV06', None),
            E_host=res.E_host, host_term=round(float(host_term), 4),
            D_mpc=res.D_mpc, n_window=res.n_window, K_syst=res.K_systematic,
            rms=res.rms, color_leverage=res.color_leverage,
            slope_share=res.slope_syst_share,
            mu_lcdm=round(float(mu_lcdm(m['z'])), 4)))
        print(f'[{k+1}/{len(meta)}] {m["snid"]} z={m["z"]:.3f} '
              f'{"FAIL:" + str(res.failed) if res.failed else res.mode}',
              flush=True)
    return finish(out, tag)


def finish(out, tag):
    # ---- summary ----
    okr = [r for r in out if not r['failed'] and r.get('mu') is not None]
    fails = [r for r in out if r['failed']]
    resid = np.array([r['mu'] - r['mu_lcdm'] for r in okr])
    zs = np.array([r['z'] for r in okr])
    modes = {}
    for r in okr:
        modes[r['mode']] = modes.get(r['mode'], 0) + 1
    print(f'fit: {len(okr)}; failed: {len(fails)}; modes: {modes}')
    off = float(np.median(resid)) if len(resid) else np.nan
    print(f'median Hubble offset {off:+.3f} mag; rms about median '
          f'{np.std(resid - off):.3f}')
    for zlo, zhi in ((0.05, 0.2), (0.2, 0.35), (0.35, 0.5), (0.5, 0.65)):
        mb = (zs >= zlo) & (zs < zhi)
        if mb.sum():
            print(f'  z {zlo}-{zhi}: n={mb.sum()} '
                  f'rms={np.std(resid[mb] - off):.3f}')
    ks = [r['K_syst'] for r in okr if r.get('K_syst') is not None]
    if ks:
        print(f'K_syst median {np.median(ks):.3f}, 90th {np.percentile(ks, 90):.3f}')
    eb = [(r['z'], r['eB_fit']) for r in okr if r.get('eB_fit')]
    if eb:
        eb = np.array(eb)
        for zlo, zhi in ((0.05, 0.2), (0.2, 0.35), (0.35, 0.5), (0.5, 0.65)):
            mb = (eb[:, 0] >= zlo) & (eb[:, 0] < zhi)
            if mb.sum():
                print(f'  eB_BV06 median z {zlo}-{zhi}: '
                      f'{np.median(eb[mb, 1]):.3f} mag (n={mb.sum()})')
    census = {}
    for r in fails:
        census[r['failed']] = census.get(r['failed'], 0) + 1
    print('failure census:', census)
    # ---- outputs ----
    with open(os.path.join(HERE, f'des_validation{tag}.csv'), 'w',
              newline='') as f:
        fields = ['snid', 'z', 'failed', 'mode', 'mu', 'emu', 'eB_fit', 'E_host',
                  'host_term', 'D_mpc', 'n_window', 'K_syst', 'rms',
                  'color_leverage', 'slope_share', 'mu_lcdm', 'note']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, '') for k in fields})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1]))
    zg = np.linspace(0.04, 0.7, 120)
    a1.plot(zg, [mu_lcdm(x) + off for x in zg], 'k-', lw=1,
            label=f'LCDM (H0={H0}, Om={OM}) + {off:+.2f}')
    for r in okr:
        c = {'L': 'tab:blue', 'Q': 'tab:orange', 'S': 'tab:red',
             'R': 'tab:green'}.get(r['mode'], 'k')
        a1.errorbar(r['z'], r['mu'], yerr=r['emu'] or 0.15, fmt='o', color=c,
                    ms=4, alpha=0.8, elinewidth=0.7)
        a2.errorbar(r['z'], r['mu'] - r['mu_lcdm'] - off, yerr=r['emu'] or 0.15,
                    fmt='o', color=c, ms=4, alpha=0.8, elinewidth=0.7)
    okz = np.argsort([r['z'] for r in okr])
    for k, i in enumerate(okz):
        r = okr[i]
        c = {'L': 'tab:blue', 'Q': 'tab:orange', 'S': 'tab:red',
             'R': 'tab:green'}.get(r['mode'], 'k')
        dx, dy = [(3, 5), (3, -10), (-3, 5), (-3, -10)][k % 4]
        a2.annotate(r['snid'], (r['z'], r['mu'] - r['mu_lcdm'] - off),
                    textcoords='offset points', xytext=(dx, dy),
                    ha='left' if dx > 0 else 'right', fontsize=5.0, color=c)
    for mname, c in (('L', 'tab:blue'), ('Q', 'tab:orange'), ('S', 'tab:red')):
        a1.plot([], [], 'o', color=c, label=f'mode {mname} ({modes.get(mname, 0)})')
    a2.axhline(0, color='k', lw=0.7)
    a1.set_ylabel('mu (CMAGIC)')
    a2.set_ylabel('residual [mag]')
    a2.set_xlabel('heliocentric redshift')
    a1.legend(fontsize=8)
    a1.set_title(f'DES-SN5YR SNe Ia: CMAGIC high-z chain, '
                 f'{"assisted" if tag else "blind"} '
                 f'({len(okr)} fit / {len(out)}; labels = DES SNID)')
    a1.grid(alpha=0.3); a2.grid(alpha=0.3)
    a2.set_ylim(-1, 1)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, f'des_hubble{tag}.png'), dpi=170)
    print(f'wrote des_hubble{tag}.png, des_validation{tag}.csv')


if __name__ == '__main__':
    main()
