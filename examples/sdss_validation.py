"""SDSS-II SN survey validation of the high-z CMAGIC chain (COOKBOOK section 10).

Data: the public SDSS-II Supernova Survey data release (Sako et al. 2018,
PASP 130, 064002), obtained from the release mirror
  https://portal.nersc.gov/project/dessn/SDSS/dataRelease/
    sdsssn_master.dat2                (master table: types, redshifts)
    SDSS_dataRelease-snana.tar.gz     (SNANA-format ugriz light curves)

The repository ships a cached subset (examples/data/sdss/) sufficient to re-run
this script: spectroscopically confirmed SNe Ia at z = 0.05-0.35 with >= 4
epochs in the rest-frame CMAGIC window range. SDSS photometry is public; the
subset CSV carries its provenance in the header. To rebuild the cache from the
full release, set SDSS_RAW=/path/to/downloaded/release and run with
--rebuild-cache.

Outputs (written next to this script): sdss_hubble.png, sdss_validation.csv,
and a printed summary (N fit, mode census, Hubble rms by z bin, template
systematic, failure census by gate).
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cmagic

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'data', 'sdss')
N_TARGET = 60
ZMIN, ZMAX = 0.05, 0.35
H0, OM = 72.0, 0.3
C_KMS = 299792.458


def mu_lcdm(z):
    zz = np.linspace(0., z, 300)
    ez = np.sqrt(OM * (1 + zz) ** 3 + (1 - OM))
    dc = C_KMS / H0 * np.trapezoid(1. / ez, zz)
    return 5 * np.log10((1 + z) * dc) + 25.


# ---------------------------------------------------------------- cache building
def rebuild_cache(raw_dir):
    from astropy.io import fits
    os.makedirs(CACHE, exist_ok=True)
    master = {}
    with open(os.path.join(raw_dir, 'sdsssn_master.dat2')) as fh:
        hdr = fh.readline().split()
        icid, icls, iz = hdr.index('CID'), hdr.index('Classification'), \
            hdr.index('zspecHelio')
        for ln in fh:
            p = ln.split()
            if len(p) <= max(icid, icls, iz):
                continue
            master[p[icid]] = (p[icls], p[iz])
    fh_head = fits.open(os.path.join(
        raw_dir, 'SDSS_dataRelease-snana', 'SDSS_allCandidates+BOSS',
        'SDSS_allCandidates+BOSS_HEAD.FITS'))[1].data
    fh_phot = fits.open(os.path.join(
        raw_dir, 'SDSS_dataRelease-snana', 'SDSS_allCandidates+BOSS',
        'SDSS_allCandidates+BOSS_PHOT.FITS'))[1].data
    rows_out, meta_out = [], []
    for i in range(len(fh_head)):
        snid = str(fh_head['SNID'][i]).strip()
        cls, zs = master.get(snid, ('', ''))
        if cls != 'SNIa':                      # spectroscopically confirmed only
            continue
        try:
            z = float(zs)
        except ValueError:
            continue
        if not (ZMIN <= z <= ZMAX):
            continue
        pmin, pmax = int(fh_head['PTROBS_MIN'][i]) - 1, int(fh_head['PTROBS_MAX'][i])
        ph = fh_phot[pmin:pmax]
        pk = float(fh_head['PEAKMJD'][i])
        ok = (ph['FLUXCAL'] > 0) & (ph['FLUXCALERR'] > 0) & \
             (ph['FLUXCAL'] / ph['FLUXCALERR'] > 3.)
        band = np.char.strip(ph['FLT'].astype(str))
        keep = ok & np.isin(band, ('g', 'r', 'i', 'z'))
        if keep.sum() < 8:
            continue
        rest = (ph['MJD'][keep] - pk) / (1 + z)
        nwin = len({round(float(x)) for x in rest if 3. < x < 30.})
        if nwin < 4:
            continue
        meta_out.append(dict(snid=snid, z=z, mwebv=float(fh_head['MWEBV'][i]),
                             peakmjd=pk, n_phot=int(keep.sum())))
        for j in np.where(keep)[0]:
            f = float(ph['FLUXCAL'][j]); ef = float(ph['FLUXCALERR'][j])
            rows_out.append(dict(snid=snid, mjd=float(ph['MJD'][j]),
                                 band=str(band[j]),
                                 mag=round(27.5 - 2.5 * np.log10(f), 4),
                                 emag=round(1.0857 * ef / f, 4)))
        if len(meta_out) >= N_TARGET:
            break
    hdr = ('# SDSS-II Supernova Survey photometry subset (public data release:\n'
           '# Sako et al. 2018, PASP 130, 064002;\n'
           '# https://portal.nersc.gov/project/dessn/SDSS/dataRelease/).\n'
           '# Spectroscopically confirmed SNe Ia, z 0.05-0.35; SNANA FLUXCAL\n'
           '# converted to AB magnitudes (27.5 - 2.5 log10 FLUXCAL), SNR > 3.\n')
    with open(os.path.join(CACHE, 'sdss_ia_subset_phot.csv'), 'w') as f:
        f.write(hdr)
        f.write('snid,mjd,band,mag,emag\n')
        for r in rows_out:
            f.write(f"{r['snid']},{r['mjd']:.3f},{r['band']},{r['mag']},{r['emag']}\n")
    with open(os.path.join(CACHE, 'sdss_ia_subset_meta.csv'), 'w') as f:
        f.write(hdr)
        f.write('snid,z,mwebv,peakmjd,n_phot\n')
        for r in meta_out:
            f.write(f"{r['snid']},{r['z']},{r['mwebv']},{r['peakmjd']},{r['n_phot']}\n")
    print(f'cache rebuilt: {len(meta_out)} SNe, {len(rows_out)} photometry rows')


def load_cache():
    meta, phot = [], {}
    with open(os.path.join(CACHE, 'sdss_ia_subset_meta.csv')) as fh:
        for r in csv.DictReader(ln for ln in fh if not ln.startswith('#')):
            r['z'] = float(r['z']); r['mwebv'] = float(r['mwebv'])
            r['peakmjd'] = float(r['peakmjd'])
            meta.append(r)
    with open(os.path.join(CACHE, 'sdss_ia_subset_phot.csv')) as fh:
        for r in csv.DictReader(ln for ln in fh if not ln.startswith('#')):
            phot.setdefault(r['snid'], []).append(r)
    return meta, phot


# ---------------------------------------------------------------- the blind run
def main():
    if '--rebuild-cache' in sys.argv:
        rebuild_cache(os.environ['SDSS_RAW'])
    meta, phot = load_cache()
    print(f'{len(meta)} spectroscopically confirmed SDSS-II Ia in the cache')
    FILTERS = {'g': 'sdssg', 'r': 'sdssr', 'i': 'sdssi', 'z': 'sdssz'}
    out = []
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
            out.append(dict(snid=m['snid'], z=m['z'], failed='exception',
                            note=str(e)[:60]))
            continue
        out.append(dict(snid=m['snid'], z=m['z'], failed=res.failed,
                        mode=res.mode, mu=res.mu, emu=None if res.eD_mpc is None
                        else round(5 / np.log(10) * res.eD_mpc / res.D_mpc, 4),
                        D_mpc=res.D_mpc, n_window=res.n_window,
                        K_syst=res.K_systematic, rms=res.rms,
                        color_leverage=res.color_leverage,
                        slope_share=res.slope_syst_share,
                        mu_lcdm=round(float(mu_lcdm(m['z'])), 4)))
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
    print(f'median Hubble offset {off:+.3f} mag (grey zero point; see cookbook), '
          f'rms about the median {np.std(resid - off):.3f} mag')
    for zlo, zhi in ((0.05, 0.15), (0.15, 0.25), (0.25, 0.35)):
        mband = (zs >= zlo) & (zs < zhi)
        if mband.sum():
            print(f'  z {zlo}-{zhi}: n={mband.sum()} rms='
                  f'{np.std(resid[mband] - off):.3f} mag')
    ks = [r['K_syst'] for r in okr if r.get('K_syst') is not None]
    if ks:
        print(f'template-insensitivity spread (cookbook 8.1): median {np.median(ks):.3f}'
              f' mag, 90th pct {np.percentile(ks, 90):.3f}')
    census = {}
    for r in fails:
        census[r['failed']] = census.get(r['failed'], 0) + 1
    print('failure census:', census)
    # ---- outputs ----
    with open(os.path.join(HERE, 'sdss_validation.csv'), 'w', newline='') as f:
        fields = ['snid', 'z', 'failed', 'mode', 'mu', 'emu', 'D_mpc', 'n_window',
                  'K_syst', 'rms', 'color_leverage', 'slope_share', 'mu_lcdm',
                  'note']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, '') for k in fields})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1]))
    zg = np.linspace(0.03, 0.4, 100)
    a1.plot(zg, [mu_lcdm(x) + off for x in zg], 'k-', lw=1,
            label=f'LCDM (H0={H0}, Om={OM}) + {off:+.2f}')
    for r in okr:
        c = {'L': 'tab:blue', 'Q': 'tab:orange', 'S': 'tab:red',
             'R': 'tab:green'}.get(r['mode'], 'k')
        a1.errorbar(r['z'], r['mu'], yerr=r['emu'] or 0.15, fmt='o', color=c,
                    ms=4, alpha=0.8)
        a2.errorbar(r['z'], r['mu'] - r['mu_lcdm'] - off, yerr=r['emu'] or 0.15,
                    fmt='o', color=c, ms=4, alpha=0.8)
    for mname, c in (('L', 'tab:blue'), ('Q', 'tab:orange'), ('S', 'tab:red')):
        a1.plot([], [], 'o', color=c, label=f'mode {mname} ({modes.get(mname, 0)})')
    a2.axhline(0, color='k', lw=0.7)
    a1.set_ylabel('mu (CMAGIC)'); a2.set_ylabel('residual [mag]')
    a2.set_xlabel('heliocentric redshift')
    a1.legend(fontsize=8)
    a1.set_title(f'SDSS-II SNe Ia: CMAGIC high-z chain, blind '
                 f'({len(okr)} fit / {len(out)})')
    a1.grid(alpha=0.3); a2.grid(alpha=0.3)
    a2.set_ylim(-1, 1)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'sdss_hubble.png'), dpi=170)
    print('wrote sdss_hubble.png, sdss_validation.csv')


if __name__ == '__main__':
    main()
