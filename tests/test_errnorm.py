"""v0.3.1: per-fit chi2/dof error normalization (the PI's rule, named case 5103)."""
import csv
import os

import numpy as np
import pytest

import cmagic
from cmagic import core, standardize_sample

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_rescaling_overdispersed_fit():
    # synthetic mode-L fit whose points scatter 3x their quoted errors:
    # errors must scale by sqrt(chi2/dof)
    rng = np.random.default_rng(7)
    t0 = 55000.
    ph = np.array([5., 8., 11., 14., 17., 20., 24.])
    c = 0.25 + 0.03 * ph
    B = 15. + 1.94 * (c - 0.6) + rng.standard_normal(len(ph)) * 0.06
    tB = t0 + ph
    fit = core.cmagic_fit(tB, B, np.full(len(ph), 0.02), tB, B - c,
                          np.full(len(ph), 0.02), 0., t0, 1.1, mode='L',
                          b_max=14.6, v_max=14.6, k_in_synthesis=True)
    assert fit['failed'] is None
    assert fit['chi2dof'] > 2.
    assert abs(fit['err_scale'] - np.sqrt(fit['chi2dof'])) < 5e-3
    # the scaled error must match the naive error times the scale
    naive = np.sqrt(1. / np.sum(1 / 0.02 ** 2 * np.ones(len(ph))))
    assert abs(fit['eB_BV06'] - fit['err_scale'] * naive) < 5e-3


def test_underdispersed_not_shrunk():
    # chi2/dof < 1 must NOT shrink errors (s = max(1, sqrt(chi2/dof)))
    t0 = 55000.
    ph = np.array([4., 7., 10., 13., 16., 19., 22., 25., 28.])
    c = 0.20 + 0.03 * ph
    B = 15. + 1.94 * (c - 0.6)          # perfect line, huge quoted errors
    tB = t0 + ph
    fit = core.cmagic_fit(tB, B, np.full(len(ph), 0.1), tB, B - c,
                          np.full(len(ph), 0.1), 0., t0, 1.1, mode='L',
                          b_max=14.6, v_max=14.6, k_in_synthesis=True)
    assert fit['err_scale'] == 1.0


@pytest.mark.filterwarnings('ignore')
def test_5103_regression():
    sncosmo = pytest.importorskip('sncosmo')
    rows = [r for r in csv.DictReader(
        (ln for ln in open(os.path.join(HERE, 'examples', 'data', 'sdss',
                                        'sdss_ia_subset_phot.csv'))
         if not ln.startswith('#'))) if r['snid'] == '5103']
    table = dict(mjd=np.array([float(r['mjd']) for r in rows]),
                 band=np.array([r['band'] for r in rows]),
                 mag=np.array([float(r['mag']) for r in rows]),
                 emag=np.array([float(r['emag']) for r in rows]))
    res = cmagic.distance(photometry=table, z=0.161989,
                          filters={'g': 'sdssg', 'r': 'sdssr', 'i': 'sdssi',
                                   'z': 'sdssz'},
                          engine='auto', ebv_mw=0.0426, n_mc=60, seed=42)
    assert res.failed is None and res.mode == 'L'
    assert res.chi2dof is not None and res.chi2dof > 1.5      # inconsistent points
    assert res.err_scale > 1.2                                 # error enlarged


def test_sample_inflation_shrinks_after_rescaling():
    # per-fit rescaled emu (err_scale applied) must reduce the sample-level
    # mode-L inflation factor relative to unscaled emu
    fn = os.path.join(HERE, 'examples', 'sdss_salt3_comparison.csv')
    rows, rows_scaled = [], []
    rng = np.random.default_rng(3)
    for r in csv.DictReader(open(fn)):
        if not (r['cm_mu'] and not r['cm_failed'] and r['x1'] and r['c']
                and not r['salt3_cut']):
            continue
        base = dict(z=float(r['z']), mu=float(r['cm_mu']),
                    emu=float(r['cm_emu'] or 0.15), mode=r['cm_mode'],
                    x1=float(r['x1']), c=float(r['c']))
        rows.append(dict(base))
        b2 = dict(base); b2['emu'] = base['emu'] * 1.4   # emulated per-fit rescale
        rows_scaled.append(b2)
    S0 = standardize_sample(rows, covariates=('x1', 'c'))
    S1 = standardize_sample(rows_scaled, covariates=('x1', 'c'))
    assert S1.inflation['L'] < S0.inflation['L']
