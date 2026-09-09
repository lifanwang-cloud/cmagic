"""v0.3 tests: correlation removal, boundary gate, mode-S inflation, precedence."""
import csv
import os

import numpy as np
import pytest

import cmagic
from cmagic import standardize_sample

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _sdss_rows():
    fn = os.path.join(HERE, 'examples', 'sdss_salt3_comparison.csv')
    rows = []
    for r in csv.DictReader(open(fn)):
        if not (r['cm_mu'] and not r['cm_failed'] and r['x1'] and r['c']
                and not r['salt3_cut']):
            continue
        rows.append(dict(z=float(r['z']), mu=float(r['cm_mu']),
                         emu=float(r['cm_emu'] or 0.15), mode=r['cm_mode'],
                         x1=float(r['x1']), c=float(r['c'])))
    return rows


def test_correlation_removed_on_sdss():
    rows = _sdss_rows()
    assert len(rows) > 20
    S = standardize_sample(rows, covariates=('x1', 'c'), fit_mode_offsets=True)
    # the PI's acceptance criterion: post-fit correlations ~ 0 (weighted r is
    # zero by construction of the WLS fit)
    assert abs(S.post_correlations['x1']) < 0.05
    assert abs(S.post_correlations['c']) < 0.05
    assert 'beta_C' in S.coefficients and 'delta_S' in S.coefficients
    # post rms must not exceed pre rms
    assert S.stats_post['all']['rms'] <= S.stats_pre['all']['rms'] + 1e-6


def test_mode_s_chi2_unity_post_inflation():
    rows = _sdss_rows()
    S = standardize_sample(rows, covariates=('x1', 'c'))
    for k, v in S.stats_post.items():
        if k.endswith('_inflated'):
            assert abs(v['chi2_dof'] - 1.0) < 0.35


def test_precedence_no_double_color_correction():
    rows = [dict(z=0.05 + 0.01 * i, mu=36.7 + 5 * np.log10((0.05 + 0.01 * i) / 0.05),
                 emu=0.1, mode='L', x1=0.5 * (i % 4 - 1.5), c=0.05 * (i % 3 - 1),
                 host_term=0.2) for i in range(12)]
    S_c = standardize_sample(rows, covariates=('x1', 'c'),
                             fit_mode_offsets=False, cosmology='binned')
    S_no = standardize_sample(rows, covariates=('x1',),
                              fit_mode_offsets=False, cosmology='binned')
    # with the fitted color term, the internal host correction is backed out of
    # the base mu (mu + host_term); without it, mu is used as delivered
    d = np.array(S_c.mu_corr) - np.array(S_no.mu_corr)
    # difference must carry the 0.2 host term (up to the different fitted M0/beta)
    assert abs(np.median(d) - 0.2) < 0.1


def test_boundary_gate_stretch():
    sncosmo = pytest.importorskip('sncosmo')
    from cmagic import highz
    z = 0.08
    m = sncosmo.Model(source='hsiao')
    m.set(z=z, t0=55000., amplitude=1e-8)
    rows = []
    for p in (-4., 0., 4., 8., 12., 16., 20.):
        for b in ('sdssg', 'sdssr', 'sdssi'):
            t = 55000. + p * (1 + z) / 1.45     # time-compressed: s_true ~ 0.69
            mag = float(m.bandmag(b, 'ab', 55000. + p * (1 + z)))
            rows.append(dict(mjd=t, band=b[-1], mag=mag, emag=0.02, bandpass=b))
    syn = highz.synthesize(rows, z, engine='hsiao', n_mc=0)
    assert not syn['ok']
    assert 'dm15_provenance' in syn['gates']
