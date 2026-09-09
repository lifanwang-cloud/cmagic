"""0.2.0 tests: Mode S covariance algebra, synthesis round-trip, coverage gate."""
import numpy as np
import pytest

sncosmo = pytest.importorskip('sncosmo')
import cmagic
from cmagic import core, highz

BETA = 1.94


def test_mode_s_covariance_algebra():
    # analytic check: two window points, hand-built covariance
    # m_i = (1-b) B_i + b V_i + 0.6 b  -> Var(m) = (1-b)^2 VB + b^2 VV + 2b(1-b) CBV
    b = BETA
    VB, VV, CBV = 0.04 ** 2, 0.03 ** 2, 0.5 * 0.04 * 0.03
    cov = np.zeros((6, 6))
    for i in range(3):
        cov[i, i] = VB
        cov[3 + i, 3 + i] = VV
        cov[i, 3 + i] = cov[3 + i, i] = CBV
    var_m_expected = (1 - b) ** 2 * VB + b ** 2 * VV + 2 * b * (1 - b) * CBV
    # note the sign: Cov(B,V) > 0 REDUCES Var(m) because (1-b) < 0 for b = 1.94
    assert 2 * b * (1 - b) * CBV < 0
    # run Mode S on two synthetic points at colors 0.5 and 0.8
    t0 = 55000.
    # third point beyond the color turnback anchors t(B-V max); the window keeps
    # the first two
    tB = t0 + np.array([8., 12., 24.])
    B = np.array([15.0, 15.5, 16.3]); V = np.array([14.5, 14.7, 15.2])
    fit = core.cmagic_fit(tB, B, np.full(3, np.sqrt(VB)), tB, V,
                          np.full(3, np.sqrt(VV)), 0.0, t0, 1.1, mode='S',
                          cov_BV=cov, sigma_beta=0.16, sigma_int=0.0,
                          k_in_synthesis=True, b_max=14.9, v_max=14.9)
    assert fit['failed'] is None
    # GLS variance must exceed the best single point's slope-free variance floor
    c = B - V
    lever = (c - 0.6)
    Vij = np.full((2, 2), np.nan)
    for i in range(2):
        for j in range(2):
            Vij[i, j] = (var_m_expected if i == j else 0.) \
                + 0.16 ** 2 * lever[i] * lever[j]
    one = np.ones(2)
    var_expected = 1. / float(one @ np.linalg.inv(Vij) @ one)
    assert abs(fit['eB_BV06'] - np.sqrt(var_expected)) < 2e-3
    assert fit['slope_syst_share'] > 0.
    assert fit['n_window'] == 2


def test_synthesis_round_trip_z005():
    # redshift the Hsiao template, observe it in griz, synthesize rest B,V and
    # compare with the template's own rest-frame mags: agreement < 0.01 mag
    z = 0.05
    m = sncosmo.Model(source='hsiao')
    m.set(z=z, t0=55000., amplitude=1e-8)
    rows = []
    rng = np.random.default_rng(1)
    for p in (-5., 0., 6., 10., 15., 20., 25.):
        for b in ('sdssg', 'sdssr', 'sdssi', 'sdssz'):
            t = 55000. + p * (1 + z)
            mag = float(m.bandmag(b, 'ab', t))
            rows.append(dict(mjd=t, band=b[-1], mag=mag, emag=0.01, bandpass=b,
                             source=''))
    syn = highz.synthesize(rows, z, engine='hsiao', n_mc=20, seed=3)
    assert syn['ok']
    rest = sncosmo.Model(source='hsiao')
    rest.set(z=0., t0=0., amplitude=1e-8)
    for p, Bs, Vs in zip(syn['phases'], syn['B'], syn['V']):
        assert abs(Bs - float(rest.bandmag('bessellb', 'vega', p))) < 0.01
        assert abs(Vs - float(rest.bandmag('bessellv', 'vega', p))) < 0.01


def test_rest_band_uncovered_gate():
    # at z = 1.2, rest V redshifts to 12100 A: sdss griz cannot cover it
    rows = [dict(mjd=55000. + i, band='r', mag=23., emag=0.05, bandpass='sdssr',
                 source='') for i in range(6)]
    syn = highz.synthesize(rows, 1.2, engine='hsiao', n_mc=0)
    assert not syn['ok']
    assert 'rest_band_uncovered' in syn['gates']
    r = cmagic.distance(photometry=dict(
        mjd=np.array([55000. + i for i in range(6)], float),
        band=np.array(['r'] * 6), mag=np.full(6, 23.), emag=np.full(6, 0.05)),
        z=1.2, filters={'r': 'sdssr'})
    assert r.failed == 'rest_band_uncovered'


def test_template_insensitivity_crossing():
    # cookbook 8.1: generate data with SALT3-NIR (template B), synthesize with the
    # Hsiao engine (template A) morphed to the observed colors; B_BV0.6 must be
    # recovered within the photometric error scale
    from cmagic import core
    z = 0.12
    gen = sncosmo.Model(source='salt3-nir')
    gen.set(z=z, t0=56000., x1=0.3, c=0.05)
    gen.set_source_peakabsmag(-19.3, 'bessellb', 'vega')
    rows = []
    for p in (-4., 0., 5., 8., 11., 14., 17., 20., 24.):
        for b in ('sdssg', 'sdssr', 'sdssi', 'sdssz'):
            t = 56000. + p * (1 + z)
            try:
                mag = float(gen.bandmag(b, 'ab', t))
            except Exception:
                continue
            rows.append(dict(mjd=t, band=b[-1], mag=mag, emag=0.02, bandpass=b))
    # truth: the generator's own rest-frame locus, fitted the same way
    rest = sncosmo.Model(source='salt3-nir')
    rest.set(z=0., t0=0., x0=gen.parameters[2], x1=0.3, c=0.05)
    ph = np.array([-4., 0., 5., 8., 11., 14., 17., 20., 24.])
    Bt = np.array([float(rest.bandmag('bessellb', 'vega', p)) for p in ph])
    Vt = np.array([float(rest.bandmag('bessellv', 'vega', p)) for p in ph])
    truth = core.cmagic_fit(56000. + ph, Bt, np.full(len(ph), .02),
                            56000. + ph, Vt, np.full(len(ph), .02),
                            0., 56000., 1.05, mode='auto', b_max=Bt.min(),
                            v_max=Vt.min(), k_in_synthesis=True)
    from cmagic import highz
    syn = highz.synthesize(rows, z, engine='hsiao', n_mc=0)
    assert syn['ok']
    g = np.isfinite(syn['B'])
    fit = core.cmagic_fit(56000. + syn['phases'][g], syn['B'][g],
                          np.full(int(g.sum()), .02),
                          56000. + syn['phases'][g], syn['V'][g],
                          np.full(int(g.sum()), .02),
                          0., 56000., 1.05, mode='auto', b_max=syn['b_max'],
                          v_max=syn['v_max'], k_in_synthesis=True)
    assert fit['failed'] is None and truth['failed'] is None
    assert abs(fit['B_BV06_raw'] - truth['B_BV06_raw']) < 0.05
