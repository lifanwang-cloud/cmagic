"""End-to-end regression on one cached SDSS-II object (public photometry;
Sako et al. 2018). Seeded, deterministic; target from the frozen validation run."""
import csv
import os

import numpy as np
import pytest

pytest.importorskip('sncosmo')
import cmagic

CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     'examples', 'data', 'sdss')


def test_sdss_3901_mode_l():
    rows = [r for r in csv.DictReader(
        (ln for ln in open(os.path.join(CACHE, 'sdss_ia_subset_phot.csv'))
         if not ln.startswith('#'))) if r['snid'] == '3901']
    assert len(rows) > 20
    table = dict(mjd=np.array([float(r['mjd']) for r in rows]),
                 band=np.array([r['band'] for r in rows]),
                 mag=np.array([float(r['mag']) for r in rows]),
                 emag=np.array([float(r['emag']) for r in rows]))
    res = cmagic.distance(photometry=table, z=0.0628607,
                          filters={'g': 'sdssg', 'r': 'sdssr', 'i': 'sdssi',
                                   'z': 'sdssz'},
                          engine='auto', ebv_mw=0.0242, n_mc=120, seed=42)
    assert res.failed is None
    assert res.mode == 'L'
    assert abs(res.mu - 37.008) < 0.05
    assert res.K_systematic is not None and res.K_systematic < 0.15
    assert res.gates['rest_band_uncovered'] == 'pass'
