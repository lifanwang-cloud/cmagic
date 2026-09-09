"""End-to-end regressions on the two published-photometry demo objects.
Targets are this package's own validated outputs (consistent with the archival
CMAGIC chain: SN 2004dt ~68 Mpc; SN 2006X Mode R bracketing the 16.9 Mpc
Cepheid distance between the standard and its anomalous dust law), enforced at
+-0.05 mag in distance modulus."""
import os
import cmagic

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'examples', 'data')


def test_sn2004dt_mode_l():
    r = cmagic.distance(photometry=os.path.join(DATA, 'sn2004dt_loss_bv.csv'),
                        z=0.0194, ebv_mw=0.0218)
    assert r.failed is None
    assert r.mode == 'L'
    assert abs(r.mu - 34.175) <= 0.05
    assert r.beta_free is not None and abs(r.beta_free - 1.955) <= 0.1
    assert r.gates['rms_gate'] == 'pass'


def test_sn2006x_mode_r():
    common = dict(photometry=os.path.join(DATA, 'sn2006x_csp_bv.csv'),
                  z=0.0053, ebv_mw=0.0225)
    r31 = cmagic.distance(**common, rb_true=3.1, e_true='auto')
    r25 = cmagic.distance(**common, rb_true=2.48, e_true='auto')
    assert r31.failed is None and r31.mode == 'R'
    assert abs(r31.mu - 30.586) <= 0.05
    assert abs(r25.mu - 31.370) <= 0.05
    # the two dust laws must bracket the Cepheid distance to M100
    assert r31.D_mpc < 16.9 < r25.D_mpc


def test_named_failure_no_peak():
    # only post-maximum tail data: dm15 unmeasurable -> named gate, no number
    import numpy as np
    rows = dict(mjd=np.arange(53270., 53290.),
                band=np.array(['B', 'V'] * 10),
                mag=np.linspace(16, 17, 20),
                emag=np.full(20, 0.02))
    r = cmagic.distance(photometry=rows, z=0.02)
    assert r.failed == 'dm15_provenance'
    assert r.D_mpc is None
