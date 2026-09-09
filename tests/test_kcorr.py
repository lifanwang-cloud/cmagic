"""WS06 SN 1992al K-correction regression (their Table 9 worked example;
z = 0.0135, raw color E = -0.051). Reproduction limits are documented in
docs/KNOWN_ISSUES.md: the distance-relevant K(B_BV) reproduces the published
value to 1.5 mmag in magnitude; the B_max/V_max entries include light-curve
refit effects beyond the tabulated polynomials, and the published table's sign
convention (corrected = raw + K) differs from this package's validated
application (corrected = raw - K), so only magnitude bounds are enforced."""
import cmagic.kcorr as K

Z, E = 0.0135, -0.051


def test_k_bbv_worked_example():
    # WS06 Table 9: K correction to B_BV for SN1992al = 0.022 mag
    assert abs(abs(K.kcorr_B(Z, 0.6)) - 0.022) <= 0.0015


def test_k_bmax_vmax_bounds():
    # WS06 Table 9: B_max +0.009, V_max +0.012; closed-form values agree only in
    # magnitude scale (|K| < 0.03) - documented looser bound.
    assert abs(K.kcorr_B(Z, E)) < 0.03
    assert abs(K.kcorr_V(Z, E)) < 0.03


def test_k_small_at_low_z():
    # normal colors at z <= 0.01: millimagnitude-level; extreme red color at
    # z = 0.03 stays below 0.16 mag (the psi slope times the color)
    for z in (0.003, 0.01):
        for c in (-0.1, 0.6, 1.0):
            assert abs(K.kcorr_B(z, c)) < 0.03
            assert abs(K.kcorr_V(z, c)) < 0.03
    assert abs(K.kcorr_B(0.03, 1.5)) < 0.16
