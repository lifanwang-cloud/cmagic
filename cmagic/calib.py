"""Standard-candle calibration: W03 Table 3 linear law M_BV(dm15) on color-restricted
subsamples (H0 = 65 scale), and the WS06 bilinear refinement (kink at dm15 = 1.1; zero
point anchored to W03 Table 3 at dm15 = 1.1 - WS06 do not tabulate one)."""
import numpy as np

# (color_cut, a, b, sigma): M_BV = a + b (dm15 - 1.1); W03 Table 3, R_B = 3.3, H0 = 65,
# 3-sigma-rejection columns of the 'BV' rows.
W03_TABLE3 = [(0.05, -19.35, 0.60, 0.07), (0.20, -19.38, 0.72, 0.11),
              (0.50, -19.44, 0.87, 0.14)]
E0_COEF = (-0.118, 0.249)          # W03 Eq 8: E0 = a + b (dm15 - 1.1)
WS06_BILINEAR_BV = (0.18, 0.25)    # alpha, alpha' about the dm15 = 1.1 kink
BETA_POPULATION = 1.94             # W03 K-corrected mean slope
DM15_CALIB_MAX = 1.7               # beyond: 91bg-like, outside both calibrations


def E0(dm15):
    """W03 Eq 8 intrinsic color locus."""
    return E0_COEF[0] + E0_COEF[1] * (dm15 - 1.1)


def M_BV_w03(dm15, color_corr):
    """W03 Table 3 candle magnitude (H0=65). Returns (M, sigma, bin_label)."""
    cc = color_corr if np.isfinite(color_corr) else 0.5
    row = next((r for r in W03_TABLE3 if cc <= r[0]), W03_TABLE3[-1])
    return row[1] + row[2] * (dm15 - 1.1), row[3], f'<={row[0]:.2f}'


def M_BV_ws06(dm15):
    """WS06 bilinear refinement, anchored to W03 Table 3 at dm15 = 1.1."""
    a, ap = WS06_BILINEAR_BV
    return W03_TABLE3[0][1] + a * (dm15 - 1.1) + ap * abs(dm15 - 1.1)


def rescale_mu(mu65, h0):
    """Explicit H0 rescale of a H0=65-scale distance modulus."""
    return mu65 + 5 * np.log10(65. / h0)
