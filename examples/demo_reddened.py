"""Mode R on a heavily reddened supernova: the branch never reaches B-V = 0.6, so
the pipeline measures at the shifted target c* = 0.6 + E_guess and never
extrapolates; the host dust law R_B is the one open parameter.
Data: SN 2006X (CSP DR3; Krisciunas et al. 2017, AJ 154, 211; VizieR J/AJ/154/211).
SN 2006X has a documented anomalous dust law R_B ~ 2.5 (e.g. Wang X. et al. 2008)."""
import os
import cmagic

here = os.path.dirname(os.path.abspath(__file__))
common = dict(
    photometry=os.path.join(here, 'data', 'sn2006x_csp_bv.csv'),
    z=0.0053, ebv_mw=0.0225,
)
r31 = cmagic.distance(**common, rb_true=3.1, e_true='auto')
r25 = cmagic.distance(**common, rb_true=2.48, e_true='auto',
                      panel=os.path.join(here, 'sn2006x_panel.png'))
print('mode', r31.mode, ' E_guess =', r31.E_guess, ' c* =', r31.c_star,
      ' m* =', r31.m_star)
print('D(R_B = 3.1)  =', r31.D_mpc, 'Mpc')
print('D(R_B = 2.48) =', r25.D_mpc, 'Mpc   (Cepheid distance to M100: 16.9 Mpc)')
