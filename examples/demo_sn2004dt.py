"""A clean linear-mode CMAGIC distance from published photometry.
Data: SN 2004dt (KAIT/LOSS; Ganeshalingam et al. 2010, ApJS 190, 418; VizieR
J/ApJS/190/418). The peak is covered, so t_bmax and dm15 are measured from the
light curve automatically."""
import os
import cmagic

here = os.path.dirname(os.path.abspath(__file__))
result = cmagic.distance(
    photometry=os.path.join(here, 'data', 'sn2004dt_loss_bv.csv'),
    z=0.0194,                # heliocentric redshift of SN 2004dt
    ebv_mw=0.0218,           # Schlafly & Finkbeiner Milky-Way E(B-V)
    panel=os.path.join(here, 'sn2004dt_panel.png'),
)
print(result)
print(f'  mode {result.mode}, window {result.window} d, n={result.n_nights}, '
      f'rms={result.rms}')
print(f'  beta_free = {result.beta_free} +- {result.ebeta_free}')
print(f'  D = {result.D_mpc} +- {result.eD_mpc} Mpc  (mu = {result.mu})')
