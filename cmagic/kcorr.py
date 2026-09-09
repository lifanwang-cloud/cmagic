"""WS06 K-correction scheme (Wang, Strovink et al. 2006, ApJ 641, 50, Appendix D).
Corrections to B and V magnitudes are linear in B-V color, with intercept mu(z) and
slope psi(z) given by cubic polynomials in redshift (WS06 Table 10). Convention here:
corrected = raw - K (the convention validated by this package's end-to-end regression
tests; see docs/KNOWN_ISSUES.md for the reproduction limits of the WS06 worked example)."""

_COEF = {'muB': (-0.97119, 28.9034, -187.314), 'psiB': (3.86474, -23.3410, 165.731),
         'muV': (-1.17487, 16.5470, -259.775), 'psiV': (3.04291, -10.4068, 149.714)}


def _cubic(c, z):
    return c[0] * z + c[1] * z ** 2 + c[2] * z ** 3


def mu_B(z):
    return _cubic(_COEF['muB'], z)


def psi_B(z):
    return _cubic(_COEF['psiB'], z)


def kcorr_B(z, color):
    """K correction to a B magnitude at B-V = color."""
    return _cubic(_COEF['muB'], z) + _cubic(_COEF['psiB'], z) * color


def kcorr_V(z, color):
    """K correction to a V magnitude at B-V = color."""
    return _cubic(_COEF['muV'], z) + _cubic(_COEF['psiV'], z) * color
