"""Sample-level fitted standardization (cookbook section 6.1, v0.3).

The PI's rule: "the correct way is to always fit the x1 and c correction so that
any correlation is removed." Given a SAMPLE of CMagicResult objects, fit

    resid_i = M0 + delta_S [mode_i == S] + alpha_C x1_i + beta_C c_i

to the Hubble residuals by robust weighted least squares, correct each object's
mu, propagate the coefficient covariance, inflate per-mode errors so the
post-fit chi2/dof is ~1 per mode, and report the post-fit residual correlations
against the covariates (zero by construction of the fit).

PRECEDENCE (no double color correction): when 'c' is among the fitted
covariates, the internal W03 Eqs 7-9 host-extinction correction
(-(R_B - beta) E_host) is REMOVED from the base mu before fitting - the fitted
beta_C c term supersedes it. The internal correction applies only when the user
excludes 'c' from the covariates. Single-object calls are never corrected:
fitted corrections are a sample-level operation, like SALT3 training.
"""
from dataclasses import dataclass, field

import numpy as np

C_KMS = 299792.458


def _mu_lcdm(z, h0=72.0, om=0.3):
    zz = np.linspace(0., z, 300)
    ez = np.sqrt(om * (1 + zz) ** 3 + (1 - om))
    dc = C_KMS / h0 * np.trapezoid(1. / ez, zz)
    return 5 * np.log10((1 + z) * dc) + 25.


@dataclass
class SampleStandardization:
    coefficients: dict = field(default_factory=dict)      # name -> (value, err)
    coefficient_cov: np.ndarray | None = None
    mu_corr: list = field(default_factory=list)           # per object
    emu_corr: list = field(default_factory=list)
    used: list = field(default_factory=list)              # indices fitted
    dropped: list = field(default_factory=list)           # (index, reason)
    inflation: dict = field(default_factory=dict)         # mode -> factor
    stats_pre: dict = field(default_factory=dict)
    stats_post: dict = field(default_factory=dict)
    post_correlations: dict = field(default_factory=dict)  # covariate -> r


def _extract(res, external, i):
    """(z, mu, emu, mode, x1, c, host_term) from a CMagicResult or a dict row."""
    if isinstance(res, dict):
        z = float(res['z']); mu = float(res['mu']); emu = float(res.get('emu', 0.15))
        mode = res.get('mode', 'L')
        x1 = res.get('x1'); c = res.get('c')
        host = float(res.get('host_term', 0.))
    else:
        prov = res.provenance or {}
        z = prov.get('true_z', prov.get('z'))
        mu = res.mu
        emu = (5 / np.log(10)) * res.eD_mpc / res.D_mpc \
            if (res.eD_mpc and res.D_mpc) else 0.15
        mode = res.mode
        x1, c = res.x1, res.c
        beta_used = prov.get('beta_used', prov.get('beta_fixed', 1.94))
        rb = prov.get('rb_host', 3.1)
        host = (rb - beta_used) * (res.E_host or 0.)
    if external is not None:
        ext = external[i] if not isinstance(external, dict) else \
            external.get(getattr(res, 'snid', None) or
                         (res.get('snid') if isinstance(res, dict) else None))
        if ext is not None:
            x1, c = ext
    return z, mu, emu, mode, x1, c, host


def standardize_sample(results, covariates=('x1', 'c'), fit_mode_offsets=True,
                       cosmology=None, h0=72.0, om=0.3, external_covariates=None,
                       clip_sigma=3.0, inflate_by_mode=True):
    """Fit and remove sample-level shape/color/mode terms. `results`: a list of
    CMagicResult (or dict rows with z, mu, emu, mode, x1, c). `cosmology`: None
    (flat LCDM h0/om), a callable mu(z), or 'binned' (subtract per-0.05-z-bin
    medians). Returns SampleStandardization."""
    n = len(results)
    Z, MU, EMU, MODE, X1, C, HOST = [], [], [], [], [], [], []
    S = SampleStandardization()
    fit_c = 'c' in covariates
    for i, r in enumerate(results):
        try:
            z, mu, emu, mode, x1, c, host = _extract(r, external_covariates, i)
        except Exception as e:
            S.dropped.append((i, f'extract: {e}'))
            Z.append(np.nan); MU.append(np.nan); EMU.append(np.nan)
            MODE.append(''); X1.append(np.nan); C.append(np.nan); HOST.append(0.)
            continue
        if mu is None or z is None:
            S.dropped.append((i, 'no mu/z (failed fit)'))
            Z.append(np.nan); MU.append(np.nan); EMU.append(np.nan)
            MODE.append(mode or ''); X1.append(np.nan); C.append(np.nan)
            HOST.append(0.)
            continue
        # precedence: fitted color term supersedes the internal host correction
        mu_base = mu + host if fit_c else mu
        Z.append(z); MU.append(mu_base); EMU.append(max(emu, 0.02))
        MODE.append(mode); HOST.append(host)
        X1.append(np.nan if x1 is None else float(x1))
        C.append(np.nan if c is None else float(c))
    Z, MU, EMU = map(np.asarray, (Z, MU, EMU))
    X1 = np.asarray(X1); C = np.asarray(C)
    MODE = np.asarray(MODE)
    if cosmology == 'binned':
        mu_ref = np.full(n, np.nan)
        for zlo in np.arange(0., 2., 0.05):
            m = (Z >= zlo) & (Z < zlo + 0.05) & np.isfinite(MU)
            if m.sum():
                mu_ref[m] = np.median(MU[m])
    elif callable(cosmology):
        mu_ref = np.array([cosmology(z) if np.isfinite(z) else np.nan for z in Z])
    else:
        mu_ref = np.array([_mu_lcdm(z, h0, om) if np.isfinite(z) else np.nan
                           for z in Z])
    resid = MU - mu_ref
    cols, names = [np.ones(n)], ['M0']
    if fit_mode_offsets:
        cols.append((MODE == 'S').astype(float)); names.append('delta_S')
    if 'x1' in covariates:
        cols.append(X1); names.append('alpha_C')
    if fit_c:
        cols.append(C); names.append('beta_C')
    X = np.column_stack(cols)
    # drop degenerate (zero-variance) covariate columns - a sample where every
    # object shares the same x1 (or c, or mode) cannot constrain that term
    keepcol = [0] + [j for j in range(1, X.shape[1])
                     if np.nanstd(X[:, j]) > 1e-9]
    dropped_terms = [names[j] for j in range(1, X.shape[1]) if j not in keepcol]
    X = X[:, keepcol]
    names = [names[j] for j in keepcol]
    for nm in dropped_terms:
        S.dropped.append((-1, f'degenerate covariate column: {nm}'))
    ok = np.isfinite(resid) & np.all(np.isfinite(X), axis=1) & np.isfinite(EMU)
    for i in np.where(~ok & np.isfinite(MU))[0]:
        S.dropped.append((int(i), 'missing covariate'))
    keep = ok.copy()
    coef = np.zeros(X.shape[1]); cov = np.eye(X.shape[1])
    for _ in range(5):                       # robust: iterative sigma clip
        W = 1 / EMU[keep] ** 2
        Xk = X[keep]
        A = Xk.T @ (Xk * W[:, None])
        cov = np.linalg.inv(A)
        coef = cov @ (Xk.T @ (resid[keep] * W))
        r_all = resid - X @ coef
        sig = float(np.std(r_all[keep]))
        new = ok & (np.abs(r_all) < clip_sigma * max(sig, 0.02))
        if (new == keep).all() or new.sum() < X.shape[1] + 2:
            break
        keep = new
    r_post = resid - X @ coef
    S.coefficients = {nm: (round(float(v), 4), round(float(np.sqrt(cov[j, j])), 4))
                      for j, (nm, v) in enumerate(zip(names, coef))}
    S.coefficient_cov = cov
    S.used = [int(i) for i in np.where(keep)[0]]
    # per-mode chi2 inflation so post chi2/dof ~ 1 per mode
    infl = {}
    for mo in set(MODE[keep]):
        m = keep & (MODE == mo)
        chi2 = float(np.mean((r_post[m] / EMU[m]) ** 2))
        infl[mo] = round(max(np.sqrt(chi2), 1.0), 3) if inflate_by_mode else 1.0
    S.inflation = infl
    # corrected values + propagated errors
    for i in range(n):
        if not np.isfinite(MU[i]) or not np.all(np.isfinite(X[i])):
            S.mu_corr.append(None); S.emu_corr.append(None)
            continue
        prop = float(np.sqrt(X[i] @ cov @ X[i]))
        f = infl.get(MODE[i], 1.0)
        S.mu_corr.append(round(float(MU[i] - X[i] @ coef), 4))
        S.emu_corr.append(round(float(np.hypot(EMU[i] * f, prop)), 4))
    # stats + post-fit correlations
    def _st(r, m):
        return dict(rms=round(float(np.std(r[m])), 4),
                    chi2_dof=round(float(np.mean((r[m] / EMU[m]) ** 2)), 2),
                    n=int(m.sum()))
    S.stats_pre = {'all': _st(resid - np.median(resid[keep]), keep)}
    S.stats_post = {'all': _st(r_post, keep)}
    emu_inf = EMU * np.array([infl.get(m, 1.) for m in MODE])
    for mo in set(MODE[keep]):
        m = keep & (MODE == mo)
        S.stats_pre[f'mode_{mo}'] = _st(resid - np.median(resid[keep]), m)
        S.stats_post[f'mode_{mo}'] = _st(r_post, m)
        S.stats_post[f'mode_{mo}_inflated'] = dict(
            chi2_dof=round(float(np.mean((r_post[m] / emu_inf[m]) ** 2)), 2),
            n=int(m.sum()))
    for nm, arr in (('x1', X1), ('c', C)):
        if nm in covariates:
            m = keep & np.isfinite(arr)
            if m.sum() > 3:
                # weighted correlation: zero by construction of the WLS fit
                w = 1 / EMU[m] ** 2
                xw = arr[m] - np.average(arr[m], weights=w)
                rw = r_post[m] - np.average(r_post[m], weights=w)
                num = np.average(xw * rw, weights=w)
                den = np.sqrt(np.average(xw ** 2, weights=w)
                              * np.average(rw ** 2, weights=w))
                S.post_correlations[nm] = round(float(num / den), 4)
                S.post_correlations[nm + '_unweighted'] = round(float(
                    np.corrcoef(arr[m], r_post[m])[0, 1]), 4)
            else:
                S.post_correlations[nm] = None
    return S
