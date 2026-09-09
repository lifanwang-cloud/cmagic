"""The CMAGIC measurement engine: window selection, Modes L/Q/R, K corrections,
extinction, standardization. See docs/METHOD.md. All functions operate on nightly
B, V arrays; magnitudes Vega-consistent."""
import numpy as np

from . import kcorr, calib

WINDOWS = {'normal': (7., 30.), '91T': (10., 29.)}
GATE_NAMES = ('window_empty', 'too_few_nights', 'no_bracket', 'rms_gate',
              'beta_gate', 'dm15_provenance', 'photometry_defect')


def _gates_init():
    return {g: 'pass' for g in GATE_NAMES}


def cmagic_fit(tB, mB, eB, tV, mV, eV, z, t_bmax, dm15,
               subtype='normal', mode='auto', ebv_mw=0.0, rb_mw=4.15, rb_host=3.1,
               beta_fixed=calib.BETA_POPULATION, b_max=None, v_max=None,
               e_true=None, rb_true=None, h0=72., standardization='w03',
               source='unspecified', cov_BV=None, sigma_beta=0.16, sigma_int=0.08,
               k_in_synthesis=False):
    """One CMAGIC measurement. Returns a provenance dict; 'failed' is None or the
    first fatal gate name; 'gates' is every gate evaluated (pass/fail)."""
    P = dict(mode_requested=mode, subtype=subtype, z=z, dm15=float(dm15),
             ebv_mw=ebv_mw, rb_mw=rb_mw, rb_host=rb_host, beta_fixed=beta_fixed,
             h0=h0, standardization=standardization, source=source,
             gates=_gates_init(), failed=None, flags=[], ok=False)

    def fail(gate, note=''):
        P['gates'][gate] = 'fail'
        if P['failed'] is None:
            P['failed'] = gate
        if note:
            P['flags'].append(note)
        return P

    win = WINDOWS['91T' if subtype == '91T' else 'normal']
    if len(tB) < 3 or len(tV) < 3:
        return fail('too_few_nights', '<3 nightly points in B or V')
    ph = (tB - t_bmax) / (1 + z)
    Vi = np.interp(tB, tV, mV)
    near = np.array([np.min(np.abs(tV - t)) < 2.5 for t in tB])
    cand = (ph > 0) & near & np.isfinite(mB) & np.isfinite(Vi)
    o = np.argsort(ph[cand])
    php = ph[cand][o]; colp = (mB - Vi)[cand][o]
    Bp = mB[cand][o]; eBp = eB[cand][o]
    P['locus'] = dict(phases=[round(float(x), 2) for x in php],
                      colors=[round(float(x), 3) for x in colp],
                      B=[round(float(x), 3) for x in Bp])
    if len(php) < 1:
        return fail('too_few_nights', 'no paired post-maximum nights')
    if len(php) < 3 and mode in ('L', 'Q', 'R'):
        return fail('too_few_nights', '<3 paired post-maximum nights')
    sm = (np.convolve(colp, np.ones(3) / 3, 'same') if len(colp) >= 3
          else colp.copy())
    if len(colp) >= 3:
        sm[0] = np.mean(colp[:2]); sm[-1] = np.mean(colp[-2:])
    i_turn = int(np.argmax(sm >= sm.max() - 0.02))
    t_bvmax = float(php[i_turn])
    P['t_BVmax'] = round(t_bvmax, 1)
    lo = win[0] / dm15
    hi = min(win[1] / dm15, t_bvmax - 3.)
    P['window'] = [round(lo, 1), round(hi, 1)]
    if hi <= lo:
        return fail('window_empty',
                    f'window [{lo:.1f},{hi:.1f}] closed by the color maximum '
                    f'at +{t_bvmax:.1f} d')
    sel = (php >= lo) & (php <= hi)
    n = int(sel.sum())
    P.update(n_sel=n,
             ph_span=[round(float(php[sel].min()), 1),
                      round(float(php[sel].max()), 1)] if n else None,
             col_span=[round(float(colp[sel].min()), 3),
                       round(float(colp[sel].max()), 3)] if n else None,
             sel=dict(phases=[round(float(x), 2) for x in php[sel]],
                      colors=[round(float(x), 3) for x in colp[sel]],
                      B=[round(float(x), 3) for x in Bp[sel]]))
    sel_idx = np.where(cand)[0][o][sel] if cov_BV is not None else None
    if n < 1:
        return fail('too_few_nights', 'no nights in the window')
    if n < 4 and mode in ('L', 'Q', 'R'):
        return fail('too_few_nights', f'{n} nights in the window (mode {mode})')
    c = colp[sel]; B = Bp[sel]; w = 1 / eBp[sel] ** 2
    span = float(c.max() - c.min())
    gap = float(np.max(np.diff(np.sort(c)))) if n > 1 else 99.
    P.update(col_span_mag=round(span, 3), max_color_gap=round(gap, 3))
    brackets06 = (c.min() - 0.05 <= 0.6 <= c.max() + 0.05)
    reddened = c.min() > 0.6
    m = mode
    if m == 'auto':
        if n >= 4:
            m = 'L' if brackets06 else ('R' if reddened else 'Q')
        else:
            m = 'R' if reddened else 'S'
    P['mode'] = m

    def wlsq(deg, x0):
        A = np.vander(c - x0, deg + 1, increasing=True)
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        coef = cov @ (A.T @ (B * w))
        res = B - A @ coef
        return coef, cov, float(np.sqrt(np.mean(res ** 2)))

    def free_diag(x0):
        coef, cov, rms_free = wlsq(1, x0)
        bf, ebf = float(coef[1]), float(np.sqrt(cov[1, 1]))
        valid = bool(n >= 5 and span >= 0.25 and gap <= 0.4 and rms_free <= 0.15
                     and 1.5 <= bf <= 2.5)
        P.update(beta_free=round(bf, 3), ebeta_free=round(ebf, 3),
                 rms_free=round(rms_free, 4), beta_free_valid=valid)
        if not valid:
            P['gates']['beta_gate'] = 'fail'   # diagnostic-level, not fatal
        return bf

    beta_used = beta_fixed
    if m == 'L':
        if not brackets06:
            return fail('no_bracket', 'branch does not span B-V = 0.6')
        b06_f = float(np.average(B - beta_fixed * (c - 0.6), weights=w))
        rms_f = float(np.sqrt(np.mean((B - (b06_f + beta_fixed * (c - 0.6))) ** 2)))
        P.update(B_BV06_raw=round(b06_f, 4), rms=round(rms_f, 4),
                 eB_BV06=round(max(float(np.sqrt(1 / w.sum())),
                                   rms_f / np.sqrt(n)), 4))
        free_diag(0.6)
        if rms_f > 0.15:
            return fail('rms_gate', f'forced-slope rms {rms_f:.3f} > 0.15')
        raw06 = b06_f
    elif m == 'Q':
        if not brackets06:
            return fail('no_bracket', 'branch does not span B-V = 0.6')
        coef, cov, rms_q = wlsq(2, 0.6)
        P.update(B_BV06_raw=round(float(coef[0]), 4),
                 eB_BV06=round(float(np.sqrt(cov[0, 0])), 4), rms=round(rms_q, 4),
                 quad_coefs=[round(float(x), 4) for x in coef])
        if rms_q > 0.15:
            return fail('rms_gate', f'quadratic rms {rms_q:.3f} > 0.15')
        raw06 = float(coef[0])
        beta_used = float(coef[1])
    elif m == 'S':
        # ---- Mode S (sparse; COOKBOOK section 9): fixed-slope per-point estimates
        # combined by GLS with the full covariance ----
        mi = B - beta_fixed * (c - 0.6)
        nS = len(mi)
        if cov_BV is not None and sel_idx is not None:
            nn = cov_BV.shape[0] // 2
            J = np.zeros((nS, 2 * nn))
            for a_, i_ in enumerate(sel_idx):
                J[a_, i_] = 1. - beta_fixed
                J[a_, nn + i_] = beta_fixed
            Cm = J @ cov_BV @ J.T
        else:
            phs_sel = php[sel]
            eV_near = np.array([eV[np.argmin(np.abs((tV - t_bmax) / (1 + z) - p))]
                                for p in phs_sel]) if len(eV) else np.zeros(nS)
            Cm = np.diag((1 - beta_fixed) ** 2 * eBp[sel] ** 2
                         + beta_fixed ** 2 * eV_near ** 2)
        Vslope = sigma_beta ** 2 * np.outer(c - 0.6, c - 0.6)
        Vij = Cm + Vslope + np.eye(nS) * sigma_int ** 2
        one = np.ones(nS)
        Vinv = np.linalg.inv(Vij)
        var = 1. / float(one @ Vinv @ one)
        mhat = float(var * (one @ Vinv @ mi))
        Vij_noslope = Cm + np.eye(nS) * sigma_int ** 2
        Vinv0 = np.linalg.inv(Vij_noslope)
        var0 = 1. / float(one @ Vinv0 @ one)
        P.update(B_BV06_raw=round(mhat, 4), eB_BV06=round(float(np.sqrt(var)), 4),
                 rms=round(float(np.std(mi)), 4) if nS > 1 else 0.0,
                 n_window=nS, color_leverage=round(float(abs(np.mean(c) - 0.6)), 3),
                 slope_syst_share=round(float(max(var - var0, 0.) / var), 3),
                 sigma_beta=sigma_beta, sigma_int=sigma_int)
        P['flags'].append(f'modeS: {nS} point(s), color leverage '
                          f'{P["color_leverage"]}, slope-systematic share '
                          f'{P["slope_syst_share"]}')
        raw06 = mhat
    else:   # Mode R
        if b_max is None or v_max is None or not np.isfinite(b_max + v_max):
            return fail('no_bracket', 'Mode R requires B_max and V_max for E_guess')
        e_guess = float(b_max - v_max)
        cstar = 0.6 + e_guess
        P.update(E_guess=round(e_guess, 3), c_star=round(cstar, 3))
        if not (c.min() - 0.05 <= cstar <= c.max() + 0.05):
            return fail('no_bracket',
                        f'branch does not span the shifted target c*={cstar:.2f}')
        bstar = float(np.average(B - beta_fixed * (c - cstar), weights=w))
        rms_f = float(np.sqrt(np.mean((B - (bstar + beta_fixed * (c - cstar))) ** 2)))
        P.update(m_star=round(bstar, 4), rms=round(rms_f, 4), beta_used=beta_fixed,
                 eB_BV06=round(max(float(np.sqrt(1 / w.sum())),
                                   rms_f / np.sqrt(n)), 4))
        free_diag(cstar)
        if rms_f > 0.15:
            return fail('rms_gate', f'forced-slope rms {rms_f:.3f} > 0.15')
        raw06 = None
    # ---- K corrections (corrected = raw - K); skipped when the synthesis already
    # performed the cross-filter K-correction (high-z path) ----
    E_raw = float(b_max - v_max) if (b_max is not None and v_max is not None
                                     and np.isfinite(b_max + v_max)) else 0.0
    if k_in_synthesis:
        K_bbv = K_bmax = K_dm15 = 0.0
        P['flags'].append('K-corrections inside the template synthesis')
    else:
        K_bbv = kcorr.kcorr_B(z, 0.6)
        K_bmax = kcorr.kcorr_B(z, E_raw)
        c15 = float(np.interp(15., php, colp))
        K_dm15 = kcorr.psi_B(z) * (c15 - E_raw)
    P.update(K_B_BV=round(K_bbv, 4), K_Bmax=round(K_bmax, 4),
             K_dm15=round(K_dm15, 4))
    dm15_k = dm15 - K_dm15
    if m in ('L', 'Q', 'S'):
        b06_k = raw06 - K_bbv
        b06_mw = b06_k - (rb_mw - beta_used) * ebv_mw
        bint = b06_mw - 0.6 * beta_used
        calE0 = calib.E0(dm15_k)
        if b_max is not None and v_max is not None and np.isfinite(b_max + v_max):
            bmax_c = (b_max - K_bmax) - rb_mw * ebv_mw
            vmax_c = v_max - (rb_mw - 1.) * ebv_mw
            calE = (bmax_c - bint) / beta_used
            E_host = max(calE - calE0, 0.)
            color_corr = bmax_c - vmax_c - E_host
        else:
            E_host, color_corr = 0., np.nan
            P['flags'].append('no_peak_coverage: host extinction assumed 0')
        b06_final = bint - (rb_host - beta_used) * E_host
        P.update(E_host=round(float(E_host), 4), B_BV_intercept=round(bint, 4),
                 B_BV_final=round(b06_final, 4),
                 color_corr=round(float(color_corr), 4)
                 if np.isfinite(color_corr) else None)
    else:
        if e_true is None:
            P['flags'].append('modeR_deferred: E_true not supplied; m*, c*, E_guess '
                              'delivered, no distance')
            P['ok'] = True
            return P
        if e_true == 'auto':
            bint_star = (P['m_star'] - K_bbv) - beta_used * cstar
            bmax_c = (b_max - K_bmax) - rb_mw * ebv_mw
            calE = (bmax_c - bint_star) / beta_used
            e_true = max(calE - calib.E0(dm15_k) - ebv_mw, 0.)
            P['E_true_auto'] = round(float(e_true), 4)
        rb_use = rb_true if rb_true is not None else rb_host
        b06 = (P['m_star'] - K_bbv) - beta_used * (cstar - 0.6 - e_true) \
            - rb_use * e_true
        b06_final = b06 - 0.6 * beta_used - (rb_mw - beta_used) * ebv_mw
        E_host = float(e_true)
        color_corr = float(P['E_guess'] - e_true - ebv_mw)
        P.update(E_host=round(E_host, 4), rb_used=rb_use,
                 B_BV_final=round(b06_final, 4), color_corr=round(color_corr, 4))
    if dm15_k >= calib.DM15_CALIB_MAX:
        P['flags'].append('uncalibrated: dm15 >= 1.7 outside the M_BV calibration')
    if standardization == 'ws06':
        M = calib.M_BV_ws06(dm15_k); sig = 0.10; binlab = 'ws06_bilinear'
    else:
        M, sig, binlab = calib.M_BV_w03(
            dm15_k, color_corr if (color_corr is not None
                                   and np.isfinite(color_corr)) else 0.5)
    mu65 = b06_final - M
    mu = calib.rescale_mu(mu65, h0)
    D = float(10 ** ((mu - 25) / 5))
    P.update(M_BV=round(float(M), 4), calib_bin=binlab, sigma_calib=sig,
             mu65=round(float(mu65), 4), mu=round(float(mu), 4),
             D_mpc=round(D, 3),
             eD_mpc=round(D * np.log(10) / 5
                          * float(np.sqrt(sig ** 2 + P.get('eB_BV06', 0.05) ** 2
                                          + ((rb_host - beta_used) * 0.08) ** 2
                                          * (1 if P.get('E_host') is not None
                                             else 0))), 3),
             ok=True)
    return P
