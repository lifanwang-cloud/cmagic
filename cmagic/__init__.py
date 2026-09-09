"""cmagic - Type Ia supernova distances from the color-magnitude intercept.
Wang et al. 2003, ApJ 590, 944; Wang, Strovink et al. 2006, ApJ 641, 50.
Public API: cmagic.distance(...) -> CMagicResult. See README.md and docs/METHOD.md."""
from dataclasses import dataclass, field

import numpy as np

from . import io, core, calib, kcorr           # noqa: F401
from .core import cmagic_fit                   # noqa: F401
from .sample import standardize_sample, SampleStandardization  # noqa: F401

__version__ = '0.3.1'
__all__ = ['distance', 'CMagicResult', 'cmagic_fit', 'standardize_sample',
           'SampleStandardization']


@dataclass
class CMagicResult:
    """The result of cmagic.distance(). `failed` is None or the named gate."""
    D_mpc: float | None = None
    eD_mpc: float | None = None
    mu: float | None = None
    B_BV06: float | None = None
    eB_BV06: float | None = None
    mode: str | None = None
    beta_free: float | None = None
    ebeta_free: float | None = None
    E_host: float | None = None
    E_guess: float | None = None
    window: list | None = None
    n_nights: int | None = None
    rms: float | None = None
    gates: dict = field(default_factory=dict)
    flags: list = field(default_factory=list)
    failed: str | None = None
    source: str = ''
    t_bmax: float | None = None
    dm15: float | None = None
    m_star: float | None = None
    c_star: float | None = None
    x1: float | None = None
    c: float | None = None
    chi2dof: float | None = None
    err_scale: float | None = None
    template: str | None = None
    engine: str | None = None
    K_systematic: float | None = None
    n_window: int | None = None
    color_leverage: float | None = None
    slope_syst_share: float | None = None
    provenance: dict = field(default_factory=dict)

    def __repr__(self):
        if self.failed:
            return f'CMagicResult(failed={self.failed!r}, gates={self.gates})'
        return (f'CMagicResult(D_mpc={self.D_mpc}, eD_mpc={self.eD_mpc}, '
                f'mu={self.mu}, mode={self.mode!r})')


def distance(photometry, z, t_bmax=None, dm15=None, ebv_mw=0.0, rb_host=3.1,
             h0=72.0, panel=None, system='vega', source_column=None,
             subtype='normal', rb_mw=4.15, e_true=None, rb_true=None,
             b_max=None, v_max=None, standardization='w03', mode='auto',
             filters=None, engine='auto', bandpass_column=None, n_mc=200,
             seed=0, sigma_int=0.08, sigma_beta=0.16):
    """A CMAGIC distance from a B, V light curve - or, at high redshift / with
    non-BV filters, from multi-band photometry via template synthesis (see
    docs/CMAGIC_COOKBOOK.md Part II). Returns a CMagicResult; on failure
    `result.failed` names the gate and no degraded number is returned."""
    rows_all = io.load_photometry(photometry, source_column=source_column,
                                  keep_all_bands=True,
                                  bandpass_column=bandpass_column)
    if filters:
        for r in rows_all:
            r['bandpass'] = filters.get(r['band'], r.get('bandpass', ''))
    obs_bands = sorted(set(r['band'] for r in rows_all))
    nonbv = any(b not in ('B', 'V') for b in obs_bands)
    if (z > 0.1) or (nonbv and (filters or bandpass_column)):
        return _distance_highz(rows_all, z, dm15=dm15, ebv_mw=ebv_mw,
                               rb_host=rb_host, h0=h0, panel=panel,
                               engine=engine, n_mc=n_mc, seed=seed,
                               sigma_int=sigma_int, sigma_beta=sigma_beta,
                               subtype=subtype, rb_mw=rb_mw,
                               standardization=standardization, mode=mode,
                               e_true=e_true, rb_true=rb_true, system=system)
    rows = io.load_photometry(photometry, source_column=source_column)
    src_used, prov_extra = 'single-source (as supplied)', {}
    if source_column and any(r['source'] for r in rows):
        guess = np.median([r['mjd'] for r in rows]) if rows else 0.
        t_ref = t_bmax if t_bmax is not None else guess

        def chk(rr):
            tb, mb, _ = io.nightly(rr, 'B')
            tv, mv, _ = io.nightly(rr, 'V')
            return len(tb) >= 4 and len(tv) >= 4
        src_used, rows_sel, prov_extra = io.select_source(rows, chk)
        if rows_sel is None:
            res = CMagicResult(failed='too_few_nights',
                               gates={g: 'pass' for g in core.GATE_NAMES})
            res.gates['too_few_nights'] = 'fail'
            res.flags = ['no single source passes the coverage check']
            return res
        rows = rows_sel
    tB, mB, eB = io.nightly(rows, 'B')
    tV, mV, eV = io.nightly(rows, 'V')
    # ---- peak / dm15 provenance ----
    pk = None
    if t_bmax is None or dm15 is None or b_max is None or v_max is None:
        pk = io.measure_peak(tB, mB, tV, mV, z,
                             t0_guess=t_bmax if t_bmax is not None else None)
    if t_bmax is None:
        t_bmax = pk['t_bmax'] if pk else None
    if dm15 is None:
        dm15 = pk['dm15'] if (pk and np.isfinite(pk['dm15'])) else None
    if b_max is None and pk:
        b_max = pk['b_max']
    if v_max is None and pk and np.isfinite(pk['v_max']):
        v_max = pk['v_max']
    if t_bmax is None or dm15 is None or not np.isfinite(dm15) or dm15 <= 0:
        res = CMagicResult(failed='dm15_provenance',
                           gates={g: 'pass' for g in core.GATE_NAMES})
        res.gates['dm15_provenance'] = 'fail'
        res.flags = ['t_bmax/dm15 neither supplied nor measurable from the '
                     'light curve (peak not covered)']
        res.provenance = dict(peak_fit=pk, system=system, **prov_extra)
        return res
    # ---- measurement, with the rerun-style auto fallback L -> Q -> R ----
    kw = dict(subtype=subtype, ebv_mw=ebv_mw, rb_mw=rb_mw, rb_host=rb_host,
              b_max=b_max, v_max=v_max, e_true=e_true, rb_true=rb_true, h0=h0,
              standardization=standardization, source=src_used)
    P = core.cmagic_fit(tB, mB, eB, tV, mV, eV, z, t_bmax, dm15, mode=mode, **kw)
    if (P.get('mode') == 'R' and P.get('ok') and 'modeR_deferred' in
            ';'.join(P.get('flags', [])) and e_true is None):
        P2 = core.cmagic_fit(tB, mB, eB, tV, mV, eV, z, t_bmax, dm15, mode='R',
                             **{**kw, 'e_true': 'auto'})
        if P2.get('ok'):
            P = P2
    elif (not P.get('ok') and P.get('mode') == 'L'
          and P['gates'].get('no_bracket') == 'fail' and mode == 'auto'):
        P2 = core.cmagic_fit(tB, mB, eB, tV, mV, eV, z, t_bmax, dm15, mode='Q',
                             **kw)
        if P2.get('ok'):
            P = P2
    res = CMagicResult(
        D_mpc=P.get('D_mpc'), eD_mpc=P.get('eD_mpc'), mu=P.get('mu'),
        B_BV06=P.get('B_BV_final', P.get('B_BV06_raw')),
        eB_BV06=P.get('eB_BV06'), mode=P.get('mode'),
        beta_free=P.get('beta_free'), ebeta_free=P.get('ebeta_free'),
        E_host=P.get('E_host'), E_guess=P.get('E_guess'),
        window=P.get('window'), n_nights=P.get('n_sel'), rms=P.get('rms'),
        gates=P.get('gates', {}), flags=P.get('flags', []),
        failed=P.get('failed'), source=str(src_used),
        chi2dof=P.get('chi2dof'), err_scale=P.get('err_scale'),
        t_bmax=t_bmax, dm15=float(dm15), m_star=P.get('m_star'),
        c_star=P.get('c_star'),
        provenance={**P, 'peak_fit': pk, 'system': system, **prov_extra})
    if panel:
        from .panel import make_panel
        make_panel(P, panel)
    return res


def _distance_highz(rows_all, z, dm15=None, ebv_mw=0.0, rb_host=3.1, h0=72.0,
                    panel=None, engine='auto', n_mc=200, seed=0, sigma_int=0.08,
                    sigma_beta=0.16, subtype='normal', rb_mw=4.15,
                    standardization='w03', mode='auto', e_true=None,
                    rb_true=None, system='ab'):
    panel_path = panel
    """High-redshift chain: MW de-reddening of the observed bands, template
    synthesis of rest-frame B, V (COOKBOOK section 8), then the Part I fit on the
    synthesized light curve with the synthesis covariance (Mode S when sparse)."""
    from . import highz

    def _fail(gate, note, prov):
        res = CMagicResult(failed=gate, gates={g: 'pass' for g in core.GATE_NAMES})
        res.gates[gate] = 'fail'
        res.flags = [note]
        res.provenance = prov
        if panel_path:
            from .panel import make_panel
            make_panel(dict(failed=gate, flags=[note], gates=res.gates,
                            source=f'synthesis:{engine}', ok=False),
                       panel_path, title=f'CMAGIC high-z FAILED: {gate} (z={z})')
        return res

    rows = [dict(r) for r in rows_all if r.get('bandpass')]
    if not rows:
        return _fail('rest_band_uncovered',
                     'no filter identifiers: pass filters={band: sncosmo_name} '
                     'or bandpass_column=', {})
    # first-order MW de-reddening of the observed bands (CCM, observer frame)
    if ebv_mw:
        import sncosmo as _snc
        for r in rows:
            leff = float(_snc.get_bandpass(r['bandpass']).wave_eff)
            r['mag'] -= float(highz._ccm_shape(np.array([leff]))[0]) * ebv_mw
    eng = 'hsiao' if engine == 'auto' else engine
    syn = highz.synthesize(rows, z, engine=eng, n_mc=n_mc, seed=seed)
    if not syn.get('ok'):
        return _fail(syn['gates'][0] if syn.get('gates') else 'rest_band_uncovered',
                     syn.get('note', 'synthesis failed'), dict(synthesis=str(syn)))
    K_syst = None   # filled by the template-insensitivity test after the fit
    dm15_use = dm15 if dm15 is not None else syn.get('dm15_rest')
    if dm15_use is None or not np.isfinite(dm15_use) or dm15_use <= 0:
        return _fail('dm15_provenance',
                     'rest-frame decline rate neither supplied nor recovered by '
                     'the synthesis', dict(synthesis_params=syn.get('params')))
    t0 = syn['t0']
    tB = t0 + np.asarray(syn['phases'])          # rest-frame clock; z=0 below
    P = core.cmagic_fit(tB, np.asarray(syn['B']), np.asarray(syn['eB']),
                        tB, np.asarray(syn['V']), np.asarray(syn['eV']),
                        0.0, t0, float(dm15_use), subtype=subtype, mode=mode,
                        ebv_mw=0.0, rb_mw=rb_mw, rb_host=rb_host,
                        b_max=syn.get('b_max'), v_max=syn.get('v_max'),
                        e_true=e_true, rb_true=rb_true, h0=h0,
                        standardization=standardization,
                        source=f'synthesis:{eng}', cov_BV=syn.get('cov'),
                        sigma_beta=sigma_beta, sigma_int=sigma_int,
                        k_in_synthesis=True)
    if 'rest_band_uncovered' not in P['gates']:
        P['gates']['rest_band_uncovered'] = 'pass'
    # ---- template-insensitivity test (COOKBOOK section 8.1): offset templates,
    # each morphed to the same observed colors; K_systematic = B_BV0.6 spread ----
    variant_b06 = {}
    if P.get('failed') is None and P.get('B_BV06_raw') is not None:
        try:
            for name, (phv, Bv, Vv) in highz.insensitivity_variants(
                    rows, z, syn, seed=seed).items():
                gv = np.isfinite(Bv) & np.isfinite(Vv)
                if gv.sum() < 1:
                    continue
                Pv = core.cmagic_fit(
                    t0 + phv[gv], Bv[gv], np.full(gv.sum(), 0.03),
                    t0 + phv[gv], Vv[gv], np.full(gv.sum(), 0.03),
                    0.0, t0, float(dm15_use), subtype=subtype,
                    mode=P.get('mode', 'auto'), ebv_mw=0.0, rb_mw=rb_mw,
                    rb_host=rb_host, b_max=syn.get('b_max'),
                    v_max=syn.get('v_max'), h0=h0,
                    standardization=standardization, source=f'variant:{name}',
                    sigma_beta=sigma_beta, sigma_int=sigma_int,
                    k_in_synthesis=True)
                if Pv.get('B_BV06_raw') is not None:
                    variant_b06[name] = float(Pv['B_BV06_raw'])
        except Exception:
            pass
        if variant_b06:
            vals = list(variant_b06.values()) + [float(P['B_BV06_raw'])]
            K_syst = round((max(vals) - min(vals)) / 2., 4)
    res = CMagicResult(
        D_mpc=P.get('D_mpc'), eD_mpc=P.get('eD_mpc'), mu=P.get('mu'),
        B_BV06=P.get('B_BV_final', P.get('B_BV06_raw')),
        eB_BV06=P.get('eB_BV06'), mode=P.get('mode'),
        beta_free=P.get('beta_free'), ebeta_free=P.get('ebeta_free'),
        E_host=P.get('E_host'), E_guess=P.get('E_guess'),
        window=P.get('window'), n_nights=P.get('n_sel'), rms=P.get('rms'),
        gates=P.get('gates', {}), flags=P.get('flags', []),
        failed=P.get('failed'), source=f'synthesis:{eng}',
        chi2dof=P.get('chi2dof'), err_scale=P.get('err_scale'),
        t_bmax=t0, dm15=float(dm15_use), m_star=P.get('m_star'),
        c_star=P.get('c_star'),
        x1=syn.get('x1'), c=syn.get('c'),
        template=eng, engine=eng, K_systematic=K_syst,
        n_window=P.get('n_window', P.get('n_sel')),
        color_leverage=P.get('color_leverage'),
        slope_syst_share=P.get('slope_syst_share'),
        provenance={**P, 'synthesis': {k: v for k, v in syn.items()
                                       if k not in ('cov',)},
                    'template_variants_B_BV06': variant_b06,
                    'system': system, 'true_z': z})
    if panel:
        from .panel import make_panel
        make_panel(P, panel, title=f'CMAGIC high-z (engine {eng}, z={z})')
    return res
