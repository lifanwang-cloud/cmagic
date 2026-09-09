"""cmagic - Type Ia supernova distances from the color-magnitude intercept.
Wang et al. 2003, ApJ 590, 944; Wang, Strovink et al. 2006, ApJ 641, 50.
Public API: cmagic.distance(...) -> CMagicResult. See README.md and docs/METHOD.md."""
from dataclasses import dataclass, field

import numpy as np

from . import io, core, calib, kcorr           # noqa: F401
from .core import cmagic_fit                   # noqa: F401

__version__ = '0.1.0'
__all__ = ['distance', 'CMagicResult', 'cmagic_fit']


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
    provenance: dict = field(default_factory=dict)

    def __repr__(self):
        if self.failed:
            return f'CMagicResult(failed={self.failed!r}, gates={self.gates})'
        return (f'CMagicResult(D_mpc={self.D_mpc}, eD_mpc={self.eD_mpc}, '
                f'mu={self.mu}, mode={self.mode!r})')


def distance(photometry, z, t_bmax=None, dm15=None, ebv_mw=0.0, rb_host=3.1,
             h0=72.0, panel=None, system='vega', source_column=None,
             subtype='normal', rb_mw=4.15, e_true=None, rb_true=None,
             b_max=None, v_max=None, standardization='w03', mode='auto'):
    """A CMAGIC distance from a B, V light curve. See README.md for the interface
    contract and docs/METHOD.md for the method. Returns a CMagicResult; on failure
    `result.failed` names the gate and no degraded number is returned."""
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
        t_bmax=t_bmax, dm15=float(dm15), m_star=P.get('m_star'),
        c_star=P.get('c_star'),
        provenance={**P, 'peak_fit': pk, 'system': system, **prov_extra})
    if panel:
        from .panel import make_panel
        make_panel(P, panel)
    return res
