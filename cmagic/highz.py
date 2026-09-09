"""High-redshift extension (docs/CMAGIC_COOKBOOK.md Part II, sections 8-8.1):
template synthesis of rest-frame B, V light curves from observer-frame multi-band
photometry.

Engines:
  'hsiao'     - Hsiao spectral template; global (t0, stretch) from the multi-band
                shape, per-epoch grey + CCM-shaped color warp fitted to all observed
                bands of that epoch (the warp is what makes the K-correction a
                color correction; cookbook section 8.1).
  'salt3-nir' - SALT3-NIR (Pierel et al. 2022); (x1, c) determined ITERATIVELY
                (fit template -> synthesize rest LCs -> refit -> loop; gate
                `kcorr_no_convergence`), per cookbook section 8.1.

The template-insensitivity test (section 8.1) is provided by `insensitivity_variants`:
recompute the synthesis with deliberately offset templates, each morphed to the same
observed colors; the caller reports the B_BV0.6 spread as K_systematic.

Requires the optional dependency sncosmo (pip install cmagic[highz]).
"""
import numpy as np

try:
    import sncosmo as _snc
except ImportError:                                    # pragma: no cover
    _snc = None

_REST_B, _REST_V = 4400., 5500.
_MSG = ("cmagic.highz requires the optional dependency 'sncosmo' "
        "(pip install cmagic[highz] or pip install sncosmo)")


def _require():
    if _snc is None:
        raise ImportError(_MSG)


def register_filter(name, wave_A, trans):
    """Register a user filter (JWST NIRCam etc.) for `filters=` mappings.
    wave_A in Angstrom, trans dimensionless throughput."""
    _require()
    bp = _snc.Bandpass(np.asarray(wave_A, float), np.asarray(trans, float),
                       name=name)
    _snc.register(bp, name, force=True)
    return name


def _ccm_shape(lam_A, rv=3.1):
    """CCM89 A_lambda/E(B-V), optical/NIR, unit E(B-V)."""
    x = 1e4 / np.asarray(lam_A, float)
    a = np.ones_like(x); b = np.zeros_like(x)
    m = (x >= 1.1) & (x <= 3.3)
    y = x[m] - 1.82
    a[m] = (1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4
            + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7)
    b[m] = (1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4
            - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7)
    m2 = x < 1.1
    a[m2] = 0.574 * np.maximum(x[m2], 0.3) ** 1.61
    b[m2] = -0.527 * np.maximum(x[m2], 0.3) ** 1.61
    return rv * a + b


def _epochs(rows, gap=0.6):
    rr = sorted(rows, key=lambda r: r['mjd'])
    out = []
    for r in rr:
        if out and r['mjd'] - out[-1][-1]['mjd'] < gap:
            out[-1].append(r)
        else:
            out.append([r])
    return out


def _magsys_of(b):
    return 'vega' if b.startswith(('bessell', 'csp', '2mass')) else 'ab'


def _coverage_gate(bandnames, z):
    lo = min(_snc.get_bandpass(b).minwave() for b in bandnames)
    hi = max(_snc.get_bandpass(b).maxwave() for b in bandnames)
    missing = [rb for rb in (_REST_B, _REST_V)
               if not (lo <= rb * (1 + z) <= hi)]
    return missing, (lo, hi)


def _sn_table(rows):
    from astropy.table import Table
    data = dict(time=[], band=[], flux=[], fluxerr=[], zp=[], zpsys=[])
    for r in rows:
        f = 10 ** (-0.4 * (r['mag'] - 25.))
        data['time'].append(r['mjd']); data['band'].append(r['bandpass'])
        data['flux'].append(f)
        data['fluxerr'].append(f * r['emag'] * np.log(10) / 2.5)
        data['zp'].append(25.); data['zpsys'].append(_magsys_of(r['bandpass']))
    return Table(data)


# ---------------------------------------------------------------- warp synthesis
def _band_dust_coeffs(model_obs, model_rest, bands, t0):
    """Band-integrated CCM tilt coefficients: the synthetic response of each band
    (through its actual transmission, on the current template SED) to E(B-V) = 0.1,
    in mag per unit E(B-V). More faithful than a lambda_eff point evaluation."""
    ccm = {}
    E = 0.1
    dusty = _snc.Model(source=model_obs.source,
                       effects=[_snc.CCM89Dust()], effect_names=['w'],
                       effect_frames=['obs'])
    pars = dict(zip(model_obs.param_names, model_obs.parameters))
    dusty.set(**{k: v for k, v in pars.items() if k in dusty.param_names})
    dusty.set(webv=E)
    for b in bands:
        try:
            ccm[b] = float(dusty.bandmag(b, _magsys_of(b), t0)
                           - model_obs.bandmag(b, _magsys_of(b), t0)) / E
        except Exception:
            ccm[b] = float(_ccm_shape(np.array(
                [float(_snc.get_bandpass(b).wave_eff)]))[0])
    dr = _snc.Model(source=model_rest.source,
                    effects=[_snc.CCM89Dust()], effect_names=['w'],
                    effect_frames=['obs'])
    pr = dict(zip(model_rest.param_names, model_rest.parameters))
    dr.set(**{k: v for k, v in pr.items() if k in dr.param_names})
    dr.set(webv=E)
    out = {}
    for name, b in (('B', 'bessellb'), ('V', 'bessellv')):
        try:
            out[name] = float(dr.bandmag(b, 'vega', 0.)
                              - model_rest.bandmag(b, 'vega', 0.)) / E
        except Exception:
            out[name] = float(_ccm_shape(np.array(
                [_REST_B if name == 'B' else _REST_V]))[0])
    return ccm, out['B'], out['V']


def warped_rest_BV(rows, z, model_obs, model_rest, t0, s=1.0, mag_lookup=None):
    """Per-epoch grey + CCM-tilt warp of an arbitrary base template to the observed
    photometry, then rest-frame B, V synthesis (cookbook sections 8 and 8.1: 'any
    reasonable template, morphed to the observed colors'). The tilt coefficients are
    band-integrated on the current template. Returns (phases_rest, B, V) - NaN where
    an epoch cannot be warped."""
    if mag_lookup is None:
        mag_lookup = lambda r: r['mag']      # noqa: E731
    eps = _epochs(rows)
    bands = sorted(set(r['bandpass'] for r in rows))
    ccm, ccmB, ccmV = _band_dust_coeffs(model_obs, model_rest, bands, t0)
    tmid = np.array([np.mean([r['mjd'] for r in e]) for e in eps])
    ph = (tmid - t0) / (1 + z)
    B, V, warps = [], [], []
    for e, p_ in zip(eps, ph):
        A, y = [], []
        for r in e:
            try:
                msyn = model_obs.bandmag(r['bandpass'], _magsys_of(r['bandpass']),
                                         t0 + (r['mjd'] - t0) / s)
            except Exception:
                continue
            A.append([1., ccm[r['bandpass']]])
            y.append(mag_lookup(r) - msyn)
        if not y:
            B.append(np.nan); V.append(np.nan); warps.append((np.nan, np.nan))
            continue
        A = np.array(A); y = np.array(y)
        if len(y) == 1:
            g, k = float(y[0]), 0.
        else:
            sol, *_ = np.linalg.lstsq(A, y, rcond=None)
            g, k = float(sol[0]), float(sol[1])
        warps.append((g, k))
        try:
            mB = model_rest.bandmag('bessellb', 'vega', p_ / s)
            mV = model_rest.bandmag('bessellv', 'vega', p_ / s)
        except Exception:
            B.append(np.nan); V.append(np.nan)
            continue
        B.append(mB + g + k * ccmB)
        V.append(mV + g + k * ccmV)
    return ph, np.array(B), np.array(V), warps


def _hsiao_pair(z, t0):
    mo = _snc.Model(source='hsiao'); mo.set(z=z, t0=t0, amplitude=1.)
    mr = _snc.Model(source='hsiao'); mr.set(z=0., t0=0., amplitude=1.)
    return mo, mr


def _salt_pair(z, t0, x1, c, x0=1e-10):
    mo = _snc.Model(source='salt3-nir')
    mo.set(z=z, t0=t0, x0=x0, x1=x1, c=c)
    mr = _snc.Model(source='salt3-nir')
    mr.set(z=0., t0=0., x0=x0, x1=x1, c=c)
    return mo, mr


# ---------------------------------------------------------------- main synthesis
def synthesize(rows, z, engine='hsiao', n_mc=200, seed=0, t0_guess=None):
    """Synthesize rest-frame nightly B, V (dicts rows: mjd, mag, emag, bandpass).
    Returns dict(ok, gates, phases[rest d], B, V, eB, eV, cov, t0, params, b_max,
    v_max, dm15_rest, iterations[for salt3-nir])."""
    _require()
    rng = np.random.default_rng(seed)
    P = dict(engine=engine, z=z, gates=[], ok=False, n_mc=n_mc)
    rows = [r for r in rows if r.get('bandpass')]
    bands = sorted(set(r['bandpass'] for r in rows))
    if not bands:
        P['gates'].append('rest_band_uncovered')
        P['note'] = 'no bandpass identifiers on the photometry'
        return P
    missing, span = _coverage_gate(bands, z)
    P['filter_span_A'] = [round(span[0]), round(span[1])]
    if missing:
        P['gates'].append('rest_band_uncovered')
        P['note'] = (f'rest {missing} x (1+z) outside the observed filter span '
                     f'{P["filter_span_A"]} A - synthesis would extrapolate')
        return P
    tab = _sn_table(rows)
    if t0_guess is None:
        t0_guess = float(tab['time'][np.argmax(tab['flux'])])

    if engine == 'salt3-nir':
        # ---- iterative (x1, c): fit -> synthesize -> refit (cookbook 8.1) ----
        x1, c = 0.0, 0.0
        trace = []
        t0 = t0_guess
        converged = False
        for it in range(8):
            model = _snc.Model(source='salt3-nir')
            model.set(z=z, x1=x1, c=c)
            res, fitted = _snc.fit_lc(tab, model, ['t0', 'x0'],
                                      bounds={'t0': (t0_guess - 20,
                                                     t0_guess + 20)})
            pars = dict(zip(res.param_names, res.parameters))
            t0, x0 = pars['t0'], pars['x0']
            mo, mr = _salt_pair(z, t0, x1, c, x0)
            ph, B, V, _ = warped_rest_BV(rows, z, mo, mr, t0)
            gsel = np.isfinite(B) & np.isfinite(V)
            if gsel.sum() < 3:
                break
            # refit the template to the synthesized rest-frame light curves
            from astropy.table import Table
            fB = 10 ** (-0.4 * (B[gsel] - 25.))
            fV = 10 ** (-0.4 * (V[gsel] - 25.))
            tsyn = Table(dict(
                time=np.concatenate([t0 + ph[gsel], t0 + ph[gsel]]),
                band=['bessellb'] * int(gsel.sum()) + ['bessellv'] * int(gsel.sum()),
                flux=np.concatenate([fB, fV]),
                fluxerr=np.concatenate([0.03 * fB, 0.03 * fV]),
                zp=[25.] * 2 * int(gsel.sum()),
                zpsys=['vega'] * 2 * int(gsel.sum())))
            m2 = _snc.Model(source='salt3-nir'); m2.set(z=0.)
            try:
                r2, f2 = _snc.fit_lc(tsyn, m2, ['t0', 'x0', 'x1', 'c'],
                                     bounds={'t0': (t0 - 10, t0 + 10),
                                             'x1': (-4, 4), 'c': (-0.4, 1.0)})
                p2 = dict(zip(r2.param_names, r2.parameters))
            except Exception:
                break
            dx1, dc = p2['x1'] - x1, p2['c'] - c
            trace.append(dict(it=it, x1=round(float(p2['x1']), 3),
                              c=round(float(p2['c']), 4),
                              dx1=round(float(dx1), 3), dc=round(float(dc), 4)))
            x1, c = float(p2['x1']), float(p2['c'])
            if abs(dx1) < 0.05 and abs(dc) < 0.01:
                converged = True
                break
        P['iterations'] = trace
        if not converged:
            P['gates'].append('kcorr_no_convergence')
            P['note'] = f'(x1, c) not converged after {len(trace)} iterations'
            return P
        if abs(x1) > 3.9 or c <= -0.39 or c >= 0.99:
            P['gates'].append('dm15_provenance')
            P['note'] = (f'template parameter at bound (x1={x1:.2f}, c={c:.3f}): '
                         'no trustworthy shape/window')
            P['iterations'] = trace
            return P
        P['params'] = dict(t0=float(t0), x0=float(x0), x1=x1, c=c)
        P['x1'] = round(x1, 3); P['c'] = round(c, 4)
        mo, mr = _salt_pair(z, t0, x1, c, x0)
        ph, B0, V0, _ = warped_rest_BV(rows, z, mo, mr, t0)
        s = 1.0
        draws = []
        for _ in range(n_mc if n_mc > 1 else 0):
            pert = {id(r): r['mag'] + rng.standard_normal() * r['emag']
                    for r in rows}
            _, Bd, Vd, _ = warped_rest_BV(rows, z, mo, mr, t0,
                                          mag_lookup=lambda r: pert[id(r)])
            draws.append(np.concatenate([Bd, Vd]))
        b_max = float(mr.bandmag('bessellb', 'vega', 0.))
        v_max = float(mr.bandmag('bessellv', 'vega', 0.))
        dm15 = float(mr.bandmag('bessellb', 'vega', 15.) - b_max)
        # tie the absolute level to the near-peak warp (mr has arbitrary x0)
        eps0 = _epochs(rows)
        i_near = int(np.argmin(np.abs(ph)))
        off = B0[i_near] - float(mr.bandmag('bessellb', 'vega', ph[i_near]))
        b_max += off; v_max += off
    elif engine == 'hsiao':
        mo, mr = _hsiao_pair(z, 0.)
        ccm0, _, _ = _band_dust_coeffs(mo, mr, bands, t0_guess)
        eps_all = _epochs(rows)

        def total_resid(t0_, s_):
            # rms of the per-epoch residuals AFTER the grey + tilt warp: the
            # warp-invariant part of the template mismatch, which is what (t0, s)
            # must minimize (epochs with < 3 bands carry no residual)
            mo.set(t0=t0_)
            tot, npt = 0., 0
            for e in eps_all:
                A, y = [], []
                for r in e:
                    try:
                        msyn = mo.bandmag(r['bandpass'],
                                          _magsys_of(r['bandpass']),
                                          t0_ + (r['mjd'] - t0_) / s_)
                    except Exception:
                        continue
                    A.append([1., ccm0[r['bandpass']]])
                    y.append(r['mag'] - msyn)
                if len(y) < 3:
                    continue
                A = np.array(A); y = np.array(y)
                sol, *_ = np.linalg.lstsq(A, y, rcond=None)
                res = y - A @ sol
                tot += float(res @ res); npt += len(y)
            return np.sqrt(tot / max(npt, 1))
        t0s = np.arange(t0_guess - 8, t0_guess + 8.01, 2.0)
        ss = np.arange(0.75, 1.3, 0.05)
        grid = [(total_resid(tt, s_), tt, s_) for tt in t0s for s_ in ss]
        _, t0, s = min(grid)
        for step in (1.0, 0.5):
            cand = [(total_resid(tt, s), tt) for tt in
                    np.arange(t0 - 2 * step, t0 + 2 * step + 0.01, step)]
            _, t0 = min(cand)
        if s <= 0.76 or s >= 1.27:
            P['gates'].append('dm15_provenance')
            P['note'] = (f'stretch at the grid bound (s={s:.2f}): no trustworthy '
                         'shape/window')
            return P
        P['params'] = dict(t0=float(t0), stretch=float(s))
        mo.set(t0=t0)
        ph, B0, V0, warps0 = warped_rest_BV(rows, z, mo, mr, t0, s)
        # covariates: x1-equivalent from the stretch (approximate linear mapping
        # s ~ 0.98 + 0.091 x1, Guy et al. 2010) and the near-peak color-warp tilt
        # (E(B-V)-like amplitude) as the c proxy - documented approximations
        P['x1'] = round((s - 0.98) / 0.091, 3)
        ks = [w[1] for w, p_ in zip(warps0, ph)
              if np.isfinite(w[1]) and abs(p_) < 10.]
        P['c'] = round(float(np.mean(ks)), 4) if ks else None
        draws = []
        for _ in range(n_mc if n_mc > 1 else 0):
            pert = {id(r): r['mag'] + rng.standard_normal() * r['emag']
                    for r in rows}
            _, Bd, Vd, _ = warped_rest_BV(rows, z, mo, mr, t0, s,
                                          mag_lookup=lambda r: pert[id(r)])
            draws.append(np.concatenate([Bd, Vd]))
        i_near = int(np.argmin(np.abs(ph)))
        off = B0[i_near] - float(mr.bandmag('bessellb', 'vega', ph[i_near] / s))
        b_max = float(mr.bandmag('bessellb', 'vega', 0.)) + off
        v_max = float(mr.bandmag('bessellv', 'vega', 0.)) + off
        dm15 = float(mr.bandmag('bessellb', 'vega', 15. / s)
                     - mr.bandmag('bessellb', 'vega', 0.))
    else:
        raise ValueError(f'unknown engine {engine!r}')

    good = np.isfinite(B0) & np.isfinite(V0)
    D = np.array([d for d in draws if np.isfinite(d).all()]) if draws else \
        np.zeros((0, 2 * len(B0)))
    n = len(B0)
    cov = np.cov(D.T) if len(D) >= 8 else np.diag(np.full(2 * n, 0.05 ** 2))
    P.update(ok=True, phases=ph[good], B=B0[good], V=V0[good],
             eB=np.sqrt(np.diag(cov)[:n])[good],
             eV=np.sqrt(np.diag(cov)[n:])[good],
             cov=cov[np.ix_(np.concatenate([good, good]),
                            np.concatenate([good, good]))],
             t0=float(t0), b_max=float(b_max), v_max=float(v_max),
             dm15_rest=float(dm15) if np.isfinite(dm15) else None,
             n_epochs=int(good.sum()))
    return P


# ------------------------------------------------- template-insensitivity test
def insensitivity_variants(rows, z, base, seed=0):
    """Cookbook section 8.1: recompute the synthesis with deliberately offset
    templates (x1 +/- 0.5, c +/- 0.05) and with the Hsiao template, each morphed to
    the same observed colors. Returns {variant: (phases, B, V)} for the caller to
    turn into a B_BV0.6 spread (K_systematic)."""
    _require()
    rows = [r for r in rows if r.get('bandpass')]
    t0 = base['t0']
    par = base.get('params', {})
    x1 = par.get('x1', 0.); c = par.get('c', 0.); x0 = par.get('x0', 1e-10)
    out = {}
    variants = [('x1+0.5', x1 + 0.5, c), ('x1-0.5', x1 - 0.5, c),
                ('c+0.05', x1, c + 0.05), ('c-0.05', x1, c - 0.05)]
    for name, xx1, cc in variants:
        try:
            mo, mr = _salt_pair(z, t0, float(np.clip(xx1, -4, 4)),
                                float(np.clip(cc, -0.4, 1.0)), x0)
            out[name] = warped_rest_BV(rows, z, mo, mr, t0)[:3]
        except Exception:
            continue
    try:
        mo, mr = _hsiao_pair(z, t0)
        mo.set(t0=t0)
        s = base.get('params', {}).get('stretch', 1.0)
        out['hsiao'] = warped_rest_BV(rows, z, mo, mr, t0, s)[:3]
    except Exception:
        pass
    return out
