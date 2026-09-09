"""Photometry input, nightly binning, source selection, and peak measurement.
Accepted photometry forms (docs/METHOD.md section 7 provenance rules):
  - CSV path: columns mjd, band, mag, emag (optional: a source column named by
    source_column); '#' comment lines ignored;
  - numpy structured array with those fields;
  - dict of arrays with those keys.
Bands: rows with band 'B' or 'V' are used; everything else is ignored."""
import csv as _csv
import re as _re
import numpy as np


def load_photometry(photometry, source_column=None):
    """Normalize any accepted input into a list of row dicts."""
    rows = []
    if isinstance(photometry, (str, bytes)):
        with open(photometry) as fh:
            lines = [ln for ln in fh if not ln.lstrip().startswith('#')]
        for r in _csv.DictReader(lines):
            rows.append(r)
    elif isinstance(photometry, np.ndarray) and photometry.dtype.names:
        for i in range(len(photometry)):
            rows.append({k: photometry[k][i] for k in photometry.dtype.names})
    elif isinstance(photometry, dict):
        keys = list(photometry)
        n = len(photometry[keys[0]])
        for i in range(n):
            rows.append({k: photometry[k][i] for k in keys})
    else:
        raise TypeError('photometry must be a CSV path, structured array, or dict '
                        'of arrays')
    out = []
    for r in rows:
        try:
            m = float(r['mag']); e = float(r['emag']); t = float(r['mjd'])
        except (KeyError, TypeError, ValueError):
            continue
        b = str(r['band']).strip()
        if b not in ('B', 'V') or not (np.isfinite(m) and np.isfinite(e)):
            continue
        if e > 0.5 or e < 0:
            continue
        out.append(dict(mjd=t, band=b, mag=m, emag=max(e, 0.01),
                        source=str(r.get(source_column, '')) if source_column
                        else ''))
    return out


def nightly(rows, band):
    """Nightly weighted means; within-night median clip; 3-sigma running-median
    photometric clip (the only rejection in the pipeline; METHOD section 2)."""
    rr = sorted([r for r in rows if r['band'] == band], key=lambda r: r['mjd'])
    if not rr:
        return np.array([]), np.array([]), np.array([])
    nights = []
    for r in rr:
        if nights and r['mjd'] - nights[-1][-1]['mjd'] < 0.5:
            nights[-1].append(r)
        else:
            nights.append([r])
    t, m, e = [], [], []
    for n in nights:
        w = np.array([1 / r['emag'] ** 2 for r in n])
        mm = np.array([r['mag'] for r in n])
        if len(n) > 1:
            med = np.median(mm)
            thr = max(0.15, 3 * float(np.median([r['emag'] for r in n])))
            keep = np.abs(mm - med) <= thr
            if keep.sum() == 0:
                keep = np.abs(mm - med) == np.min(np.abs(mm - med))
            mm, w = mm[keep], w[keep]
        t.append(np.mean([r['mjd'] for r in n]))
        m.append(float(np.sum(mm * w[:len(mm)]) / np.sum(w[:len(mm)])))
        e.append(max(1 / np.sqrt(np.sum(w[:len(mm)])), 0.01))
    t, m, e = map(np.array, (t, m, e))
    if len(t) > 6:
        good = np.ones(len(t), bool)
        for i in range(len(t)):
            j = np.argsort(np.abs(t - t[i]))[:7]
            if np.abs(m[i] - np.median(m[j])) > 0.2:
                good[i] = False
        t, m, e = t[good], m[good], e[good]
    return t, m, e


# ---------------- source selection (METHOD section 7) ----------------
_CSP_DR3 = ('154/211', '154..211', '170905146')
_CFA = ('700/331', '700..331', '200/12', '200...12', '331..815')
_LOSS = ('190/418', '190..418', '425.1789')


def bibkey(catalog):
    c = str(catalog).strip()
    if c.startswith('OSC:'):
        m = _re.match(r'OSC:\d{4}([A-Za-z&+.]+?)\.+(\d+)[.]+([A-Za-z]?\d+)[A-Za-z]$', c)
        if m:
            return (m.group(1).replace('.', '').replace('&', '+').upper(),
                    m.group(2), m.group(3).upper())
        m = _re.match(r'OSC:\d{4}arXiv(\d+)', c)
        return ('ARXIV', m.group(1), '') if m else None
    p = c.split('/')
    return (p[1].replace('&', '+').upper(), p[2], p[3].upper()) \
        if len(p) >= 4 and p[0] == 'J' else None


def source_rank(name):
    n = str(name)
    if 'uvot' in n.lower() or 'swift' in n.lower():
        return 99
    if any(k in n for k in _CSP_DR3) or 'csp' in n.lower():
        return 0
    if any(k in n for k in _CFA) or 'cfa' in n.lower():
        return 1
    if any(k in n for k in _LOSS) or 'kait' in n.lower() or 'loss' in n.lower():
        return 2
    return 3


def select_source(rows, check):
    """Best single source whose rows pass check(rows); rows must carry 'source'.
    OSC rows duplicating a VizieR catalog of the same publication are dropped first."""
    viz = {bibkey(r['source']) for r in rows
           if r['source'] and not r['source'].startswith('OSC:') and bibkey(r['source'])}
    groups = {}
    ndup = 0
    for r in rows:
        s = r['source']
        if s.startswith('OSC:') and bibkey(s) in viz:
            ndup += 1
            continue
        groups.setdefault(s, []).append(r)
    ranked = sorted(groups, key=lambda s: (source_rank(s), -len(groups[s])))
    ok = [s for s in ranked if source_rank(s) < 99 and check(groups[s])]
    if not ok:
        return None, None, dict(n_osc_duplicates_dropped=ndup)
    return ok[0], groups[ok[0]], dict(n_osc_duplicates_dropped=ndup)


# ---------------- peak measurement ----------------
def measure_peak(tB, mB, tV, mV, z, t0_guess=None):
    """t_Bmax, B_max, V_max, dm15 from quartic fits near the peak (no extrapolation
    beyond the data). Returns dict or None when the peak is not covered."""
    if len(tB) < 5:
        return None
    if t0_guess is None:
        t0_guess = tB[np.argmin(mB)]
    m = np.abs(tB - t0_guess) < 10.
    if m.sum() < 5 or tB[m].min() > t0_guess - 2 or tB[m].max() < t0_guess + 2:
        return None
    p = np.polyfit(tB[m] - t0_guess, mB[m], 4)
    lo = max(-8., float((tB[m] - t0_guess).min()))
    hi = min(8., float((tB[m] - t0_guess).max()))
    if hi - lo < 3.:
        return None
    xs = np.linspace(lo, hi, 801)
    ys = np.polyval(p, xs)
    i = int(np.argmin(ys))
    if xs[i] < lo + 0.5 or xs[i] > hi - 0.5:
        return None
    t_bmax = float(t0_guess + xs[i]); b_max = float(ys[i])
    t15 = t_bmax + 15. * (1 + z)
    m2 = np.abs(tB - t15) < 7.
    dm15 = np.nan
    if m2.sum() >= 3:
        p2 = np.polyfit(tB[m2] - t15, mB[m2], min(2, m2.sum() - 1))
        dm15 = float(np.polyval(p2, 0.)) - b_max
    v_max = np.nan
    mv = np.abs(tV - t_bmax) < 7.
    if mv.sum() >= 3:
        pv = np.polyfit(tV[mv] - t_bmax, mV[mv], min(3, mv.sum() - 1))
        v_max = float(np.polyval(pv, 0.))
    return dict(t_bmax=t_bmax, b_max=b_max, v_max=v_max, dm15=dm15)
