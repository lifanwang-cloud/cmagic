#!/usr/bin/env python
"""Per-object CMAGIC diagnostic panels for the SDSS-II validation subset.
Identical chain to sdss_validation.py (same cache, filters, seeds), with
panel= set. Writes examples/figures/sdss_panels/<snid>_<mode|gate>.png."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sdss_validation as V
import cmagic

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures", "sdss_panels")
os.makedirs(OUT, exist_ok=True)

meta, phot = V.load_cache()
FILTERS = {'g': 'sdssg', 'r': 'sdssr', 'i': 'sdssi', 'z': 'sdssz'}
count = 0
for m in meta:
    rows = phot[m['snid']]
    table = dict(mjd=np.array([float(r['mjd']) for r in rows]),
                 band=np.array([r['band'] for r in rows]),
                 mag=np.array([float(r['mag']) for r in rows]),
                 emag=np.array([float(r['emag']) for r in rows]))
    tag = "err"
    try:
        res = cmagic.distance(photometry=table, z=m['z'], filters=FILTERS,
                              engine='auto', ebv_mw=m['mwebv'], n_mc=120, seed=42,
                              h0=V.H0, panel=os.path.join(OUT, f"tmp_{m['snid']}.png"))
        tag = res.mode if res.failed is None else f"FAIL-{res.failed}"
    except Exception as e:
        print(m['snid'], "exception:", str(e)[:60]); continue
    src = os.path.join(OUT, f"tmp_{m['snid']}.png")
    dst = os.path.join(OUT, f"{m['snid']}_{tag}.png")
    if os.path.exists(src):
        os.replace(src, dst); count += 1
print(f"{count} panels in {OUT}")
