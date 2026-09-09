#!/usr/bin/env python
"""Per-object CMAGIC diagnostic panels for the DES-SN5YR validation subset.
Identical chain to des_validation.py --assisted (same cache, filters, seeds,
release t_bmax + external dm15), with panel= set.
Writes examples/figures/des_panels/<snid>_<mode|gate>.png."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import des_validation as V
import cmagic

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures", "des_panels")
os.makedirs(OUT, exist_ok=True)

meta, phot = V.load_cache()
only = sys.argv[1:] or None
if only:
    meta = [m for m in meta if m['snid'] in only]
count = 0
for k, m in enumerate(meta):
    rows = phot[m['snid']]
    table = dict(mjd=np.array([float(r['mjd']) for r in rows]),
                 band=np.array([r['band'] for r in rows]),
                 mag=np.array([float(r['mag']) for r in rows]),
                 emag=np.array([float(r['emag']) for r in rows]))
    tag = "err"
    try:
        res = cmagic.distance(
            photometry=table, z=m['z'], filters=V.FILTERS, engine='auto',
            ebv_mw=m['mwebv'], n_mc=100, seed=42, h0=V.H0,
            t_bmax=m['pkmjd'],
            dm15=V.dm15_from_release_salt3(m['x1_des'], m['c_des']),
            panel=os.path.join(OUT, f"tmp_{m['snid']}.png"))
        tag = res.mode if res.failed is None else f"FAIL-{res.failed}"
    except Exception as e:
        print(m['snid'], "exception:", str(e)[:60], flush=True)
        continue
    src = os.path.join(OUT, f"tmp_{m['snid']}.png")
    dst = os.path.join(OUT, f"{m['snid']}_{tag}.png")
    if os.path.exists(src):
        os.replace(src, dst)
        count += 1
    print(f"[{k+1}/{len(meta)}] {m['snid']} -> {tag}", flush=True)
print(f"{count} panels in {OUT}")
