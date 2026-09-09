# cmagic

**Type Ia supernova distances from the color–magnitude intercept (CMAGIC).**

CMAGIC measures the B-band magnitude of a Type Ia supernova at a fixed color
(B−V = 0.6 mag) — a standard candle with roughly half the extinction sensitivity of
peak-magnitude methods and an intrinsic dispersion of ≈ 0.08–0.16 mag
(Wang et al. 2003, ApJ 590, 944; Wang, Strovink et al. 2006, ApJ 641, 50).

**You do not need to know how CMAGIC works to get a distance.** Give the pipeline a B and V
light curve and a redshift; it selects the fitting window, chooses the measurement mode,
applies K-corrections and extinction corrections, standardizes, and returns a distance with
its error — or an explicit, named reason why your data cannot yield one. The full method is
documented in [docs/METHOD.md](docs/METHOD.md).

## Installation

```bash
pip install git+https://github.com/lifanwang-cloud/cmagic.git
# or, from a clone:
git clone https://github.com/lifanwang-cloud/cmagic.git
cd cmagic
pip install .
```

Requirements: Python ≥ 3.10, `numpy`, `matplotlib` (for diagnostic panels only).
No other dependencies.

## Quick start: a distance in five lines

```python
import cmagic

result = cmagic.distance(
    photometry="my_sn_photometry.csv",   # columns: mjd, band, mag, emag  (B and V rows)
    z=0.0055,                            # heliocentric redshift
)
print(result.D_mpc, "+/-", result.eD_mpc, "Mpc   [mode:", result.mode, "]")
```

That is the entire required interface. `photometry` may be a CSV path, a numpy structured
array, or a dict of arrays. Everything else is optional:

```python
result = cmagic.distance(
    photometry=phot,
    z=0.0055,
    t_bmax=53286.5,        # MJD of B maximum, if you know it (else measured from the data)
    dm15=1.17,             # decline rate, if you know it (else measured from the data)
    ebv_mw=0.023,          # Milky-Way E(B-V) (else 0; pass your Schlafly & Finkbeiner value)
    rb_host=3.1,           # host dust law (override for known-anomalous objects)
    h0=72.0,               # output distance scale
    panel="my_sn.png",     # write the diagnostic figure
)
```

## What you get back

`result` is a dataclass with, among others:

| field | meaning |
|---|---|
| `D_mpc`, `eD_mpc`, `mu` | the distance, its error, the distance modulus |
| `B_BV06`, `eB_BV06` | the standard-candle magnitude at B−V = 0.6 |
| `mode` | `"L"` (linear), `"Q"` (quadratic), or `"R"` (reddened-target) — chosen automatically |
| `beta_free`, `ebeta_free` | the free color–magnitude slope (diagnostic; 1.94 is the population value) |
| `E_host`, `E_guess` | internal reddening estimates |
| `window`, `n_nights`, `rms` | the fitted region and fit quality |
| `gates` | every validity gate evaluated, pass/fail |
| `flags` | e.g. `uncalibrated` for 91bg-like objects (measurement delivered, standardization outside the calibration) |

If the data cannot support a distance, `result.failed` is the named gate
(`window_empty`, `too_few_nights`, `no_bracket`, `dm15_provenance`, …) — see
[docs/METHOD.md §7](docs/METHOD.md). The pipeline never returns a silently degraded number.

## Three rules about input photometry

1. **One instrument per fit.** Do not mix observatories in a band — mixed compilations
   inflate the fit rms several-fold and bias the extinction estimator. If your table carries
   multiple sources, pass `source_column=` and the pipeline will select the best single one.
2. B and V should be nightly (or denser) from a few days past maximum through
   about a month past maximum; pre-maximum data are used only to locate B_max and the peak color.
3. Magnitudes in the Johnson–Cousins (Vega) system, or CSP natural (pass `system="csp"`).

## Examples and validation

- `examples/demo_sn2004dt.py` — a clean linear-mode distance from published LOSS photometry (Ganeshalingam et al. 2010).
- `examples/demo_reddened.py` — Mode R on a heavily reddened supernova (never extrapolates;
  the host dust law is exposed as the one open parameter).
- `tests/` — the validation suite: the Wang–Strovink K-correction worked example and
  end-to-end regression tests on the bundled demo objects. Run with `pytest`.

## Citing

If you use this package, cite Wang et al. 2003 (ApJ 590, 944) and
Wang, Strovink et al. 2006 (ApJ 641, 50), plus this repository.

## License

MIT (see LICENSE).
