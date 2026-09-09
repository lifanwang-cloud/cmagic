# The CMAGIC Method

**C**olor–**MAG**nitude **I**ntercept **C**alibration: Type Ia supernova distances from the
magnitude at a fixed color.

This document describes the method implemented by this package, in enough detail to use it
critically. It is based on the two founding papers —
Wang, Goldhaber, Aldering & Perlmutter (2003, ApJ 590, 944; hereafter **W03**) and
Wang, Strovink, et al. (2006, ApJ 641, 50; hereafter **WS06**) — together with a reformulation
(three measurement modes and explicit failure gates) developed for archival applications.
Users who only want a distance can read §1 and then the README; the pipeline makes every choice
described here automatically and reports what it chose.

---

## 1. The physical hypothesis

After maximum light, a Type Ia supernova's B-band magnitude and its B−V color evolve
*together*: on the color–magnitude diagram (B versus B−V), the supernova travels a remarkably
straight line for much of the first month past maximum, as the ejecta expand and cool. The
CMAGIC hypothesis is:

> **The B magnitude at a fixed color, B−V = 0.6 mag — written B_BV0.6 — is a calibratable
> standard candle**, in exactly the sense that the peak magnitude is the standard candle of
> light-curve fitters such as SALT3.

Why measure brightness at fixed *color* rather than at fixed *time* (peak)?

1. **The locus is straight and slowly traversed.** The linear color–magnitude relation
   (slope β_BV = 1.94 ± 0.16 after K-corrections; W03) is measured from many nights of data,
   not from the single moment of peak, so B_BV0.6 is statistically robust and does not require
   dense coverage of maximum light.
2. **Reduced extinction sensitivity.** Dust moves a supernova *along* a vector of slope
   R_B ≈ 4.1 in the same diagram. The extinction correction to B_BV0.6 is therefore
   (R_B − β_BV)·E(B−V) ≈ 2.2·E rather than R_B·E ≈ 4.1·E — roughly half the sensitivity of a
   peak-magnitude method to reddening errors (W03).
3. **Smaller intrinsic dispersion.** In the low-extinction calibration sample, B_BV0.6
   standardizes to σ ≈ 0.07–0.10 mag, versus ≈ 0.10 mag for B_max (W03, Table 2).
4. **Independence.** The residuals of CMAGIC and peak-magnitude distances correlate only
   weakly (Pearson ≈ 0.15; WS06), so the two are usefully independent distance estimators
   from the same photometry.

The distance chain is then conceptually identical to any standard-candle method:

```
B_BV0.6  →  K-correction  →  extinction correction  →  M_BV(Δm15) calibration  →  μ, D
```

---

## 2. The color–magnitude diagram and the linear region

Plot nightly B against nightly B−V. The locus has four parts:

- a **pre-maximum branch** (distinct; hysteresis relative to the post-max branch),
- a **post-maximum shoulder** in the first days after B maximum (curved; W03 exclude a
  two-day gap around the transition),
- the **linear region** — the CMAGIC branch proper, typically spanning B−V ≈ 0.2–1.0 for an
  unreddened object,
- the **color turnback**, when B−V reaches its maximum (t ≈ +30 d for normal decliners) and
  turns blueward; beyond it the relation is no longer the same line.

**The fiducial window is decline-rate-scaled** (WS06):

| subtype | window (days past B max) |
|---|---|
| normal | [ 7/Δm₁₅ , 30/Δm₁₅ ] |
| 91T-like | [ 10/Δm₁₅ , 29/Δm₁₅ ] |

For Δm₁₅ = 1.1 this is +6.4…+27.3 d. The package intersects this window with an object-specific
safety exit at t(B−V max) − 3 d, computed from the object's own smoothed color curve, which
protects objects with unusually early or late color maxima. Within the window, nights are
combined by weighted averaging; no residual-based point rejection is applied beyond a
conservative 3σ photometric clip.

Two practical caveats established on archival data:

- **The window inherits Δm₁₅ quality.** A poorly measured decline rate (e.g., a light-curve
  fit with the stretch parameter at its bound) produces a nonsense window; the pipeline gates
  on Δm₁₅ provenance.
- **Slow decliners.** For Δm₁₅ ≲ 0.9 the 7/Δm₁₅ entry can admit the tail of the post-max
  shoulder, mildly biasing a free-slope fit. The calibrated chain (fixed slope, §3) is robust
  to this at the few-percent level in distance; the free slope is reported as a diagnostic,
  not used in the distance.

---

## 3. The three measurement modes

The goal is always the same number — the magnitude at B−V = 0.6 — but real light curves come
in three conditions. The pipeline chooses the mode automatically and records the choice.

### Mode L — linear (the classical method)

When the linear region is well populated, fit a weighted straight line to B versus B−V in the
window and evaluate it at B−V = 0.6:

```
B(c) = B_BV0.6 + β_BV (c − 0.6)
```

The **calibrated chain fixes β_BV = 1.94** (the K-corrected population slope; W03) — this is
what the M_BV calibration of §6 was constructed with. The free slope β_free is also fit and
reported with its error as a *diagnostic*: values far from 1.94 on clean, single-source
photometry indicate either peculiarity or time-variable extinction (see §7's SN 2006X note).
A free-slope value is considered valid only if the fit has ≥ 5 nights, color span ≥ 0.25 mag,
no internal color gap > 0.4 mag, rms ≤ 0.15 mag, and β_free within the WS06 validity cut
[1.5, 2.5].

### Mode Q — quadratic interpolation

When the linear region is not obvious — sparse windows, curved (parabolic) loci, and
91bg-like objects for which W03 showed the linear model fails — B_BV0.6 is still perfectly
well defined *as the magnitude where the locus crosses B−V = 0.6*. Mode Q fits

```
B(c) = a₀ + a₁ (c − 0.6) + a₂ (c − 0.6)²
```

over the post-maximum branch and reports a₀ ± σ. **Interpolation only**: the observed colors
must bracket 0.6 (or approach it within 0.05 mag), otherwise the mode fails with a named gate
(`no_bracket`). For 91bg-like objects (Δm₁₅ ≥ 1.7) the *measurement* is delivered but the
*standardization* is flagged `uncalibrated`: the M_BV(Δm₁₅) calibration of §6 does not extend
to them.

### Mode R — heavily reddened objects

When extinction pushes the entire observed branch redward of B−V = 0.6, extrapolating the fit
back to 0.6 is forbidden — extrapolation across a magnitude or more of color is precisely the
systematic that corrupts naive fits. Instead (PI prescription):

1. Estimate the reddening from the peak color: **E_guess = B_max − V_max** (the intrinsic
   peak color of a normal Ia is ≈ 0; optionally refine by subtracting the intrinsic color
   locus ℰ₀(Δm₁₅) of Eq. W03-8 below).
2. Measure the magnitude at the **shifted target color c\* = 0.6 + E_guess**, by Mode L or
   Mode Q, strictly by interpolation within the observed branch. Record (m\*, c\*, E_guess,
   β_used).
3. When the true reddening E_true is later determined (from §5, from external data, or from a
   dust-law analysis), apply the deferred correction:

```
B_BV0.6 = m* − β_used (c* − 0.6 − E_true) − 0.6 β_used − R_B E_true + 0.6 β_used
        = m* − β_used (c* − 0.6) + (β_used − R_B) E_true
```

   For a linear locus this is algebraically identical to the classical chain (measure at
   observed 0.6, correct by (R_B − β)E) — Mode R changes *where the locus is evaluated*, never
   the physics — but it avoids the extrapolation entirely. The host dust law R_B is
   configurable per object (default 3.1; e.g., SN 2006X requires its known anomalous
   R_B ≈ 2.5, and Mode R makes the dust law the *only* open parameter for such objects).

---

## 4. K-corrections (the WS06 scheme)

Five raw outputs receive K-corrections: B_max, V_max, Δm₁₅, B_BV0.6, and β_BV. Following
WS06 (their Appendix D): each correction is **linear in the supernova's raw color** E, with an
intercept μᵏ(z) and slope ψᵏ(z) that are cubic polynomials in redshift, calibrated from a
spectrophotometric library of reddened/unreddened light-curve pairs. The package evaluates the
polynomials at the object's z and raw color. At z ≤ 0.03 the corrections are ≤ 0.01 mag; they
are always applied and always quoted. (Reproduction note: the distance-relevant K(B_BV)
reproduces the WS06 worked example to ~1 mmag; the published per-band table entries for
B_max/V_max include light-curve-refit effects beyond the tabulated polynomials, an ambiguity
≤ 0.01 mag in μ at these redshifts, documented in the validation suite.)

---

## 5. Extinction

Milky-Way extinction is removed first (Cardelli law, R_B = 4.15, with E(B−V)_MW from the
Schlafly & Finkbeiner maps). The **host** component uses CMAGIC's internal estimator (W03
Eqs. 7–9): the color-excess measure

```
ℰ = (B_max − B_BV0.6) / β_BV                                (W03 Eq. 7)
ℰ₀ = (−0.118 ± 0.013) + (0.249 ± 0.043)(Δm₁₅ − 1.1)         (W03 Eq. 8)
E_host = ℰ − ℰ₀                                              (W03 Eq. 9)
```

with ℰ₀ the intrinsic (dust-free) locus of the same quantity as a function of decline rate.
The host correction to the candle is A_BV = (R_B − β_BV)·E_host with the host-law
R_B = 3.1 by default (override per object when the dust law is known to be anomalous).
E_host absorbs *color calibration errors* of the photometry as if they were dust — the reason
the pipeline enforces single-source photometry (§7).

---

## 6. Standardization and the distance

The absolute calibration is the W03 Table 3 linear law M_BV(Δm₁₅) = a + b·Δm₁₅, constructed
on color-restricted subsamples (B_max − V_max ≤ 0.05 / 0.2 / 0.5) with dispersions
0.07–0.16 mag, on the H₀ = 65 scale (rescaled explicitly; the H₀ dependence is a single
coherent 5 log(H₀/65) shift). The WS06 refinement — a bilinear decline-rate law with a kink at
Δm₁₅ = 1.1, which beats a quadratic at 43:1 odds — is available as an option (zero point
anchored to W03 Table 3 at Δm₁₅ = 1.1). Then

```
μ = B_BV0.6(corrected) − M_BV(Δm₁₅),      D = 10^{(μ−25)/5} Mpc.
```

Error budget per object: the fit error on B_BV0.6 (typically 0.02–0.06 mag on good archival
photometry), the calibration dispersion (0.07–0.16 mag, dominant), the extinction-estimator
error propagated through (R_B − β), and the K-correction ambiguity (≤ 0.01 mag at z ≤ 0.03).

---

## 7. Failure semantics and provenance (what makes this pipeline safe to use blind)

Every fit either returns a distance **or a named gate**, never a silently degraded number:

| gate | meaning |
|---|---|
| `window_empty` | the decline-rate window contains no valid nights (e.g., 91bg-like: the color maximum arrives so early that no linear branch exists — that *is* the physical result) |
| `too_few_nights` | < 4 nights in the window |
| `no_bracket` | Mode Q/R interpolation impossible: the branch never crosses the target color |
| `rms_gate` / `beta_gate` | Mode L quality gates failed |
| `dm15_provenance` | decline rate unmeasured or from a bound-hit fit: no trustworthy window |
| `photometry_defect` | upstream data flagged (duplicated columns, cross-instrument artifacts) |

Provenance rules learned from archival forensics, enforced automatically:

- **One photometric source per fit** (never mix observatories in a band; preference order
  CSP > CfA > KAIT/LOSS > other; Swift/UVOT optical bands are never used). Mixed compilations
  inflate the CMAGIC rms three- to five-fold and bias the internal extinction estimator.
- **De-duplication**: rows duplicating the same publication via different aggregators are
  dropped by bibcode matching.
- Every result carries: source used, nights selected, phase and color spans, mode, β_free ± σ,
  rms, all gates evaluated, and the diagnostic panel (locus, window boundaries, fit, color-curve
  inset, annotation block).

**A diagnostic worth knowing about**: on clean single-source photometry, a *valid* free slope
significantly different from 1.94 carries physics. SN 2006X — with documented time-variable
circumstellar absorption — shows β_free = 2.25–2.28 (rms at calibration quality), intermediate
between the intrinsic slope and its anomalous dust vector R_B ≈ 2.5: time-variable extinction
slides points along the dust vector. The free slope is CMAGIC's built-in dust-variability alarm.

---

## References

- Wang, L., Goldhaber, G., Aldering, G., & Perlmutter, S. 2003, ApJ, 590, 944
- Wang, L., Strovink, M., Conley, A., Goldhaber, G., Kowalski, M., Perlmutter, S., &
  Siegrist, J. 2006, ApJ, 641, 50
- Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
- Schlafly, E. F., & Finkbeiner, D. P. 2011, ApJ, 737, 103
