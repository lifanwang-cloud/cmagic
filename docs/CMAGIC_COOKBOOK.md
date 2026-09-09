# The CMAGIC Cookbook

**C**olor–**MAG**nitude **I**ntercept **C**alibration: Type Ia supernova distances from the
magnitude at a fixed color.

**Contents** — Part I (low redshift): §1 the hypothesis; §2 the linear region; §3 modes
L/Q/R; §4 K-corrections; §5 extinction; §6 standardization; §7 failure gates and
provenance. Part II (high redshift): §8 cross-filter template synthesis; §8.1 the
iterative template fit and the insensitivity of B_BV0.6; §9 Mode S (sparse, full
covariance); §10 validation (SDSS-II) and application scope; §10.1 the SDSS-II worked
example (numbers).

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

### 6.1 Sample-level fitted corrections (v0.3)

The PI's rule, verbatim: **"the correct way is to always fit the x1 and c correction so
that any correlation is removed."** Whenever a SAMPLE of CMAGIC distances is assembled,
the standardization of §6 gains fitted shape and color terms:

```
resid_i = M0 + delta_S·[mode_i = S] + alpha_C·x1_i + beta_C·c_i
```

fit by robust weighted least squares against the Hubble residuals
(`cmagic.standardize_sample`; covariates from the object's own engine — the SALT3-NIR
iterative (x1, c), or the Hsiao stretch mapped through s ≈ 0.98 + 0.091 x1 with the
near-peak color-warp tilt as the c proxy — or supplied externally, e.g. from an
independent SALT3 run). The mode-offset term delta_S repairs the L-versus-S zero-point
split found in the §10.1 comparison; coefficient covariance is propagated into the
corrected errors, and per-mode errors are inflated so the post-fit chi²/dof is unity per
mode. The post-fit weighted residual correlations against x1 and c are zero by
construction — that is the acceptance criterion, and the test suite asserts it.

**Precedence (no double color correction).** When beta_C is fitted, the internal W03
Eqs 7–9 host-extinction correction is backed out of the base magnitudes — the fitted
color term empirically absorbs what the internal E_host estimator measures (and, on
synthesized high-z peaks, what it misses). E_host remains reported for diagnostics, but
its correction applies only when the user excludes 'c' from the covariates.

**Single-object caveat.** A single supernova is never corrected by fitted terms:
`distance()` output is estimator-pure, and the fitted corrections are a sample-level
operation — exactly as SALT3's alpha/beta are training products, not per-object physics.

**SDSS-II numbers (v0.3, 28-object fit set; 6 objects newly gated by the
boundary-hit rule, snids 1794/2017/2031/2440/2635/2992, and 2030 excluded for its
bound-hit external covariate):** M0 = +0.142 ± 0.037, delta_S = −0.101 ± 0.062,
alpha_C = −0.055 ± 0.038, beta_C = +0.953 ± 0.348. Pre → post: rms 0.246 → 0.213
(mode L 0.181 → 0.165, mode S 0.292 → 0.263); chi²/dof 2.85 → 2.09, and 1.0 per mode
after inflation (factors L 1.28, S 1.63). Post-fit weighted correlations: r(x1) = 0.00,
r(c) = 0.00 (unweighted −0.06, +0.01).

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

---

# Part II — CMAGIC at High Redshift

At z ≲ 0.1, single-filter K-corrections suffice and Part I applies as written. At higher
redshift the observer-frame filters no longer sample the rest-frame B and V bands at all:
**cross-filter K-corrections** are required — effectively interpolations of the observed
multi-band photometry, guided by a spectral template, from the observer frame to the
rest-frame B and V light curves that CMAGIC needs. This part formulates that extension,
its sparse-data estimator, and its validation.

## 8. The cross-filter problem and the template-synthesis solution

For a supernova at redshift z observed in filters {X₁ … X_n} (Rubin/LSST ugrizy, Roman
WFI, JWST NIRCam, or the standard UBVRIZ/ugriz/JHK sets), rest-frame B at rest phase p
samples the SED near 4400(1+z) Å in the observer frame — between, or beyond, the observed
bands. The classical per-band cross-filter K-correction (Kim, Goobar & Perlmutter 1996)
transfers each observed magnitude to a chosen rest band. We adopt the equivalent but more
robust **template-synthesis** formulation, which the pipeline implements end to end:

1. **Rest-frame clock.** Phases are computed as p = (t − t₀)/(1+z): cosmological time
   dilation is applied before anything else, and the CMAGIC window of §2, the decline rate
   Δm₁₅, and the color-curve turnback are all defined in *rest-frame* days.
2. **Template.** A spectral time series F_T(p, λ) of a normal Ia — either the **Hsiao
   template** (the simple default the method needs nothing fancier than) or the
   **SALT3-NIR** model (Pierel et al. 2022), whose wavelength coverage extends into the
   near-infrared and is preferred when the observed bands redshift beyond ~8500 Å rest
   (Roman and JWST at z ≳ 1, where rest-frame B–V is sampled by observer-frame YJH).
3. **Warp to the data.** The template is redshifted — F_T(p, λ(1+z))/(1+z), i.e. wavelengths
   stretched and the SED dimmed per bandpass convention — and constrained by the observed
   photometry: in Hsiao mode, a stretch (from the fitted rest-frame decline) plus a smooth
   per-epoch color warp (a grey term and an extinction-law-shaped tilt) fitted so the
   template's synthetic magnitudes reproduce every observed band at every epoch; in
   SALT3-NIR mode, the standard (t₀, x₀, x₁, c) fit. The warp is smooth in wavelength by
   construction, so it interpolates *between* the observed bands without inventing spectral
   features.
4. **Synthesize rest-frame B, V.** At each observation epoch, rest-frame B and V magnitudes
   are synthesized from the warped template through the Bessell B and V transmissions. The
   output is exactly what Part I consumes: nightly B and B−V, now in the rest frame, with
   the K-correction, the (1+z) dilation, and the wavelength stretch all inside the synthesis.
5. **Errors and covariances by Monte Carlo.** The observed photometry is perturbed by its
   errors (default 200 realizations), the warp and synthesis are redone per realization, and
   the empirical covariance of the synthesized (B_i, V_i) — including the B–V covariance
   within an epoch and correlations between epochs induced by the shared warp — is
   propagated into the CMAGIC fit. A template-choice systematic (Hsiao versus SALT3-NIR
   difference, in the bands actually used) is reported alongside.

Filter transmissions: the package ships with (or registers on demand from `sncosmo`) the
Bessell UBVRI, SDSS/LSST ugrizy, Roman WFI, and 2MASS JHK curves; any other filter — JWST
NIRCam included — is supplied as a two-column (wavelength, transmission) file via
`cmagic.highz.register_filter()`. The inputs per supernova are then exactly as the design
requires: **the light curves, the filter transmissions, and the redshift.**

## 9. Sparse light curves: the fixed-slope estimator with full covariance

High-z light curves are often sparse: one or two epochs may land in the CMAGIC linear
window. The linear fit of Mode L is then impossible — but the hypothesis of §1 does not
require fitting the slope, because the slope is a population constant. **Mode S (sparse)**:
with β fixed at 1.94, every single point in the window is an estimate of the candle,

```
m_i = B_i − β (c_i − 0.6),   c_i ≡ B_i − V_i
    = (1−β) B_i + β V_i + 0.6 β
```

with variance following from the epoch's (B, V) covariance matrix — note that B appears on
both axes of the diagram, so the B–V covariance term is *not* optional:

```
Var(m_i) = (1−β)² Var(B_i) + β² Var(V_i) + 2β(1−β) Cov(B_i, V_i)
```

Multiple sparse points are combined by generalized least squares with the full covariance

```
V_ij = Cov(m_i, m_j)  +  σ_β² (c_i − 0.6)(c_j − 0.6)  +  δ_ij σ_int²
```

where the first term carries the synthesis-induced inter-epoch correlations (§8.5), the
second propagates the population slope uncertainty σ_β = 0.16 — *fully correlated* across
the epochs of one supernova, and vanishing for points at the fiducial color — and σ_int is
the intrinsic scatter of the locus about the line (set from the low-z training sample).
The slope term makes explicit what sparse data cost: a single point at c = 0.9 pays
0.16 × 0.3 ≈ 0.05 mag of slope systematic; a point near c = 0.6 pays almost nothing —
so the *placement* of the surviving epochs, not just their number, sets the error. Mode S
reports n_window, the color leverage |c̄ − 0.6|, and the slope-systematic share of the
error budget, and it is gated exactly like the other modes (`window_empty` when nothing
survives; no silent degradation).

## 10. Validation and application scope

**Validation: the SDSS-II supernova survey** (ugriz light curves, z ≈ 0.05–0.4, published
photometry and redshifts) — deep enough that cross-filter synthesis is genuinely exercised,
low enough that results can be checked against the standard-candle expectation. The
validation applies the full chain blind (template synthesis → mode selection, with Mode S
carrying the sparse tail) and reports the Hubble diagram against ΛCDM, the residual rms by
redshift bin and by mode, the Hsiao-versus-SALT3-NIR template systematic, and the failure
census with named gates. Passing criteria are honest rather than cosmetic: residual rms
consistent with the low-z dispersion of §6 plus the propagated sparse-mode errors, and no
mode- or z-dependent bias exceeding its error.

**Application scope.** With filters registered, the identical call serves Rubin (ugrizy to
z ~ 0.4 in rest-B–V; deeper with y), Roman WFI (rest-frame B–V to z ≳ 2), and JWST NIRCam
follow-up photometry. The cookbook's rule of thumb: choose observer bands so that
4400(1+z) Å and 5500(1+z) Å each fall *inside* the covered wavelength range — synthesis
interpolates well and extrapolates badly, and the pipeline gates (`rest_band_uncovered`)
when a rest band would require extrapolating the template beyond the observed filters.

### 8.1 The iterative template fit, and why B_BV0.6 is stable against template details

The SALT3-NIR parameters cannot be fixed in one pass, because the K-correction and the
light-curve fit depend on each other. The engine iterates (PI prescription):

1. make a first guess of (x₁, c) — from the raw observed decline and color;
2. compute the K-correction / rest-frame synthesis with the template at that (x₁, c);
3. fit the light-curve template to the synthesized rest-frame light curves to derive new
   (x₁, c);
4. return to step 2; exit when (x₁, c) converge (default: |Δx₁| < 0.05 and |Δc| < 0.01,
   typically 2–4 iterations; the pipeline gates with `kcorr_no_convergence` otherwise).

Convergence is benign for a reason worth stating plainly: **the K-correction is essentially
a color correction.** Once the template — *any* reasonable Ia template — is morphed to match
the observed colors, the synthesized B magnitude at the fiducial color B−V = 0.6 is a stable
quantity, insensitive to the residual details of the template. The fiducial-color evaluation
point is what buys this: template imperfections move points *along* the color–magnitude
locus far more than they move the locus itself, and reading the locus at fixed color
cancels the along-track freedom.

This is not an assumption but a measurable property, and the pipeline measures it: the
**template-insensitivity test** recomputes the synthesis with templates of deliberately
different (x₁, c) — and with the Hsiao template in place of SALT3-NIR — each morphed to the
same observed colors, and reports the spread of the resulting B_BV0.6 as the template
systematic (superseding the simpler Hsiao-versus-SALT3-NIR delta as the quoted systematic;
the pass criterion is that this spread be small against the photometric error, which is the
quantitative form of the stability claim above).


### 10.1 Worked example: the SDSS-II validation run (v0.2.0)

Blind run of `examples/sdss_validation.py` on 60 spectroscopically confirmed SDSS-II
SNe Ia (z 0.05–0.35; public release of Sako et al. 2018 — see the script header for the
exact source URLs), engine `hsiao` primary with the §8.1 insensitivity test per object:

- **36 of 60 yield distances** — modes: L ×17, S ×19. Failure census (all named):
  `no_bracket` ×11 (blue branches never reaching B−V = 0.6 inside the window),
  `rms_gate` ×9, `window_empty` ×3, `too_few_nights` ×1.
- **Hubble residuals about ΛCDM (H₀ = 72, Ωm = 0.3)**: grey zero-point offset
  −0.17 mag (the transfer of the Vega/W03 calibration through synthesis; documented in
  KNOWN_ISSUES, not recalibrated away), rms about the median 0.44 mag overall:
  **mode L 0.24 mag**, mode S 0.56 mag; by redshift bin 0.33 / 0.53 / 0.38
  (0.05–0.15 / 0.15–0.25 / 0.25–0.35).
- **Error-budget honesty**: χ²/dof = 3.9 (L) and 9.4 (S) — the propagated errors
  underestimate the observed scatter by ≈2× (L) and ≈3× (S). The L-mode excess over the
  §6 calibration dispersion plausibly carries SMP photometry systematics and the
  host-extinction estimator applied to synthesized peak magnitudes; the S-mode excess
  concentrates at high color leverage (the §9 warning realized: single points far from
  the fiducial color). **The §10 passing criterion is NOT met at v0.2.0**; the numbers
  are reported as measured.
- **Template-insensitivity (§8.1)**: median B_BV0.6 spread over {x₁ ± 0.5, c ± 0.05,
  Hsiao} = **0.056 mag**, 90th percentile 0.091 — the stability claim holds at the
  few-hundredths level, and per-object spreads are carried as `K_systematic`.

#### SALT3 comparison (same 60-SN subset; examples/sdss_salt3_comparison.py)

SALT3 (sncosmo; Tripp with alpha = 0.14, beta = 3.1, M_B = −19.36; validity cuts
|x₁| < 3, |c| < 0.3) fits **56 of 60**; the overlap with CMAGIC's 36 is **34 objects**.
On the overlap, each method's own grey offset removed:

| method | N | grey offset | rms | χ²/dof | rms (mode-L subset) |
|---|---|---|---|---|---|
| SALT3 | 34 | +0.34 | **0.142** | 11.2* | 0.093 |
| CMAGIC | 34 | −0.20 | 0.315 | 4.3 | 0.241 |

*SALT3's χ²/dof uses its formal fit errors with no intrinsic-scatter term; adding the
usual σ_int ≈ 0.1 brings it to ≈1.3.

**Residual correlation: Pearson r = −0.01 ± 0.18** — consistent with zero, formally even
below the ~0.15 of W03/WS06 at low z, though CMAGIC's larger independent noise dilutes r
at this depth: the data are consistent with the two methods sharing little beyond the
photometry. Five of 34 objects show |Δμ| > 3σ (all with CMAGIC fainter than SALT3 at
low z, K_systematic ≤ 0.08 — so not template-choice); and the Δμ panel exposes a
**mode-dependent relative offset**: mean μ_CMAGIC − μ_SALT3 = −0.37 (mode L) versus
−0.71 (mode S), an internal L-versus-S inconsistency of ≈ 0.34 mag that is invisible in
CMAGIC's own Hubble scatter but obvious against the common SALT3 reference.

**Honest paragraph.** On sparse, cross-filter SDSS-II data at v0.2.0, CMAGIC is not
competitive with SALT3 as a distance estimator: 2–3× the residual rms on the same
objects, a lower yield (36 vs 56), and an unresolved L-versus-S zero-point split —
SALT3 was trained end-to-end for exactly this regime, while the CMAGIC chain here
transfers a low-z, Vega-calibrated, two-band construction through template synthesis.
What the comparison buys is what W03/WS06 originally claimed and the r ≈ 0 confirms:
the two estimators are statistically independent, so CMAGIC remains valuable as a
cross-check and systematics probe (its per-object K_systematic and named gates localize
what SALT3 absorbs silently), and the identified defects — the sparse-mode error model,
the synthesized-peak extinction estimator, and the mode-dependent zero point — are
concrete v0.3 targets rather than fundamental limits.


#### v0.3 update to the SALT3 comparison

With the boundary-hit gate (6 stretch-at-bound objects removed, including two of the
former >3σ outliers) and the §6.1 fitted corrections applied, the CMAGIC sample tightens
from rms 0.32 (v0.2 overlap) to **0.21 mag** overall (0.165 mode L), the L-versus-S
zero-point split is absorbed into delta_S = −0.10 ± 0.06, and the per-mode error budgets
close (chi²/dof = 1 by construction of the inflation). SALT3 remains tighter (0.14/0.09);
the gap after correction is ≈1.5× rather than ≈2–3×.
