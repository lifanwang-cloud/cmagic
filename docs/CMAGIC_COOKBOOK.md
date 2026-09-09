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
[1.5, 2.5]. **Per-fit error normalization (v0.3.1, the PI's prescription): the errors of
the CMAGIC fitting parameters must be normalized to χ²/dof = 1 per fit** — after every
weighted fit (Modes L, Q, R, and the Mode S GLS with n ≥ 2, where n points give
dof = n − 1 for the mean), the fit-parameter errors are scaled by
s = max(1, √(χ²/dof)); χ²/dof and s are stored in the result and printed on the
diagnostic panel. A single-point Mode S has no dof: it relies on the sample-level term
and says so in its flags.

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

### 5.1 The two extinction philosophies

Supernova cosmology has always had two ways to handle a red supernova, and the two teams
whose samples discovered the accelerating expansion each used one of them. The
**empirical-extinction route** — the original Riess et al. (High-z Team / MLCS) analyses —
estimates the extinction of each object from its own photometry as a color excess over an
intrinsic locus, then corrects it through a **known extinction law** (an assumed R_B).
The **Tripp route** (Tripp 1998; the standardization of Perlmutter et al.'s SCP lineage,
and of SALT2/SALT3 today) does not ask *why* an object is red: it fits one linear color
coefficient β on the whole sample and lets the regression decide the correction.

What each buys, and what each risks:

- The **empirical route** is physically explicit and per-object, and it preserves
  whatever signal remains in the residuals. Its systematic is the assumed law:
  SN 2006X's anomalous dust (R_B ≈ 2.48) corrected with the standard 3.1 produces a
  0.82 mag distance error from the law alone (§3, Mode R).
- The **Tripp route** is self-calibrating and agnostic, but its fitted β is a
  variance-weighted compromise between two physically different slopes: the dust vector
  (R_B ≈ 4.1) and the intrinsic color–luminosity slope. One β is correct only if the
  dust-to-intrinsic mixture of the sample's colors is universal. When selection changes
  that mixture with redshift — a magnitude-limited survey keeps losing dusty objects
  toward high z — the correct effective β changes along the redshift axis while the
  applied one does not, and the mismatch leaks into the Hubble diagram. That is the
  mechanism behind the SALT3 redshift drift dissected in §10.3, where the sample-fitted
  β = 2.41 ± 0.20 sits well below both the fiducial 3.1 and the dust value — the
  signature of a blended color.

CMAGIC's position in this dichotomy is distinctive on both counts:

1. **The pipeline carries both routes.** The internal W03 estimator above (E_host from
   ℰ − ℰ₀, corrected through the known law as (R_B − β_BV)·E_host) is the empirical
   route, per object. The sample-level fitted β_C·c of §6.1 is the Tripp route. The
   precedence rule of §6.1 selects which one is in force — never both — so the user
   always knows which assumption they are buying: a known dust law, or a universal
   color–luminosity relation. Disagreement between the two routes on a single object is
   itself a diagnostic; it is what flagged SN 2006X.
2. **The dust lever arm is halved by construction.** Extinction enters B_BV0.6 as
   (R_B − β_BV)·E ≈ 2.2·E rather than R_B·E ≈ 4.1·E, so whichever route errs, the
   damage to a CMAGIC distance is roughly half of what the same error does to a
   peak-magnitude method. This is also why the fitted color term on CMAGIC residuals is
   small (β_C = +0.94 ± 0.35 on the SDSS-II sample, §6.1): most of the color
   sensitivity is removed before standardization ever sees it.

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
its correction applies only when the user excludes 'c' from the covariates. This toggle
is exactly the choice between the two extinction philosophies of §5.1: covariates with
'c' = the Tripp route, without = the empirical-extinction route.

**Single-object caveat.** A single supernova is never corrected by fitted terms:
`distance()` output is estimator-pure, and the fitted corrections are a sample-level
operation — exactly as SALT3's alpha/beta are training products, not per-object physics.

**SDSS-II numbers (v0.3.1 per-fit-normalized errors; 28-object fit set; 6 objects gated
by the boundary-hit rule, snids 1794/2017/2031/2440/2635/2992, and 2030 excluded for its
bound-hit external covariate):** M0 = +0.145 ± 0.036, delta_S = −0.101 ± 0.063,
alpha_C = −0.056 ± 0.038, beta_C = +0.941 ± 0.346. Pre → post: rms 0.246 → 0.213
(mode L 0.181 → 0.165, mode S 0.292 → 0.263); chi²/dof 2.81 → 2.05, and 1.0 per mode
after inflation. Post-fit weighted correlations: r(x1) = 0.00, r(c) = 0.00.

**The per-fit normalization and the sample factors (honest accounting).** The v0.3.1
per-fit rescaling barely moves the sample-level inflation (L 1.284 → 1.270,
S 1.634 → 1.624): the per-object distance error is dominated by the calibration
dispersion (0.07–0.16) and the host-estimator term (≈0.09), so even a ×1.5 rescaling of
the fit-parameter error (e.g. snid 5103: χ²/dof = 2.2, eB_BV0.6 0.049 → 0.073) shifts
the total error by only a few thousandths of a magnitude. The residual sample factors
therefore measure genuine population/intrinsic scatter plus synthesis systematics — not
fit inconsistency, which is now normalized away per fit as prescribed.

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
- Tripp, R. 1998, A&A, 331, 815
- Riess, A. G., et al. 1998, AJ, 116, 1009
- Perlmutter, S., et al. 1999, ApJ, 517, 565
- Conley, A., et al. 2006, ApJ, 644, 1 (the blind CMAGIC cosmology at high z)
- DES Collaboration 2024, ApJ, 973, L14; Sánchez, B. O., et al. 2024, ApJ, 975, 5
  (DES-SN5YR release, github.com/des-science/DES-SN5YR)

---

# Part II — CMAGIC at High Redshift

At z ≲ 0.1, single-filter K-corrections suffice and Part I applies as written. At higher
redshift the observer-frame filters no longer sample the rest-frame B and V bands at all:
**cross-filter K-corrections** are required — effectively interpolations of the observed
multi-band photometry, guided by a spectral template, from the observer frame to the
rest-frame B and V light curves that CMAGIC needs. This part formulates that extension,
its sparse-data estimator, and its validation. There is direct precedent: Conley et al.
(2006, ApJ 644, 1) measured (Ω_m, Ω_Λ) from a blind CMAGIC analysis of 21 high-z SNe —
on "data sets not observed in a manner optimized for CMAGIC" — and confirmed the
acceleration through the color–magnitude channel.

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

### 10.2 Scatter-statistics convention

Two scatter statistics are quoted throughout, and every table states which it uses.
The **unweighted rms** is the plain standard deviation of the Hubble residuals; the
**error-weighted rms** is (Σw r²/Σw)^½ with w = 1/σ² and each method's *weighted* grey
offset removed. Fits (the sample standardization of §6.1, and every regression in the
pipeline) are always error-weighted; scatter tables in earlier sections quoted the
unweighted rms unless marked. The convention follows Wang, Strovink et al. (2006): their
§2.3 χ²-rescaling of fit errors is the per-fit √(χ²/dof) normalization of §9.1 here, and
their intrinsic-noise terms (added until χ²/dof ≈ 1) are the sample-level inflation of
§6.1. WS06's Hubble-flow variance also carries a peculiar-velocity term
(5/ln 10)·v_pec/(cz) with v_pec fitted (their 95% upper limit 486 km s⁻¹); at this
subsample's z ≥ 0.06 a 300 km s⁻¹ term contributes only 0.006–0.035 mag against
per-object errors of 0.15–0.27 mag, so it is omitted here — include it for any sample
reaching below z ≈ 0.03. On the v0.3 standardized SDSS-II sample, both are shown in
the Hubble figure: SALT3 0.143 unweighted / 0.133 weighted (n = 56); CMAGIC mode L
0.165 / 0.170 (n = 16); mode S 0.263 / 0.272 (n = 12). With the fitted corrections and
honest error inflation, mode L is within ~30% of SALT3's scatter on the same photometry;
the sparse mode remains information-limited, as its error bars now correctly reflect.


#### v0.3.1 note to the SALT3 comparison

With the per-fit error normalization and the boundary gates, the overlap is 28 objects;
CMAGIC rms 0.274 (χ²/dof 3.25 pre-standardization; grey −0.07), residual correlation
with SALT3 r = +0.01 ± 0.20, three >3σ Δμ outliers. Every figure now carries per-point
error bars (x and y), including all 60 per-object panels (failed objects get a named
failure panel), with χ²/dof and the applied error scale printed in each annotation.

### 10.3 Weighted Hubble residuals (v0.3.1)

Inverse-variance-weighted mean residuals against ΛCDM (H₀ = 72, Ω_m = 0.3), errors of the
mean scaled by √(χ²/dof) where above 1. Standardized CMAGIC (M₀ removed by the fit):
all −0.001 ± 0.038 (χ²/dof 0.98, n = 28); mode L −0.000 ± 0.043 (1.00, 16); mode S
−0.003 ± 0.083 (1.03, 12) — the per-mode χ²/dof ≈ 1 is an outcome of the v0.3.1 error
model, not a construction. By redshift bin: +0.036 ± 0.053 (z 0.05–0.15), −0.061 ± 0.071
(0.15–0.25), +0.073 ± 0.140 (0.25–0.35): no redshift trend, and z was never a fitted
covariate.

**The SALT3 redshift drift, diagnosed.** On the same photometry SALT3's weighted
residual with *fixed* fiducial Tripp coefficients (α = 0.14, β = 3.1) drifts by
+0.161 ± 0.049 (3.3σ) from the low to the high bin (+0.180 ± 0.061 on the 28-object
common subset). Most of this is the comparison's own construction, not SALT3, and none
of it is a cosmology-grade statement:

- The subsample's covariates drift with redshift exactly as magnitude-limited selection
  predicts: ⟨x1⟩ = +0.21 → +1.11 and ⟨c⟩ = +0.019 → −0.048 from the low to the high
  tercile (r(x1, z) = +0.43). Any mismatch between the fiducial coefficients and the
  ones this subsample prefers therefore maps covariate drift directly into a
  residual-vs-z slope.
- Refitting (M0, α, β) on the same 28 objects — the identical treatment CMAGIC receives
  in §6.1 — gives α = 0.097 ± 0.021, β = 2.41 ± 0.20 and cuts the drift to
  +0.098 ± 0.050 (2.0σ); the residual-vs-z slope falls from +1.03 ± 0.36 to
  +0.39 ± 0.29 mag per unit z. The coefficient-mismatch × covariate-drift product
  (+0.085 mag) accounts for the removed part.
- The symmetric comparison is then SALT3 +0.098 ± 0.050 vs CMAGIC +0.016 ± 0.103: a
  difference of 0.08 ± 0.11 — no significant method discrimination. CMAGIC's flatness
  currently has ±0.10 mag resolution and cannot claim immunity.
- The residual 2σ drift is what an uncorrected magnitude-limited subsample can show:
  no simulation-based selection/bias corrections (BBC-style) were applied here, as they
  are in every published SALT cosmology analysis, and the fits run on a magnitude-space
  conversion of a cached 60-object subset. Nothing in this table measures cosmology.
- The physics behind the mechanism — why a single fitted β cannot be universal when the
  dust-to-intrinsic color mixture drifts with redshift — is §5.1's account of the two
  extinction philosophies.

What survives as the monitoring target: CMAGIC's redshift stability at its present
±0.10 mag resolution, to be retested on a sample with per-bin errors below 0.05 mag.
Reproduce the decomposition with `examples/sdss_salt3_drift_diagnostic.py`; the
symmetric-treatment Hubble diagram (both methods' standardization fit on the same
28 objects, with the tercile-drift panel) is
`examples/make_hubble_residual_fig_symmetric.py` →
`figures/sdss_hubble_residuals_symmetric.png`. Under that treatment the weighted rms is
0.080 (SALT3 refit) vs 0.196 (CMAGIC) — both are post-fit quantities (3 and 4
parameters respectively fit on the same 28 objects), so the ~2.4× scatter gap is the
honest like-for-like number at v0.3.1. (What that gap measures — and where it goes —
is §10.4.)

### 10.4 DES-SN5YR validation: the scatter gap closes at depth

**Data.** The public DES-SN5YR release (DES Collaboration 2024; light curves Sánchez
et al. 2024; github.com/des-science/DES-SN5YR). The repository caches 122
cosmology-sample SNe Ia at zHEL = 0.05–0.65 — the range where rest-frame B and V are
synthesizable from griz — stratified in redshift (all z < 0.2, 25 per bin above;
`examples/data/des/`), SNANA SMP fluxes converted to magnitudes at SNR ≥ 3, the same
convention as the SDSS-II cache so both runs feed the identical pipeline. The release's
own SALT3 (x1, c) and bias-corrected MU ride along for cross-checks. Precedent for
CMAGIC at these redshifts: Conley et al. (2006) — see Part II's introduction.

**Yield.** Blind, the chain fits 47/122 (21 L, 26 S). Assisted with release metadata
(`t_bmax` = PKMJD; external Δm15 synthesized from the release SALT3 x1, c): 56/122
(23 L, 33 S). SALT3 fits 112/122 of the same tables. Failure census (assisted):
dm15_provenance 30, rms_gate 16, no_bracket 13, window_empty 4, too_few_nights 3.
An architectural finding: the surviving dm15_provenance failures fire *inside the
synthesis stage*, which external metadata cannot reach — passing priors into the
synthesis is a v0.4 item. **At DES depth CMAGIC's real cost is yield, not scatter.**

**Symmetric comparison (n = 53 common objects).** SALT3 Tripp coefficients refit on the
sample (α = +0.056 ± 0.029, β = +2.95 ± 0.29, σ_int = 0.192); CMAGIC standardized per
§6.1 (α_C = +0.083 ± 0.021, β_C = +0.66 ± 0.21, δ_S = −0.026 ± 0.040 — no significant
mode split here; inflation L 1.45, S 1.60; post-fit r(x1) = r(c) = 0). Weighted rms:

| estimator | weighted rms |
|---|---|
| SALT3, refit on sample (same input tables) | 0.196 |
| SALT3, release-grade (bias-corrected MU, full-flux photometry) | 0.214 |
| CMAGIC standardized, all | 0.226 |
| CMAGIC mode L only (n = 23) | **0.197** |
| CMAGIC mode S only (n = 30) | 0.254 |

CMAGIC mode L is at parity with SALT3 on DES — including with the survey's own
bias-corrected distances on the same objects. The overall gap is 1.15×, driven
entirely by mode S. Redshift drift: flat for every treatment (SALT3 fixed
+0.117 ± 0.074 → refit +0.029 ± 0.069; CMAGIC −0.017 ± 0.078; release MU +0.030) —
DES's selection control and bias corrections leave nothing like the SDSS §10.3 slope,
and CMAGIC's z-flatness is confirmed at ±0.08 resolution.

**The resolution of the §10.2 'mystery'** (CMAGIC beats peak methods at low z in
W03/WS06, loses at 2.4× on SDSS): the gap tracks the signal-to-noise contrast between
the light-curve peak and the CMAGIC window, which sits 1.5–2 mag below it. Direct
evidence: the median fit error on B_BV0.6 rises 0.055 → 0.074 → 0.106 → 0.148 mag
across the DES z bins (photon-noise floor), and SALT3's fitted σ_int jumps from 0.064
(SDSS subset) to ≈ 0.19–0.21 (DES subset) — on DES *everyone* is photometry-limited,
and the methods converge. At low z both regions have high SNR and CMAGIC's smaller
intrinsic dispersion wins (W03); at intermediate z on shallow imaging SALT3's
peak-anchored, all-epoch fit wins (SDSS, §10.3); at depth the playing field levels
(DES). The corollary is the independence budget: the CMAGIC–SALT3 residual correlation
is r = +0.01 ± 0.20 on SDSS but r = +0.73 ± 0.07 on DES — shared photon noise
dominates at depth, so CMAGIC's value as an *independent* cross-check is a property of
the intrinsic-limited regime (low z, or deep imaging of bright targets), while at the
survey limit it functions as a consistency check with different systematics, not an
independent one.

Cross-check of our quick SALT3 fits against the release's: Δc median −0.001
(rms 0.068), Δμ rms 0.156 about a +0.31 grey offset (their bias corrections and
zero-point convention). Reproduce: `examples/des_validation.py [--assisted]`,
`examples/des_salt3_comparison.py`, `examples/des_standardize_and_fig.py` →
`figures/des_hubble_residuals_symmetric.png`.
