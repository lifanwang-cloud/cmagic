# Known issues and contract notes (v0.1.0)

Deviations and limits relative to README.md / docs/METHOD.md, recorded rather than
silently patched:

1. **WS06 K-correction worked example is not fully reproducible in closed form.**
   The Table 9 (SN 1992al) output-level corrections include light-curve-refit effects
   beyond the Table 10 polynomials. This package reproduces the distance-relevant
   K(B_BV) to 1.5 mmag in magnitude; the B_max/V_max entries differ by up to
   0.03 mag and the published table's sign convention (corrected = raw + K) is
   opposite to the application validated end-to-end here (corrected = raw - K).
   At z <= 0.03 the ambiguity is <= 0.01 mag in distance modulus. Tests enforce
   the documented bounds, not the unreachable 1-mmag ideal.
2. **`photometry_defect` gate**: the automatic detector currently covers only
   aggregator duplicates (OSC rows repeating a VizieR catalog, dropped by bibcode).
   Duplicated-column or copied-row defects inside a single catalog (seen in
   archival compilations) are NOT auto-detected; inspect the diagnostic panel.
3. **`system="csp"` is bookkeeping only.** Magnitudes are used as supplied; no
   CSP-to-Johnson transformation is applied. CSP natural B, V are close enough to
   Bessell for the CMAGIC branch at the current calibration dispersion, but the
   provenance field records the declared system for the user's audit.
4. **`beta_gate` is diagnostic-level, never fatal.** The calibrated chain uses the
   fixed population slope (1.94), so an invalid free slope marks the gate 'fail'
   and withholds `beta_free_valid`, without failing the distance. (METHOD section 7
   lists it among the gates; this is the implemented reading.)
5. **`dm15_provenance` cannot see external fit quality.** The gate fires when the
   decline rate is neither supplied nor measurable from the data. If the caller
   supplies a dm15 from a bound-hit light-curve fit, the window will inherit it
   (METHOD section 2 caveat); provenance of a supplied dm15 is the caller's
   responsibility.
6. **Slow decliners** (dm15 <~ 0.9): the 7/dm15 entry can admit the tail of the
   post-maximum shoulder; the fixed-slope distance is robust at the few-percent
   level, the free slope is biased high (documented in METHOD section 2).
7. **WS06 bilinear standardization** is anchored to the W03 Table 3 zero point at
   dm15 = 1.1 because WS06 tabulate no absolute zero point (METHOD section 6 notes
   the anchoring).

## v0.2.0 (high-redshift extension)

8. **Grey Hubble zero point.** The SDSS-II validation shows a −0.17 mag grey offset of
   the synthesized-chain distance moduli against ΛCDM(H₀=72): the W03 M_BV zero point
   (Vega, low-z, H₀=65-rescaled) does not transfer exactly through template synthesis
   and the SNANA/AB calibration of the validation photometry. It is a single global
   constant (no measured z-tilt beyond the errors), documented rather than
   recalibrated away; users comparing absolute distances across the low-z and high-z
   paths should calibrate the offset on overlap objects.
9. **Error budgets underestimate the observed high-z scatter** (validation χ²/dof ≈ 4
   for mode L, ≈ 9 for mode S). Known missing terms: SMP photometry systematics,
   host-extinction estimation on synthesized peak magnitudes, and the sparse-mode
   extinction fragility at high color leverage. The Cookbook §10 passing criterion is
   therefore not met at v0.2.0; numbers are reported as measured.
10. **Monte Carlo scope.** The synthesis MC perturbs the photometry and re-derives the
    per-epoch warps at fixed (t₀, stretch) [hsiao] or fixed converged (x₁, c)
    [salt3-nir]; template-parameter uncertainty enters through the §8.1 insensitivity
    spread (`K_systematic`), not the MC covariance.
11. **MW extinction at high z** is removed from the observed bands to first order
    (band-integrated CCM at the observed wavelengths) before synthesis; SDSS-field
    E(B−V) ≲ 0.1 makes higher-order terms negligible there.
12. **`system=` remains bookkeeping**; the high-z path assumes AB for survey filters
    and Vega for bessell/csp/2mass names.
