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
