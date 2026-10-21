# Instantaneous frozen-rate source fixture

`m1_thcm1_instantaneous_sources.bin` is the retained, offline fixture for the
instantaneous frozen-rate interaction-source operation. It contains all 56
retained baseline/perturbed pairs for all three species: 168 records total, of
which 48 have identical baseline and perturbed inputs. The records come from the
fresh instantaneous-source THC_M1 campaign; its receipt, sidecar, source
snapshot, producer command, and raw rows are held in the producer evidence
outside the GRHayL checkout. Normal tests do not run THC_M1 or require
`THCM1_ROOT`.

## Operation contract

The fixture policy is `source_a1_a2_rate_normalized_v1`. Each record has 40
consumed input scalars:

1. state `N,E,F0,F1,F2` (5);
2. lapse, shift, and row-major `gammaDD` (13);
3. GRHayL coordinate velocity `vU[3]` (3);
4. scalar frozen rates and species controls (12);
5. the three number and three energy pair-emissivity slots (6), explicitly
   zero because the retained producer initializes them to zero; and
6. `sqrt_detgamma` (1), consumed by matter-coupling densitization.

Baryon density, `dt`, and `source_therm_limit` are endpoint/host controls and
are intentionally excluded from this instantaneous operation. Changes only in
those fields therefore produce identical consumed inputs and zero sensitivity.

The nine outputs are number source, energy source, three momentum sources,
matter energy coupling, and three matter momentum couplings. The evaluator
recomputes normalization from both input vectors and validates the recorded
initializer and coordinate-velocity conventions, finite `u0`, source statuses,
IDs, and paired response semantics.
Normalization products are rounded separately before addition, matching the
exporter's binary64 arithmetic even when the test compiler contracts other
expressions into fused multiply-add instructions. This affects test metadata
only, not production arithmetic or comparison tolerances.

## Retained result

The seeded invariant test preserves all 168 records and carries two kinds of
named fallback-policy classification, both as local checks rather than
cross-code agreement claims:

- 159 records retain strict paired agreement comparisons.
- Six records are explicit exact-zero-flux local policy checks: the three
  species for each of
  `rngpkt-v2-rd-radiation-anchor-energy-floor-a01` and
  `rngpkt-v2-rd-metric-anchor-offdiagonal-spd-a01`. They cover the known
  zero-Eulerian-flux, moving-fluid cases where current closure admissibility
  fallback and THC invalid trace-pressure handling select different
  `Gamma_N` behavior.
- Three records are explicit covariance-anchor fallback-policy checks: the
  three species of the
  `rngpkt-v2-rd-radiation-anchor-realizability-cone-a01:baseline:<species>`
  pairs derived from the retained pointwise anchor record. The retained source
  datum carries the closure input and the frozen rates but no serialized
  reference pressure; the named realization difference is owned by the
  pointwise fallback-policy classification. At both paired endpoints these records
  require the current identity-metric, x-aligned-flux flagged fallback with
  zero root diagnostics, the analytic Eulerian Minerbo fallback pressure and
  closure factor with valid trace, tensor symmetry, and diagonal covariance
  positivity, finite actual sources, and the established exact
  radiation/matter exchange check. They compare the named operation under its
  published closure fallback rather than the nonrealizable retained closure.

Any new or changed classification fails. These nine records are checked locally
rather than compared numerically: their stored THC outputs are not part of any
replay comparison, so a small change to them need not fail the test.

Any replacement exporter must fail closed for an incomplete
receipt, missing or empty selected source list, duplicate or unsafe selected
paths, missing artifacts, sidecar/result/command binding mismatch, incomplete
role/species coverage, duplicate IDs, packet/state/volume mismatch,
unpublished or unsuccessful producer outputs, invalid `Gamma_N`, and any
record-count or input-schema mismatch. It must export every retained
role/species record rather than selecting only passing rows. A replacement
TestData package should contain the portable payload and compact package
metadata, and takes effect only when `.github/radiation-testdata-ref` is
re-pinned to the revision that publishes it.

## Remaining issue

The named local classifications remain intentional compatibility checks. They
do not establish THC/current closure equivalence for those cases. No production
change, tolerance relaxation, or replacement of retained values is part of this
fixture.
