# Instantaneous frozen-rate source fixture

`m1_thcm1_instantaneous_sources.m1` is the retained, offline fixture for the
instantaneous frozen-rate interaction-source operation. It contains all 56
retained baseline/perturbed pairs for all three species: 168 records total.
The source campaign provenance is maintained in the separate external campaign
workspace. Its receipt, sidecar, source snapshot, producer command, and raw
rows are not copied into the public GRHayL repository. Normal tests do not run
THC_M1 or require `THCM1_ROOT`.

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
recomputes normalization from both input vectors and validates the exact
producer initializer, coordinate-velocity construction, finite `u0`, source
statuses, IDs, and paired response semantics.
Normalization products are rounded separately before addition, matching the
exporter's binary64 arithmetic even when the test compiler contracts other
expressions into fused multiply-add instructions. This affects test metadata
only, not production arithmetic or comparison tolerances.

## Retained result

The seeded invariant test preserves all 168 records:

- 162 records pass current-GRHayL/retained-source agreement checks.
- Six records are explicit local policy checks, not agreement claims. They are
  the three species for each of:
  `rngpkt-v2-rd-radiation-anchor-energy-floor-a01` and
  `rngpkt-v2-rd-metric-anchor-offdiagonal-spd-a01`.
- Those six checks cover the known zero-Eulerian-flux, moving-fluid cases where
  current closure admissibility fallback and THC invalid trace-pressure
  handling select different `Gamma_N` behavior. Any new or changed
  classification fails.

External promotion must fail closed for an incomplete receipt, missing or
empty selected source list, duplicate or unsafe selected paths, missing
artifacts, sidecar/result/command binding mismatch, incomplete role/species
coverage, duplicate IDs, packet/state/volume mismatch, unpublished or
unsuccessful producer outputs, invalid `Gamma_N`, and any record-count or
input-schema mismatch. It exports every retained role/species record; it does
not select only passing rows. The public package receives only the resulting
portable payload and compact admission metadata.

## Remaining issue

The six local classifications remain intentional compatibility checks. They do
not establish THC/current closure equivalence for those cases. No production
change, tolerance relaxation, or replacement of retained values is part of this
fixture.
