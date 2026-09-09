# Radiation M1 Compatibility Evidence Boundary

This page limits claims to the library operations present in the checkout.
The neutrino M1 implementation has one primary numerical path, with
a flagged Eulerian Minerbo admissibility fallback for finite non-PSD or
invalid exact-zero-flux closure candidates. The pointwise operations and source-update contract below describe
that path; they do not establish a downstream grid evolution.

## What is testable here

- component order `{N,E,Fx,Fy,Fz}`;
- four-point limiter, sawtooth detection, opacity suppression, and high/low
  blending;
- exactly-once face densitization and the canonical no-cap/no-separate-
  diffusion transport boundary;
- frozen-rate source branches and transactional no-update behavior;
- separate total-number and charged-current lepton exchange.

The corresponding implementations are listed in
[`TRACEABILITY.md`](../../../GRHayL/Radiation/TRACEABILITY.md).
The focused tests in `scripts/test_radiation.py` cover the zero-flux tensor
invariants, coupled pair source conservation, and provider channel mapping.
The list above is broader than that regression coverage; it does not imply
that every transport or source branch is tested.

## What is not claimed

No external-source comparison, downstream framework run, grid evolution,
schedule/AMR result, complete evolution equivalence, or physical-validation
result is established by the listed library operations. A downstream consumer owns
any such campaign and its reporting.

## Canonical-method boundary

The host must use the four-point transport operation for neutrino M1 faces and
provide uncapped metric light-cone speeds. The operation does not provide a
separate diffusion correction. The existing finite-difference/Newton source
solver remains the local implicit solver. Legacy policy fields may remain for
source or ABI compatibility, but cannot select a different numerical method.
The host must apply one limiter scalar to every coupled species, matter, and
`Y_e` increment; no host implementation is included here.
