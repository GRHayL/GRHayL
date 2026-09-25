# M1 realizability repair

The shared E/F repair keeps a radiation state inside the metric-dependent
realizability cone before closure, transport, or source evaluation. The
canonical implementation is
[ghl_m1_realizability_repair.c](../../../GRHayL/Radiation/ghl_m1_realizability_repair.c),
and the only supported policy is declared in
[ghl_m1.h](../../../GRHayL/include/ghl_m1.h).

## E/F cone

For the covariant spatial flux stored as `F_i`, define the metric norm and
reduced flux factor

$$
|F|_\gamma=\sqrt{\gamma^{ij}F_iF_j},
\qquad
f=\frac{|F|_\gamma}{E}.
$$

The configured margin is `epsilon_c`, with

$$
c=1-\epsilon_c,
\qquad
f^2\le c.
$$

Thus the permitted norm ratio is \(\sqrt{1-\epsilon_c}\), not
\(1-\epsilon_c\). The historical parameter member
`one_minus_epsilon_c_sq` stores the squared bound `1 - epsilon_c`; its name
does not change this interpretation. The metric is required to be symmetric
positive-definite and inverse-consistent before the norm is evaluated.

## Canonical repair map

First apply only the configured energy floor:

$$
E' = \max(E,E_{\rm floor}).
$$

Compute \(f'=|F|_\gamma/E'\). If \(f'\le\sqrt c\), retain the flux. If the
state violates the cone, rescale the covector by

$$
a=\frac{c}{(f')^2},
\qquad
F'_i=aF_i.
$$

The repaired reduced norm is

$$
\frac{|F'|_\gamma}{E'}=a f'=\frac{c}{f'}\le\sqrt c,
$$

so the squared inequality is satisfied. The factor is the squared-ratio
rescale `E^2(1-epsilon_c)/F^2` expressed using the post-floor energy and
metric norm. It is not the linear norm-ratio rescale used by some older
whitepaper descriptions.

The operation is component-preserving in direction and changes only the
energy floor and, when needed, the flux magnitude. It computes the original
and permitted norms and applied scale internally while validating the local
candidate. The standalone `ghl_m1_realizability_repair` API publishes only
the repaired E/F state and success/error status; it does not expose per-call
norm, scale, branch, or repaired-state diagnostics.

## Neutrino state repair

The neutrino state adds a scalar number density:

$$
U_\nu=(N,E,F_i).
$$

`ghl_m1_repair_neutrino_state` applies the independent number floor

$$
N'=\max(N,N_{\rm floor}),
$$

then delegates E/F repair to the shared map. `N_floor` is not part of the
E/F cone, and `J_floor` and `Gamma_N_floor` are validation floors for
comoving-moment/current evaluation rather than automatic additions to the
state. When `Gamma_N_floor == 0`, current number-current evaluation uses the
effective strict floor `64*DBL_EPSILON`. The endpoint backward-Euler number
formula does not inject a floor;
source-update code repairs or rejects at its explicit repair boundary.

No finite-flux cutoff is used to avoid the closure's exact-zero branch. A
zero-flux state is handled by the current closure implementation and its
observable admissibility fallback policy.

## Transaction and accounting rules

The shared repair computes a local candidate and writes it only after all
validation succeeds. A nonfinite state, invalid metric, invalid parameter
bundle, or failed norm check leaves the state unchanged. The neutrino repair
also increments caller-owned counters and absolute component-wise repair
budgets when diagnostics are supplied; those budgets are not signed physical
exchange deltas.

The local source solvers use the same boundary transactionally. After required
output pointers are valid, a source failure leaves the output state at its
source base and returns a zero exchange packet. An exhausted implicit retry
schedule is reported as terminal no-update rather than publishing a Newton
guess. The paired electron-flavor update applies this rule to both species
atomically.

## Where repair belongs

Repair is required before a state is passed to closure or a physical flux is
formed. It is not silently embedded in the neutrino implicit residual: trial
admissibility is checked by the nonlinear solve, while the final state repair
occurs at the documented endpoint boundary. A host must also account for any
repair-induced state mutation if it is enforcing coupled matter/lepton
conservation.

## Photon and historical-method boundary

The realizability cone and canonical rescale are shared M1 mathematics. Photon
equilibrium formulas, photon opacity models, and the old photon HLL/diffusion
recipe do not alter this neutrino repair map. In particular, do not replace
the current squared-ratio rule with the obsolete reduced-number-current or
HLL-era whitepaper behavior.

## Focused evidence

- [M1 error handling](../../../Unit_Tests/unit_test_m1_error_handling.c)
  covers energy/flux repair, invalid states, metric checks, and transactional
  output behavior.
- [Neutrino source-update tests](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
  covers number-floor/E/F endpoint repair and failure publication.
- [Seeded M1 invariants](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  covers repaired admissible states used by closure, moments, and sources.
- [Neutrino M1 contract](../../../wiki/gems/radiation-m1/neutrino-m1-contract.md#fixed-numerical-path)
  records the repair and fallback boundary.
