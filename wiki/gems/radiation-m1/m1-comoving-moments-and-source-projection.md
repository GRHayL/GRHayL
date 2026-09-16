# M1 comoving moments and source projection

The shared M1 closure produces an Eulerian pressure tensor `P^{ij}`. Source
terms and the neutrino number current then require the radiation moments in the
fluid frame. This leaf gives the current contractions and signs; it does not
restate the complete public API contract. The implementation is
[ghl_m1_comoving_moments.c](../../../GRHayL/Radiation/ghl_m1_comoving_moments.c),
with declarations in [ghl_m1.h](../../../GRHayL/include/ghl_m1.h#L106).

## Fluid-frame projector and moments

For the fluid four-velocity \(u^\mu\), define

$$
h^\mu{}_{\nu}=\delta^\mu{}_{\nu}+u^\mu u_\nu,
\qquad
h^\mu{}_{\nu}u^\nu=0.
$$

The comoving radiation moments are

$$
J=R^{\mu\nu}u_\mu u_\nu,
$$

$$
H^\mu=-h^\mu{}_{\alpha}R^{\alpha\beta}u_\beta,
\qquad H^\mu u_\mu=0.
$$

`J` is the fluid-frame radiation energy density. `HU[i]` and `HD[i]` in the
current `ghl_m1_comoving` record are the contravariant and covariant spatial
components of \(H^\mu\); `Hn` is its projection onto the Eulerian normal.

## E/F/P to J/H contractions

Use the E/F state with the current covariant-flux convention:

$$
F^i=\gamma^{ij}F_j,
\qquad
P_{ij}=\gamma_{ik}\gamma_{j\ell}P^{k\ell},
\qquad
F\!\cdot\!V=F_iV^i.
$$

The code evaluates

$$
J=W^2\left(E-2F_iV^i+P_{ij}V^iV^j\right),
$$

and

$$
H^i=W\left(F^i-P^{ij}V_j-JV^i\right).
$$

Lowering gives \(H_i=\gamma_{ij}H^j\). The normal projection used by the
source formulas is

$$
H_n\equiv H^\mu n_\mu=-V_iH^i.
$$

That last minus sign is important: `Hn` is not \(+V_iH^i\). The implementation
forms it explicitly after lowering `HU` with `gammaDD`; see
[ghl_m1_comoving_moments.c](../../../GRHayL/Radiation/ghl_m1_comoving_moments.c#L43).

## Radiation-side four-force projection

For a radiation-side four-force written in the shared form

$$
G^\mu=Q\,u^\mu-\kappa_{\rm tr}H^\mu,
$$

the conservative source projections are

$$
S_E=-G^\mu n_\mu=QW+\kappa_{\rm tr}H_n,
$$

$$
S_i=G^\mu\gamma_{i\mu}=QWV_i-\kappa_{\rm tr}H_i.
$$

For neutrinos, \(Q\) and \(\kappa_{\rm tr}\) are supplied by the grey rate
bundle and are defined in
[neutrino source equations](m1-neutrino-source-equations.md). The current
implementation evaluates these projections in
[ghl_m1_neutrino_sources.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c#L25).

The sign convention is that \(G^\mu\) acts on radiation. The matter source is
the equal-and-opposite conservative contribution, not another application of
the radiation source; see [exchange equations](m1-exchange-equations.md).

## Geometry sources are separate

The source projection above is an interaction source. Metric derivatives and
extrinsic curvature generate separate geometry sources in the conservative
E/F equations. Do not fold those terms into `ghl_m1_sources` and then pass the
result to the matter-coupling helper as though all terms represented
radiation-matter exchange. The separation is visible in
[ghl_m1_sources_geometry.c](../../../GRHayL/Radiation/ghl_m1_sources_geometry.c#L52)
and [ghl_m1_matter_coupling_sources.c](../../../GRHayL/Radiation/ghl_m1_matter_coupling_sources.c#L4).

## Validation and numerical meaning

Before these contractions, the current kernels validate the metric, the E/F
realizability state, and the closure tensor. A nonfinite or negative `J` is
rejected. The outputs are assembled in a local candidate and published only
after all components are finite. This makes the moments suitable for the
transactional source and current paths, but does not make an inadmissible
input state valid; repair belongs at the documented repair boundary.

The contraction is algebraic and pointwise. It does not perform EOS lookup,
rate evaluation, a grid loop, or time integration. Fluid primitives are frozen
for the local neutrino source/update call.

## Photon-whitepaper boundary

The projector, moment contractions, and source projection signs are shared
mathematics. The photon equations whitepaper's special substitution

$$
Q=\chi_{\rm abs}(J_{\rm eq}-J),
\qquad J_{\rm eq}=a_RT^4,
$$

is a photon/LTE model and is not a neutrino rule. Current neutrino code instead
uses \(Q=\eta_E-\kappa_{a,E}J\), with weak-equilibrium targets and frozen
provider coefficients. Likewise, an old HLL discretization in the methods
whitepaper is not part of this moment projection.

## Focused evidence

- [Seeded neutrino M1 invariants](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  covers closure, comoving moments, stress-energy, and source projections.
- [Neutrino source-update tests](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
  covers source-path endpoint/current behavior and transactional failures.
- [M1 error handling](../../../Unit_Tests/unit_test_m1_error_handling.c)
  covers invalid moment/closure inputs and unchanged outputs.

The complete boundary remains in the
[M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#fixed-m1-method).
