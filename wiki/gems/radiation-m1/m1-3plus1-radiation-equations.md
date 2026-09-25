# M1 3+1 radiation equations

This leaf records the shared geometric and conservative mathematics used by
the current M1 kernels. It is the E/F part of each grey neutrino state; the
transported neutrino number density `N` has an additional current described in
[number current and transport](m1-number-current-and-transport.md). The
installed header, not a whitepaper, is the API authority:
[ghl_m1.h](../../../GRHayL/include/ghl_m1.h),
[M1_INTEGRATION_CONTRACT.md](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md),
and the active [Radiation build manifest](../../../GRHayL/Radiation/make.code.defn).

## 3+1 geometry and fluid velocity

GRHayL uses the mostly-plus spacetime signature. In lapse-shift form,

$$
ds^2 = -\left(\alpha^2-\beta_k\beta^k\right)dt^2
       +2\beta_i\,dt\,dx^i+\gamma_{ij}\,dx^i dx^j,
$$

where `alpha` is the lapse, `betaU[i]` is the contravariant shift
\(\beta^i\), `gammaDD[i][j]` is \(\gamma_{ij}\), `gammaUU[i][j]` is its
inverse, and `sqrt_detgamma` is \(\sqrt{\gamma}\). The Eulerian normal is

$$
n_\mu=(-\alpha,0,0,0),\qquad
n^\mu=\left(\frac1\alpha,-\frac{\beta^i}{\alpha}\right).
$$

The fluid four-velocity is decomposed as

$$
u^\mu=W(n^\mu+V^\mu),\qquad V^\mu=(0,V^i),
$$

with

$$
V^i=\frac{v^i+\beta^i}{\alpha},\qquad
W=\frac{1}{\sqrt{1-\gamma_{ij}V^iV^j}}.
$$

Here \(v^i=u^i/u^0\) is the coordinate three-velocity stored in the fluid
primitive `vU`. The shared velocity helper computes \(V^i\), lowers it with
\(\gamma_{ij}\), and computes \(W\) from the spatial norm; it rejects a
non-timelike velocity. See
[ghl_m1_compute_eulerian_velocity](../../../GRHayL/Radiation/ghl_m1_utils.h)
and the metric field definitions in
[ghl.h](../../../GRHayL/include/ghl.h).

## Eulerian radiation moments

For a radiation stress-energy tensor \(R^{\mu\nu}\), define

$$
E=R^{\mu\nu}n_\mu n_\nu,
$$

$$
F_\mu=-\gamma_{\mu\alpha}R^{\alpha\beta}n_\beta,
\qquad
P_{\mu\nu}=\gamma_{\mu\alpha}\gamma_{\nu\beta}R^{\alpha\beta}.
$$

The spatial flux is orthogonal to \(n^\mu\). GRHayL stores the covariant
spatial components `F_i`; raising them is an explicit metric operation:

$$
F^i=\gamma^{ij}F_j,
\qquad
P^j{}_i=P^{jk}\gamma_{ki}.
$$

The decomposition reconstructed by the shared stress-energy and closure
kernels is

$$
R^{\mu\nu}=E n^\mu n^\nu+F^\mu n^\nu+n^\mu F^\nu+P^{\mu\nu}.
$$

The current closure validates the spatial trace
\(\gamma_{ij}P^{ij}=E\), symmetry, and positive semidefiniteness. The
realizability cone and its repair are documented separately in
[realizability repair](m1-realizability-repair.md).

## Conservative variables and coordinate fluxes

The local state is undensitized. The conservative quantities in the coordinate
finite-volume equations are

$$
\widetilde E=\sqrt{\gamma}\,E,
\qquad
\widetilde F_i=\sqrt{\gamma}\,F_i,
\qquad
\widetilde N=\sqrt{\gamma}\,N.
$$

For a coordinate direction \(x^d\), the shared E/F physical fluxes are

$$
f_E^{(d)}=\alpha F^d-\beta^dE,
\qquad
f_{F_i}^{(d)}=\alpha P^d{}_i-\beta^dF_i.
$$

The densitized face object is \(\sqrt{\gamma}\,f^{(d)}\). The pointwise
neutrino four-point operation performs this final face densitization exactly
once; callers using prepared volume-weighted operands use its separate
volume-weighted entry point. These are caller-visible boundaries in the
[integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#transport).

## E/F evolution equations

Writing \(S_E,S_i\) for the undensitized radiation interaction source
projections, the conservative E equation is

$$
\partial_t(\sqrt{\gamma}E)
+\partial_j\!\left[\sqrt{\gamma}(\alpha F^j-\beta^jE)\right]
=\sqrt{\gamma}\left[
\alpha P^{ij}K_{ij}-F^j\partial_j\alpha+\alpha S_E\right].
$$

The covariant momentum equation is

$$
\partial_t(\sqrt{\gamma}F_i)
+\partial_j\!\left[\sqrt{\gamma}(\alpha P^j{}_i-\beta^jF_i)\right]
=\sqrt{\gamma}\left[
-E\partial_i\alpha+F_j\partial_i\beta^j
+\frac{\alpha}{2}P^{jk}\partial_i\gamma_{jk}+\alpha S_i\right].
$$

The energy equation has two geometry terms and the momentum equation has
three. Current GRHayL computes them as densitized quantities in
[ghl_m1_sources_geometry.c](../../../GRHayL/Radiation/ghl_m1_sources_geometry.c):

$$
\widetilde S_E^{\rm geom}=\sqrt{\gamma}
\left(\alpha P^{ij}K_{ij}-F^j\partial_j\alpha\right),
$$

$$
\widetilde S_i^{\rm geom}=\sqrt{\gamma}
\left(-E\partial_i\alpha+F_j\partial_i\beta^j
+\frac{\alpha}{2}P^{jk}\partial_i\gamma_{jk}\right).
$$

Interaction sources are added separately with \(\alpha\sqrt{\gamma}\), as
shown by [ghl_m1_compute_matter_coupling_sources](../../../GRHayL/Radiation/ghl_m1_matter_coupling_sources.c)
and the neutrino source equations in the companion leaf.

## Number equation

Let \(\mathcal F_N^i\) denote the undensitized contravariant spatial number
current derived from the E/F moments. Its coordinate physical flux is

$$
f_N^{(d)}=\alpha\mathcal F_N^d-\beta^dN.
$$

The corresponding conservative equation is

$$
\partial_t(\sqrt{\gamma}N)
+\partial_i\!\left[\sqrt{\gamma}
(\alpha\mathcal F_N^i-\beta^iN)\right]
=\alpha\sqrt{\gamma}\,S_N.
$$

The current number source and the distinction between \(\mathcal F_N^i\) and
the coordinate flux are given in
[number current and transport](m1-number-current-and-transport.md). The
library supplies pointwise local quantities; the host owns reconstruction,
finite-volume divergence, time integration, boundaries, and mesh traversal.

## Photon-whitepaper boundary

The 3+1 projections, E/F/P decomposition, conservative equations, and metric
flux factors above are shared M1 mathematics and remain applicable to neutrino
transport. Photon-only material from the equations whitepaper must not be
copied into the neutrino rule: in particular, the photon equilibrium choice
\(J_{\rm eq}=a_RT^4\), photon opacity prescriptions, and photon microphysics
are not supplied by the current Radiation kernels. Neutrino `n_eq`, `J_eq`,
emissivities, and opacities arrive in a frozen provider bundle; see
[ghl_m1.h rate fields](../../../GRHayL/include/ghl_m1.h).

The older methods whitepaper's HLL flux is also not implied by these continuum
equations. The current canonical neutrino face method is four-point blended
Rusanov; generic HLL helpers do not define that path. The canonical face
operation is listed in the production-surface table of
[TRACEABILITY.md](../../../GRHayL/Radiation/TRACEABILITY.md#production-boundary).

## Focused evidence

- [Seeded M1 invariants](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  exercises closure, moments, stress-energy, geometry sources, and exchange.
- [M1 error handling](../../../Unit_Tests/unit_test_m1_error_handling.c)
  exercises invalid geometry/state and unchanged-output behavior.
- [M1 test and fixture map](../../../wiki/gems/radiation-m1/tests-and-fixtures.md)
  defines the evidence boundary; these unit tests do not prove host mesh or
  time-integrator behavior.
