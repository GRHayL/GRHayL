# M1 neutrino number current and transport

Grey neutrino M1 evolves one number density in addition to the shared E/F
energy moments. The number current is not obtained by treating `N` as a
second copy of the energy flux: current GRHayL derives it from the E/F closure,
comoving moments, and fluid velocity. The implementation is
[ghl_m1_neutrino_number_flux.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_number_flux.c),
with the state and public functions declared in
[ghl_m1.h](../../../GRHayL/include/ghl_m1.h#L1440).

## State and notation

For one species, the undensitized state is

$$
U_\nu=(N,E,F_i),
$$

where `N` is the Eulerian neutrino number density, `E` is the Eulerian energy
density, and `F_i` is the covariant Eulerian energy flux. The shared closure
and comoving calculations provide \(J\), \(H^i\), \(H_i\), \(H_n\), and the
fluid Lorentz factor \(W\). The latter is the matter Lorentz factor; it is not
the number-current normalization below.

## Current normalization and comoving number density

The current construction defines

$$
\Gamma_N=W-\frac{H_n}{J},
\qquad
n_{\rm com}=\frac{N}{\Gamma_N}.
$$

The name `Gamma_N` is intentional: it is a radiation-number current
normalization, not a fluid Lorentz factor, and a valid value may be below one.
For a nonzero number state, the implementation requires finite `J > J_floor`
and `Gamma_N > Gamma_N_floor`. It does not replace `Gamma_N` with `W` on this
path. When `Gamma_N_floor == 0`, the implementation uses the effective strict
floor `64*DBL_EPSILON`; zero initialization therefore does not select a literal
zero threshold. The zero-number branch can still return a zero current when
`Gamma_N` is singular.

## Number transport velocity and spatial current

The contravariant number transport velocity is

$$
V_N^i=\frac{WV^i+H^i/J}{\Gamma_N}.
$$

The undensitized contravariant spatial number current returned by the code is

$$
\mathcal F_N^i=N V_N^i
 =n_{\rm com}\left(WV^i+\frac{H^i}{J}\right).
$$

The implementation validates the velocity with the spatial metric, requiring
it to be causal. It evaluates the velocity algebraically before multiplying by
`N`, so a state with `N == 0` can have a finite nonzero diagnostic velocity and
still has exactly zero number flux. If `N == 0` and `Gamma_N` is singular, the
special zero-density branch returns zero current and zero velocity and uses
`W` only as a harmless stored normalization. This avoids rejecting an
otherwise valid E/F state when no particles are present.

The E/F contractions and the sign of `H_n` used here are defined in
[comoving moments and source projection](m1-comoving-moments-and-source-projection.md).

## Coordinate physical flux and conservative number equation

For a face normal to coordinate direction `d`, the undensitized coordinate
physical number flux is

$$
f_N^{(d)}=\alpha\mathcal F_N^d-\beta^dN.
$$

The conservative number equation therefore has the same lapse/shift structure
as the E/F equations:

$$
\partial_t(\sqrt{\gamma}N)
+\partial_i\!\left[\sqrt{\gamma}
(\alpha\mathcal F_N^i-\beta^iN)\right]
=\alpha\sqrt{\gamma}\,S_N.
$$

`ghl_m1_compute_neutrino_physical_number_flux` returns the undensitized
coordinate flux only. It does not apply a numerical Rusanov flux or a metric
factor. The five-component Rusanov helper uses the same supplied nonnegative
speed for `{N,E,Fx,Fy,Fz}`; the canonical four-point operation then applies
its componentwise blend and performs the one final face densitization. See
[the current flux declarations](../../../GRHayL/include/ghl_m1.h#L1440) and
[the four-point transport contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#transport).

## Floors and admissibility

`N_floor` is independent of the E/F realizability cone. The number-current
path rejects `N < N_floor`; it does not silently repair the input. The explicit
neutrino repair helper can apply `N' = max(N,N_floor)` before deriving a
current, while `J_floor` and `Gamma_N_floor` remain validation controls. The
E/F cone repair is documented in [M1 realizability repair](m1-realizability-repair.md).

## Superseded reduced-current construction

The early implementation whitepaper proposed a reduced number current based on
an Eulerian energy-flux direction, schematically `N*F^i/Fmag`. That construction
and its associated number-current HLL description are historical. Current
GRHayL uses the full E/F-derived expression for \(\Gamma_N\), \(n_{\rm com}\),
and \(V_N^i\) above, together with Rusanov transport. Do not reintroduce a
small-flux direction cutoff or infer `Gamma_N` as `W` for nonzero `N`.

## Focused evidence

- [Neutrino Rusanov/current tests](../../../Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c)
  covers number-current construction, physical number flux, closure reuse, and
  paired stored current/transport cases.
- [Seeded M1 invariants](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  covers pointwise E/F moments, source-side current use, and invariants.
- [Neutrino source-update tests](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
  covers zero-number, endpoint-current, floor, and transactional cases.
- [M1 test map](../../../wiki/gems/radiation-m1/tests-and-fixtures.md)
  records what these tests do and do not establish.
