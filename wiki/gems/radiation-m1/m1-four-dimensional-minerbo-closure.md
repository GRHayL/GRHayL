# M1 four-dimensional Minerbo closure

The current production closure maps an admissible undensitized E/F state and
the local metric/fluid velocity to a contravariant spatial pressure tensor
`P^{ij}`. Its primary construction is four-dimensional and covariant; it is
not an analytic Eulerian-only formula or an HLL recipe. The implementation
authority is
[ghl_m1_closure.c](../../../GRHayL/Radiation/ghl_m1_closure.c), with public
fields and statuses in [ghl_m1.h](../../../GRHayL/include/ghl_m1.h).

## Spacetime data used by the closure

The closure builds the ADM spacetime metric from the face/cell lapse, shift,
and spatial metric:

$$
g_{\mu\nu}=
\begin{pmatrix}
-\alpha^2+\beta_k\beta^k & \beta_j\\
\beta_i & \gamma_{ij}
\end{pmatrix},
\qquad
g^{\mu\nu}=
\begin{pmatrix}
-\alpha^{-2} & \beta^j/\alpha^2\\
\beta^i/\alpha^2 & \gamma^{ij}-\beta^i\beta^j/\alpha^2
\end{pmatrix}.
$$

The Eulerian normal is \(n_\mu=(-\alpha,0,0,0)\). For the covariant Eulerian
flux stored by GRHayL, orthogonality fixes the time component used internally:

$$
F_0=\beta^iF_i,
\qquad F_\mu=(F_0,F_i).
$$

The fluid velocity is represented as

$$
u^0=\frac{W}{\alpha},
\qquad
u^i=W\left(V^i-\frac{\beta^i}{\alpha}\right),
$$

where \(V^i=(v^i+\beta^i)/\alpha\). The spatial metric and its inverse must
be coherent; current validation requires symmetric positive-definite,
inverse-consistent `gammaDD`/`gammaUU` and a consistent positive determinant.

## Thin and thick four-dimensional tensors

For a nonzero four-flux, the current code forms a normalized flux shape
\(\widehat F_\mu=F_\mu/F_{\rm scale}\), where `F_scale` is chosen from the
largest absolute component. It then forms

$$
F_{\rm shape}^2=g^{\mu\nu}\widehat F_\mu\widehat F_\nu,
\qquad
P^{\rm thin}_{\mu\nu}=E\,
\frac{\widehat F_\mu\widehat F_\nu}{F_{\rm shape}^2}.
$$

The normalization is an implementation detail with mathematical purpose: it
avoids overflow or underflow in the direction dyad. An exactly zero Eulerian
flux produces a zero thin candidate; it is not treated by an arbitrary
small-flux cutoff.

The thick candidate is the relativistic isotropic/diffusion-limit tensor. In
the fluid frame its isotropic pressure part is

$$
K^{\mu\nu}_{\rm thick}=\frac{J}{3}
\left(g^{\mu\nu}+u^\mu u^\nu\right),
$$

and the corresponding stress tensor has the decomposition

$$
R^{\mu\nu}=J u^\mu u^\nu+H^\mu u^\nu+u^\mu H^\nu+K^{\mu\nu}_{\rm thick}.
$$

The closure implementation constructs an algebraically equivalent covariant
thick tensor directly from E/F and the fluid velocity, then projects it to
spatial `P^{ij}`. Its thick-limit comoving-energy estimate is also exposed by
the shared helper:

$$
J_{\rm thick}=\frac{3}{2W^2+1}
\left[(2W^2-1)E-2W^2F_iV^i\right].
$$

`J_thick` is a diffusion-limit diagnostic/helper, not the solved closure
moment `J` at an arbitrary \(\xi\). See
[ghl_m1_Jthick.c](../../../GRHayL/Radiation/ghl_m1_Jthick.c).

For numerical stability, the closure workspace divides the thick equations
by \(W^2\) before evaluating their differences. With
\(A=E-F_iV^i\) and \(q=1/W^2\), the equivalent expressions are

$$
\frac{J_{\rm thick}}{3}=\frac{2A-qE}{2+q},
\qquad
\widetilde H_\mu=\frac{F_\mu}{W}
-Wv_\mu\left(A+\frac{J_{\rm thick}}{3}\right).
$$

These avoid subtracting intermediate terms of order \(W^2E\). Large-energy
workspace construction uses the same equations in units of E.
Here \(\widetilde H_\mu\) denotes the Eulerian-spatial projection of the
thick-limit comoving flux, not its complete spacetime covector.

## Minerbo interpolation

The current Eddington factor is the polynomial

$$
\chi(\xi)=\frac13+
\frac{\xi^2(6-2\xi+6\xi^2)}{15},
\qquad 0\le\xi\le1.
$$

Define

$$
d_{\rm thin}=\frac{3\chi-1}{2},
\qquad
d_{\rm thick}=\frac{3(1-\chi)}{2}.
$$

The four-dimensional pressure candidate is

$$
P_{\mu\nu}(\xi)=d_{\rm thin}P^{\rm thin}_{\mu\nu}
+d_{\rm thick}P^{\rm thick}_{\mu\nu}.
$$

The code inserts this pressure into the Eulerian stress decomposition,

$$
R_{\mu\nu}(\xi)=E n_\mu n_\nu+F_\mu n_\nu+n_\mu F_\nu+P_{\mu\nu}(\xi),
$$

and computes the fluid-frame moments implied by that candidate.

## Scalar consistency root

For each trial \(\xi\), compute

$$
J(\xi)=R^{\mu\nu}(\xi)u_\mu u_\nu,
$$

$$
H_\mu(\xi)=-\left(\delta^\alpha{}_{\mu}+u^\alpha u_\mu\right)
R_{\alpha\beta}(\xi)u^\beta,
\qquad
H^2(\xi)=g^{\mu\nu}H_\mu H_\nu.
$$

The closure condition is

$$
g(\xi)=J(\xi)^2\xi^2-H^2(\xi)=0.
$$

The published diagnostic stores the normalized residual

$$
r_{\rm root}=\frac{|g(\xi)|}{s_{\rm residual}},
$$

where the implementation chooses an overflow-safe scale from the energy,
\(J^2\xi^2\), and \(|H^2|\). When `E` is large enough that `E^2` is not
representable, the code evaluates the homogeneous relation after scaling by
`E`. This preserves a meaningful residual test rather than allowing an
intermediate square to overflow.

## Root and publication behavior

The production path is a bracketed, safeguarded Brent-style scalar solve:

1. Validate the metric and E/F realizability, construct the workspace, and
   evaluate the residual at \(\xi=0\) and \(\xi=1\).
2. If an endpoint is zero or within the implementation's internal roundoff
   gate, select the corresponding endpoint and report convergence.
3. If the endpoints do not bracket a sign change, choose the endpoint with the
   smaller normalized residual and report `endpoint_fallback`.
4. Otherwise iterate inside the bracket up to
   `closure_root_max_iterations`, using the configured interval tolerance.
5. Evaluate the candidate, reject a normalized residual above
   `closure_root_residual_tolerance`, then validate its finite values, spatial
   trace, symmetry, and positive semidefiniteness before publication.

At \(\xi=0\), the evaluator computes the thick moments directly from E/F,
using the spatial projection of \(\widetilde H_\mu\) above and subtracting
the temporal contribution in \(H^2\). It does not boost the rounded pressure
entries back to obtain this endpoint residual: that operation amplifies
pressure roundoff by \(W^2\) and can change endpoint selection. The endpoint
calculation is energy-normalized internally, but its signed residual uses
the same units as the other trial points in the Brent solve. The endpoint
gate and preference for the zero endpoint are unchanged. This stabilizes
the existing selection rule; it does not assume the scalar root is unique.

The public `ghl_m1_closure` record exposes `chi`, the physical reduced flux
`xi`, root residual, iteration count, solve status, and
`four_point_compatibility`. `four_point_compatibility` is diagnostic: it is
true for the primary full-four-dimensional construction and false for the
exceptional admissibility fallback. It is not a runtime method selector.

## Exceptional admissibility fallback

If a finite primary candidate fails the pressure positive-semidefinite check,
the current code constructs an Eulerian Minerbo pressure tensor from the
metric-normalized Eulerian flux factor, validates it, and publishes it with
`solve_status = endpoint_fallback` and
`four_point_compatibility = false`. An exactly zero Eulerian flux can use this
same fallback when the moving-fluid primary candidate cannot satisfy the
trace validation. A nonfinite candidate, invalid root data, or failed fallback
is an error; the code does not publish a partially validated tensor.

This fallback is deliberately observable. The compatibility source-update
branches are fail-closed on this fallback unless
`allow_closure_fallback` explicitly permits the candidate. The default
`ghl_m1_neutrino_source_grhayl_implicit` path records the fallback diagnostic
and can publish a successful solve; it does not enforce
`allow_closure_fallback`. The relevant validation and fallback code is in
[ghl_m1_closure.c](../../../GRHayL/Radiation/ghl_m1_closure.c),
[ghl_m1_utils.h](../../../GRHayL/Radiation/ghl_m1_utils.h), and the
[neutrino source-update dispatcher](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c).

## Closure-method boundary

An analytic radiation-rest-frame boost relation or its Levermore-form rewrite
is not the current production algorithm. Current GRHayL solves the full
four-dimensional Minerbo consistency root described above. An HLL
signal-speed construction is not a closure rule and does not define current
neutrino transport.

## Focused evidence

- [Seeded M1 invariants](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  exercises closure, moments, and tensor-derived quantities.
- [M1 error handling](../../../Unit_Tests/unit_test_m1_error_handling.c)
  exercises invalid closure/metric/state boundaries and unchanged outputs.
- [M1 contract](../../../wiki/gems/radiation-m1/neutrino-m1-contract.md#fixed-numerical-path)
  records the primary path and fallback policy at the KB boundary.
