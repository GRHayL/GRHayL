# Local implicit Newton solve and failure policy

The neutrino source update solves the stiff energy/flux interaction locally.
It is a pointwise E/F solve with a post-solve number update; it is not a grid
solver and it does not refresh matter or microphysics while iterating. The
source implementation is split across
[`ghl_m1_neutrino_implicit_residual.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit_residual.c),
[`ghl_m1_neutrino_implicit_jacobian.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit_jacobian.c),
[`ghl_m1_neutrino_implicit_solve.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit_solve.c),
and the shared driver
[`ghl_m1_newton.c`](../../../GRHayL/Radiation/ghl_m1_newton.c).

## Frozen local problem

The solver receives a frozen metric, frozen fluid primitives, a validated
frozen rate bundle, timestep, baryon normalization, and one undensitized
neutrino state. The caller owns EOS and `Con2Prim`; the solver does not call
either. Rates, primitive variables, and opacity are not refreshed during a
Newton iteration or fallback substep.

The timestep is coordinate time. Every local source calculation uses

```text
dt_alpha = alpha * dt.
```

The source base is the state passed to the homogeneous solver. In the public
source dispatcher this is `state_transport`, the state after the host's
transport stage; `state_input` is retained for pre-transport validation and is
not the base for source increments.

## Residual unknowns

Newton solves the four densitized E/F unknowns

```text
U = (tilde_E, tilde_F_x, tilde_F_y, tilde_F_z)
  = sqrt(det(gamma)) * (E, F_x, F_y, F_z).
```

For base `U_base`, a trial state is undensitized, checked for E/F
admissibility, given a closure and comoving moments, and evaluated with the
frozen neutrino rates. The residual is

```text
R_0   = tilde_E   - U_base[0]
        - dt_alpha * sqrt(det(gamma)) * S_E(trial)

R_i+1 = tilde_F_i - U_base[i+1]
        - dt_alpha * sqrt(det(gamma)) * S_i(trial).
```

`N`, `Gamma_N`, the number floor, and `N_source` are deliberately absent from
this residual. The endpoint current and backward-Euler number update occur
after E/F convergence. This separation is the current implementation
contract, not the reduced number-current design described in early planning
whitepapers.

## Finite-difference Jacobian

The 4-by-4 Jacobian is evaluated numerically from the residual. For each
unknown the solver uses the shared relative/absolute finite-difference step
policy in `ghl_m1_parameters`, first tries a forward perturbation, and uses a
backward one-sided perturbation if the forward trial fails specifically due to
admissibility. Other hard errors propagate. A second inadmissible one-sided
trial becomes an invalid-Jacobian failure.

Because the rates and primitives are frozen, a perturbation changes only the
trial E/F state, closure, comoving moments, and source evaluation. It does not
invoke a provider, EOS callback, or matter recovery. The public diagnostic
helpers for residual and Jacobian inspection are declared in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h).

## Newton step and acceptance

The shared driver:

1. evaluates the residual and its mixed absolute/relative weighted merit;
2. stops when the merit is at most one;
3. builds the finite-difference Jacobian and solves the 4-by-4 linear system;
4. backtracks a trial step while its residual and admissibility are tested; and
5. optionally projects a trial E/F vector through the shared realizability
   repair before retrying the residual.

A small Newton correction alone does not establish convergence. The published
iterate must be admissible and its weighted residual merit must be at most one.
An admissible trial above that threshold can continue Newton iteration if it
improves the merit. If no trial is accepted, the driver returns a retryable
solve failure without changing the caller output.

At an exact-zero-flux starting state, the neutrino solver may form an
admissible explicit E/F predictor as the initial Newton iterate. This is only
an initial guess; it does not alter the residual base, source formula, or
convergence normalization. There is no small-flux cutoff and no relaxed
tolerance for this case.

The shared controls and diagnostics are in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h): Newton iterations,
line-search backtracks, projection use, maximum residual, and weighted merit
are reported rather than silently discarded.

## Fallback schedule and endpoint

If a Newton attempt fails with a retryable solve/admissibility status, the
solver retries the complete local interval with the fixed substep schedule

```text
1, 2, 4, 8, 16 substeps.
```

Each substep starts from the prior accepted substep state and uses the same
frozen metric, primitives, and rates. On the first successful schedule, the
converged E/F state is undensitized. The solver then:

1. derives the endpoint number-current normalization;
2. applies backward-Euler number evolution using the endpoint `Gamma_N`;
3. repairs the complete neutrino state;
4. recomputes the endpoint current and optional mean-energy diagnostics;
5. checks configured endpoint bounds; and
6. assembles radiation and equal/opposite matter exchange.

The number update is

```text
N_out = (N_base + dt_alpha * eta_N)
        / (1 + dt_alpha * kappa_a_N / Gamma_N_endpoint).
```

`Gamma_N` is a number-current normalization, not the fluid Lorentz factor.
The charged-current lepton packet is formed from the endpoint charged-current
subset after the final state is known.

## Transactional failures

Outputs are initialized to the source base and a zero exchange packet. A hard
failure—invalid inputs, invalid rates, closure failure, invalid Jacobian,
unaccepted endpoint, or disallowed fallback—does not publish a partial
candidate. Exhausting the retry schedule returns the explicit terminal
no-update status and likewise leaves the source base and exchange packet
unchanged. Diagnostics may record the failure and retry path; they do not
make a failed state publishable.

The dispatcher separately rejects
`interaction_sources_already_applied` so an explicit interaction source cannot
be silently applied a second time. Hosts must treat the state and exchange as
a pair when deciding whether to publish.

## Current source authority and evidence

The public source-update declarations and policies are in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h); source ownership is specified
by the [M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md).
The focused tests are
[`unit_test_m1_fd_jacobian.c`](../../../Unit_Tests/unit_test_m1_fd_jacobian.c),
[`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c),
and [`unit_test_m1_error_handling.c`](../../../Unit_Tests/unit_test_m1_error_handling.c).
They provide local residual, solver, endpoint, and transactional checks; they
do not establish global timestep stability or full host-framework coupling.
