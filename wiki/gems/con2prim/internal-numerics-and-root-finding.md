# Con2Prim Internal Numerics And Root-Finding

Purpose: route internal numerical helpers used by Con2Prim solvers without
changing solver support status. Public dispatch and method support remain in
[solver matrix](solver-matrix.md); recovery order remains in
[recovery flow](recovery-flow.md); limits and conversion helpers remain in
[limits and conversions](limits-and-conversions.md). Repo-local source remains
authority.

## Public API Versus Internal Helpers

Public Con2Prim entry points and diagnostics types are declared in
[`GRHayL/include/ghl_con2prim.h`](../../../GRHayL/include/ghl_con2prim.h), with
shared enums and parameter fields in
[`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h). Public callers should
route through the selector, multi-method, pre-limit, post-limit, and conversion
functions declared there.

The files below are internal numerical routes. They are built because they
appear in Con2Prim build lists, are included by solver sources, or are source
present but not built as noted below:

- Shared build list:
  [`GRHayL/Con2Prim/make.code.defn`](../../../GRHayL/Con2Prim/make.code.defn).
- Hybrid build lists:
  [`GRHayL/Con2Prim/Hybrid/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/make.code.defn),
  [`GRHayL/Con2Prim/Hybrid/Noble/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/Noble/make.code.defn),
  [`GRHayL/Con2Prim/Hybrid/Noble/Noble1D/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble1D/make.code.defn),
  [`GRHayL/Con2Prim/Hybrid/Noble/Noble2D/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble2D/make.code.defn),
  [`GRHayL/Con2Prim/Hybrid/Palenzuela1D/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/Palenzuela1D/make.code.defn), and
  [`GRHayL/Con2Prim/Hybrid/Font1D/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/Font1D/make.code.defn).
- Tabulated build lists:
  [`GRHayL/Con2Prim/Tabulated/make.code.defn`](../../../GRHayL/Con2Prim/Tabulated/make.code.defn),
  [`GRHayL/Con2Prim/Tabulated/Newman1D/make.code.defn`](../../../GRHayL/Con2Prim/Tabulated/Newman1D/make.code.defn),
  [`GRHayL/Con2Prim/Tabulated/Noble2D/make.code.defn`](../../../GRHayL/Con2Prim/Tabulated/Noble2D/make.code.defn), and
  [`GRHayL/Con2Prim/Tabulated/Palenzuela1D/make.code.defn`](../../../GRHayL/Con2Prim/Tabulated/Palenzuela1D/make.code.defn).

## Brent And Root Helper Route

[`GRHayL/Con2Prim/roots.h`](../../../GRHayL/Con2Prim/roots.h) defines the
internal `fparams_struct`, `roots_params`, `roots_info`, `swap`, `sign`, and
the `ghl_brent` prototype. `fparams_struct` carries Palenzuela root-function
state and EOS callback pointers; `roots_params` carries root-finder bounds,
tolerance, maximum iterations, iteration count, root, and residual.

[`GRHayL/Con2Prim/brent.c`](../../../GRHayL/Con2Prim/brent.c) implements
`ghl_brent`. The helper checks endpoint roots, rejects unbracketed intervals
with `ghl_error_root_not_bracketed`, records `routine_name`, `a`, and `b` in
`roots_params`, then iterates by Brent's interpolation/bisection logic until
success or `ghl_error_c2p_max_iter`.

Current built users are the Palenzuela shared solvers:
[`GRHayL/Con2Prim/Hybrid/Palenzuela1D/hybrid_Palenzuela1D.c`](../../../GRHayL/Con2Prim/Hybrid/Palenzuela1D/hybrid_Palenzuela1D.c)
and
[`GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D.c`](../../../GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D.c).
Both fill `roots_params`, call `ghl_brent`, and copy `rparams.n_iters` into
`diagnostics->n_iter` on success.

## Shared Magnetic And Momentum Contractions

[`GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c`](../../../GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c)
implements `ghl_compute_SU_Bsq_Ssq_BdotS`. It takes ADM metric data,
undensitized conservatives, and primitives with `BU`; computes raised momentum
`SU`, `Bsq`, `Ssq`, and `BdotS`; and may rescale a local copy of momentum
before contraction. It does not mutate caller conservative storage.

Built callers include
[`GRHayL/Con2Prim/guess_primitives.c`](../../../GRHayL/Con2Prim/guess_primitives.c),
Palenzuela shared solvers under
[`GRHayL/Con2Prim/Hybrid/Palenzuela1D/`](../../../GRHayL/Con2Prim/Hybrid/Palenzuela1D/)
and
[`GRHayL/Con2Prim/Tabulated/Palenzuela1D/`](../../../GRHayL/Con2Prim/Tabulated/Palenzuela1D/),
and Newman wrappers under
[`GRHayL/Con2Prim/Tabulated/Newman1D/`](../../../GRHayL/Con2Prim/Tabulated/Newman1D/).

## Noble Internal Route

[`GRHayL/Con2Prim/utils_Noble.h`](../../../GRHayL/Con2Prim/utils_Noble.h)
declares `harm_aux_vars_struct`, pressure/velocity helpers, Newton prototypes,
Noble initialization/finalization helpers, validation helpers, and residual
function prototypes.

Core Noble helper files are built from
[`GRHayL/Con2Prim/Hybrid/Noble/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/Noble/make.code.defn):

- [`initialize_Noble.c`](../../../GRHayL/Con2Prim/Hybrid/Noble/initialize_Noble.c)
  sets `harm_aux.n_iter`, copies `params->con2prim_max_iterations` and
  `params->con2prim_solver_tolerance`, forms HARM-style contractions, and
  produces scalar guesses for ordinary and entropy Noble paths.
- [`general_newton_raphson.c`](../../../GRHayL/Con2Prim/Hybrid/Noble/general_newton_raphson.c)
  runs the shared Newton loop, calls the supplied residual and validate
  callbacks, increments `harm_aux->n_iter`, and returns success,
  `ghl_error_c2p_singular`, or `ghl_error_c2p_max_iter`.
- [`validate_x.c`](../../../GRHayL/Con2Prim/Hybrid/Noble/validate_x.c)
  contains `ghl_validate_1D`, `ghl_validate_1D_entropy`, and
  `ghl_validate_2D` constraints for Newton updates.
- [`finalize_Noble.c`](../../../GRHayL/Con2Prim/Hybrid/Noble/finalize_Noble.c)
  converts solved scalars into primitives and returns whether
  `ghl_limit_utilde_and_compute_v` speed-limited the result.

When ordinary hybrid Noble finalization limits the velocity, it recomputes
`rho = D/W_final` and `w = Z/W_final^2`. This preserves the defining closure
`Z = rho h W_final^2` for the limited state instead of mixing the pre-limit
Lorentz factor with the limited density. Near the cold-pressure boundary, the
consistent closure can expose a roundoff-sized nonpositive pressure. The Noble
wrapper then returns `ghl_error_neg_pressure`, allowing configured backup
recovery; it does not clip conservative energy or impose a tabulated-EOS floor.

Hybrid Noble 1D residual files live in
[`GRHayL/Con2Prim/Hybrid/Noble/Noble1D/`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble1D/):
`func_1D.c`, `func_Z.c`, `func_rho.c`, and `func_rho2.c`. The manifest builds
all four plus `hybrid_Noble1D.c`, `hybrid_Noble1D_entropy.c`, and
`hybrid_Noble1D_entropy2.c`. The entropy2 path solves the momentum equation
directly for density. Its `Z(rho)` and analytic derivative use the active
piecewise-polytropic cold pressure, cold-energy integration constant, and
thermal Gamma. Each wrapper initializes Noble state, runs
`ghl_general_newton_raphson`, and finalizes primitives. Finalization ORs
`diagnostics->speed_limited` before the wrapper's pressure gate;
`diagnostics->n_iter` and `diagnostics->which_routine` are set only after that
gate succeeds.

Hybrid Noble 2D lives in
[`GRHayL/Con2Prim/Hybrid/Noble/Noble2D/`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble2D/).
`hybrid_Noble2D.c` uses `func_2D.c`, the same initialize/validate/Newton/finalize
helpers, ORs `speed_limited` during finalization, and records `n_iter` and
`which_routine` only after the pressure gate succeeds.

Tabulated Noble 2D lives in
[`GRHayL/Con2Prim/Tabulated/Noble2D/`](../../../GRHayL/Con2Prim/Tabulated/Noble2D/).
`tabulated_Noble2D.c` reuses `utils_Noble.h` and `ghl_general_newton_raphson`
but supplies tabulated-specific initialization, residual, table-bounds, and
finalization helpers. Finalization may update `speed_limited`; `n_iter` and
`which_routine` remain success-only writes after the pressure gate.

## Palenzuela Internal Route

[`GRHayL/Con2Prim/utils_Palenzuela1D.h`](../../../GRHayL/Con2Prim/utils_Palenzuela1D.h)
declares the shared `compute_rho_W_from_x_and_conservatives` helper and the
hybrid/tabulated Palenzuela shared solver signatures. It depends on
[`GRHayL/Con2Prim/roots.h`](../../../GRHayL/Con2Prim/roots.h) for
`fparams_struct`, `roots_params`, and `ghl_brent`.

Hybrid Palenzuela files live in
[`GRHayL/Con2Prim/Hybrid/Palenzuela1D/`](../../../GRHayL/Con2Prim/Hybrid/Palenzuela1D/).
`hybrid_Palenzuela1D_energy.c` and `hybrid_Palenzuela1D_entropy.c` call the
shared `hybrid_Palenzuela1D.c` path with an energy or entropy EOS callback and
set `diagnostics->which_routine` only after it succeeds. The shared path
computes contractions through `ghl_compute_SU_Bsq_Ssq_BdotS`, brackets the
root, calls `ghl_brent`, records `n_iter`, computes utilde, and records
`speed_limited`.

Tabulated Palenzuela files live in
[`GRHayL/Con2Prim/Tabulated/Palenzuela1D/`](../../../GRHayL/Con2Prim/Tabulated/Palenzuela1D/).
The energy and entropy wrappers mirror the hybrid success-only diagnostic
assignment, while
`tabulated_Palenzuela1D.c` adds table bounds, tabulated EOS calls, an optional
second Brent attempt from `T_min`, and the same `n_iter` and `speed_limited`
diagnostics ownership.

## Font And Newman Routes

[`GRHayL/Con2Prim/Hybrid/Font1D/`](../../../GRHayL/Con2Prim/Hybrid/Font1D/)
contains the hybrid Font path. `hybrid_Font1D.c` handles the public wrapper,
calls `hybrid_Font1D_loop.c` for the density iteration, computes the remaining
primitives, ORs `diagnostics->speed_limited` when utilde limiting is called,
and sets `diagnostics->which_routine = ghl_con2prim_id_Font1D` on success.
It resets `diagnostics->n_iter` at entry and accumulates density-loop iterations;
the shortcut reports zero only when the conservative momentum norm satisfies
`S_i gamma^ij S_j < 1e-300`.

[`GRHayL/Con2Prim/Tabulated/Newman1D/`](../../../GRHayL/Con2Prim/Tabulated/Newman1D/)
contains the tabulated Newman energy and entropy paths. Both compute
`SU/Bsq/Ssq/BdotS` through `ghl_compute_SU_Bsq_Ssq_BdotS`, use local iterative
pressure updates with a maximum step count, set `diagnostics->n_iter = step`
inside the local helper, write `speed_limited` through utilde limiting, and set
`which_routine` in the public wrapper only after the retry sequence succeeds.

## Diagnostics Ownership

Diagnostics layout and public initialization live in
[`GRHayL/include/ghl_con2prim.h`](../../../GRHayL/include/ghl_con2prim.h) and
[`GRHayL/Con2Prim/initialize_diagnostics.c`](../../../GRHayL/Con2Prim/initialize_diagnostics.c).
For this internal page, only source-proven writes are routed:

- `which_routine`: set by successful solver wrappers in Noble, Font,
  Palenzuela, and Newman paths. Multi-method backup selection is routed in
  [recovery flow](recovery-flow.md).
- `n_iter`: set from `harm_aux.n_iter` in built Noble paths, from
  `rparams.n_iters` in built Palenzuela paths, and from local `step` in Newman
  paths. Font1D reports its invocation-local accumulated density iterations.
- `speed_limited`: sticky across every attempted solver finalization or utilde
  limiting call, including attempts that later fail.

## Archival Source

`con2prim_CerdaDuran3D.cc` is source-present at
[`GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc`](../../../GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc),
with a prominent archival header. It is retained only to preserve potentially
useful numerical work. It uses identifiers incompatible with the current API;
[`GRHayL/Con2Prim/Tabulated/make.code.defn`](../../../GRHayL/Con2Prim/Tabulated/make.code.defn)
intentionally builds the shared tabulated guess helper, the
`neural_network_guess` support, and the `Newman1D`, `Noble2D`, and
`Palenzuela1D` subdirectories, but not the archival file. The
public selector in
[`GRHayL/Con2Prim/con2prim_multi_method.c`](../../../GRHayL/Con2Prim/con2prim_multi_method.c)
has no Cerda-Duran case. There is no public declaration, method ID, active
GRHayLib keyword, or test. Treat it as unbuilt, unsupported archival source.

The configured source inventory is manifest-driven and currently contains no
`.cc` entry for this file. File presence alone therefore supplies neither C++
compilation nor linkage evidence.

## Read-Only Evidence Routes

- Doxygen method overview:
  [`docs/raw/Con2Prim.dox`](../../../docs/raw/Con2Prim.dox).
- Hybrid selected-method and fixture evidence:
  [`Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`](../../../Unit_Tests/unit_test_con2prim_multi_method_hybrid.c).
- Tabulated method, backup, and fixture evidence:
  [`Unit_Tests/unit_test_con2prim_tabulated.c`](../../../Unit_Tests/unit_test_con2prim_tabulated.c).
- Hybrid failure and backup evidence:
  [`Unit_Tests/unit_test_hybrid_failure.c`](../../../Unit_Tests/unit_test_hybrid_failure.c).
- Default runner fixture routes:
  [`.github/run_tests.sh`](../../../.github/run_tests.sh).
