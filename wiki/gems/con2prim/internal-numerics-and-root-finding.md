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

Hybrid Noble 1D residual files live in
[`GRHayL/Con2Prim/Hybrid/Noble/Noble1D/`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble1D/):
`func_1D.c`, `func_Z.c`, `func_rho.c`, and source-present `func_rho2.c`.
The manifest builds the first three, not `func_rho2.c`. Built wrappers are
`hybrid_Noble1D.c` and `hybrid_Noble1D_entropy.c`; both initialize Noble
state, run `ghl_general_newton_raphson`, finalize primitives, then set
`diagnostics->speed_limited`, `diagnostics->n_iter`, and
`diagnostics->which_routine`.

Hybrid Noble 2D lives in
[`GRHayL/Con2Prim/Hybrid/Noble/Noble2D/`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble2D/).
`hybrid_Noble2D.c` uses `func_2D.c`, the same initialize/validate/Newton/finalize
helpers, and records `n_iter` and `which_routine` on success.

Tabulated Noble 2D lives in
[`GRHayL/Con2Prim/Tabulated/Noble2D/`](../../../GRHayL/Con2Prim/Tabulated/Noble2D/).
`tabulated_Noble2D.c` reuses `utils_Noble.h` and `ghl_general_newton_raphson`
but supplies tabulated-specific initialization, residual, table-bounds, and
finalization helpers before setting `n_iter` and `which_routine`.

## Palenzuela Internal Route

[`GRHayL/Con2Prim/utils_Palenzuela1D.h`](../../../GRHayL/Con2Prim/utils_Palenzuela1D.h)
declares the shared `compute_rho_W_from_x_and_conservatives` helper and the
hybrid/tabulated Palenzuela shared solver signatures. It depends on
[`GRHayL/Con2Prim/roots.h`](../../../GRHayL/Con2Prim/roots.h) for
`fparams_struct`, `roots_params`, and `ghl_brent`.

Hybrid Palenzuela files live in
[`GRHayL/Con2Prim/Hybrid/Palenzuela1D/`](../../../GRHayL/Con2Prim/Hybrid/Palenzuela1D/).
`hybrid_Palenzuela1D_energy.c` and `hybrid_Palenzuela1D_entropy.c` set
`diagnostics->which_routine`, then call the shared
`hybrid_Palenzuela1D.c` path with an energy or entropy EOS callback. The shared
path computes contractions through `ghl_compute_SU_Bsq_Ssq_BdotS`, brackets the
root, calls `ghl_brent`, records `n_iter`, computes utilde, and records
`speed_limited`.

Tabulated Palenzuela files live in
[`GRHayL/Con2Prim/Tabulated/Palenzuela1D/`](../../../GRHayL/Con2Prim/Tabulated/Palenzuela1D/).
The energy and entropy wrappers mirror the hybrid wrapper split, while
`tabulated_Palenzuela1D.c` adds table bounds, tabulated EOS calls, an optional
second Brent attempt from `T_min`, and the same `n_iter` and `speed_limited`
diagnostics ownership.

## Font And Newman Routes

[`GRHayL/Con2Prim/Hybrid/Font1D/`](../../../GRHayL/Con2Prim/Hybrid/Font1D/)
contains the hybrid Font path. `hybrid_Font1D.c` handles the public wrapper,
calls `hybrid_Font1D_loop.c` for the density iteration, computes the remaining
primitives, writes `diagnostics->speed_limited` when utilde limiting is called,
and sets `diagnostics->which_routine = ghl_con2prim_id_Font1D` on success.
Source search shows no `diagnostics->n_iter` assignment in this wrapper.

[`GRHayL/Con2Prim/Tabulated/Newman1D/`](../../../GRHayL/Con2Prim/Tabulated/Newman1D/)
contains the tabulated Newman energy and entropy paths. Both compute
`SU/Bsq/Ssq/BdotS` through `ghl_compute_SU_Bsq_Ssq_BdotS`, use local iterative
pressure updates with a maximum step count, set `diagnostics->n_iter = step`
inside the local helper, write `speed_limited` through utilde limiting, and set
`which_routine` in the public wrapper before the helper call.

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
  paths. Font1D does not show an `n_iter` write in its wrapper.
- `speed_limited`: written where solver finalization or utilde limiting calls
  `ghl_limit_utilde_and_compute_v`.

## Source-Present Drift Notes

`Noble1D_entropy2` has enum/name/declaration/source evidence in
[`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h),
[`GRHayL/Con2Prim/get_con2prim_routine_name.c`](../../../GRHayL/Con2Prim/get_con2prim_routine_name.c),
[`GRHayL/include/ghl_con2prim.h`](../../../GRHayL/include/ghl_con2prim.h), and
[`GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy2.c`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy2.c).
It is absent from
[`GRHayL/Con2Prim/Hybrid/Noble/Noble1D/make.code.defn`](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble1D/make.code.defn)
and has no selector case in
[`GRHayL/Con2Prim/con2prim_multi_method.c`](../../../GRHayL/Con2Prim/con2prim_multi_method.c).
Do not claim it as supported unless build-list and selector evidence change.
Its apparent companion `func_rho2.c` is also absent from the same manifest.
Manifest absence establishes configured-build exclusion only; it does not
establish whether maintainers intend integration, internalization, or removal.

`con2prim_CerdaDuran3D.cc` is source-present at
[`GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc`](../../../GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc),
but
[`GRHayL/Con2Prim/Tabulated/make.code.defn`](../../../GRHayL/Con2Prim/Tabulated/make.code.defn)
builds only the `Newman1D`, `Noble2D`, and `Palenzuela1D` subdirectories, and
the public selector in
[`GRHayL/Con2Prim/con2prim_multi_method.c`](../../../GRHayL/Con2Prim/con2prim_multi_method.c)
has no Cerda-Duran case. Treat this as source-present unresolved code, not
public support.

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
