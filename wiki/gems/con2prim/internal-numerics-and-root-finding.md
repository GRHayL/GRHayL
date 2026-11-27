# Con2Prim Internal Numerics And Root-Finding

> Internal Brent root-finding, Noble Newton-Raphson convergence and velocity-limit closure, Palenzuela `evolve_temp` callbacks, Font1D cycle detection, Newman iteration counters, and archival Cerda-Duran status. · Status: confirmed
> Up: [Con2Prim](index.md)

## Summary

Con2Prim solvers rely on four distinct internal numerical engines: `ghl_brent` for Palenzuela, `ghl_general_newton_raphson` for hybrid/tabulated Noble, damped fixed-point density iteration with Brent cycle detection for Font1D, and local iterative pressure updates for tabulated Newman1D. Each solver records `diagnostics->n_iter` with its own scope, and `con2prim_CerdaDuran3D.cc` remains unbuilt archival code.

## Detail

### Brent Root-Finder (`ghl_brent`) And Palenzuela Internal Route

- **`ghl_brent` (`roots.h`, `brent.c`):** records `routine_name` and caller bounds `a`, `b` in `roots_params` before any swap (preserving metadata on early returns), checks endpoint roots, rejects unbracketed intervals (`ghl_error_root_not_bracketed`), and rejects nonfinite bounds, residuals, or iterates (`ghl_error_c2p_singular`; width-based termination still requires a finite residual). Resets `rparams.n_iters = 0` per invocation; retries report only their own count.
- **Palenzuela (`utils_Palenzuela1D.h`, `Hybrid/Palenzuela1D/`, `Tabulated/Palenzuela1D/`):** wraps `ghl_compute_SU_Bsq_Ssq_BdotS`, root bracketing, `ghl_brent`, `utilde` limiting, and success-only `which_routine` writes; tabulated `tabulated_Palenzuela1D.c` adds table bounds and an optional second Brent attempt from `T_min`.
- **Tabulated Palenzuela `Y_e` and `evolve_temp` rules:** `tabulated_Palenzuela1D_energy.c` and `tabulated_Palenzuela1D_entropy.c` consume incoming `prims->Y_e` (never deriving `Y_e` from conserved `Y_e`; callers bypassing guesses must supply `Y_e`). When `params->evolve_temp` is true, callbacks invert energy/entropy for temperature; when false, both use bounded incoming temperature and the entropy callback ignores conserved entropy (though a failed-Brent `T_min` retry can still overwrite temperature). Neither Newman nor Noble shares this `evolve_temp` branch.

Evidence (descriptive behavior): `GRHayL/Con2Prim/brent.c::ghl_brent`; corroboration: `GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D.c`.

### Noble Newton Loop, Validation, And Stagnation Guards

`utils_Noble.h` and `Hybrid/Noble/` implement the shared Newton engine:
- **`initialize_Noble.c`:** copies `params->con2prim_max_iterations` and `con2prim_solver_tolerance`, rejects nonfinite/nonpositive conserved density (`ghl_error_neg_rho`), and builds scalar seeds (`Noble1D_entropy2` admits finite seeds up to `max_Lorentz_factor`; older `Noble1D_entropy` keeps its historical seed ceiling).
- **`general_newton_raphson.c` & `validate_x.c`:** runs per-invocation local budgets accumulating into `harm_aux->n_iter` (propagating nested density-solve errors before consumption). Convergence requires both relative step size (`|dx| / |x|`, absolute when `x == 0`) and recomputed residual objective `|f|` to satisfy `solver_tolerance`.
  For hybrid Noble2D's bounded squared-velocity coordinate, step size uses actual post-validation `x - x_old` and `func_2D.c` scales momentum residual by `(Bsq + Z)^2` and energy residual by `max(|Qdotn|, |Z|, Bsq)` when `ghl_validate_2D` clips at the upper bound (preventing false acceptance of a clipped iterate with large relative error; bound must match `hybrid_Noble2D.c`). Because each iteration's end-of-step residual feeds the next step, callbacks must be deterministic at a given `x` (tabulated Noble2D warm-starts from `harm_aux->temp_guess` within table root precision); an unconverged iterate with validated `x == x_old` returns `ghl_error_c2p_max_iter` immediately.

Evidence (descriptive behavior): `GRHayL/Con2Prim/Hybrid/Noble/general_newton_raphson.c::ghl_general_newton_raphson`; corroboration: `GRHayL/Con2Prim/Hybrid/Noble/validate_x.c::ghl_validate_2D`.

### Noble Closure Under Speed Limiting And Tabulated Noble2D

- **`finalize_Noble.c` & hybrid wrappers (`Hybrid/Noble/Noble1D/`, `Noble2D/`):** When ordinary hybrid Noble finalization limits velocity, it recomputes `rho = D / W_final` and enthalpy density `w = Z / W_final^2`, preserving exact `Z = rho h W_final^2` closure rather than mixing pre-limit `W` with post-limit `rho`. Near the cold-pressure boundary, subtraction amplification of roundoff in `w` can make pressure sign toolchain-sensitive and expose `ghl_error_neg_pressure`; configure a robust backup (such as `Font1D`) when marginal states must recover. `Noble1D_entropy` and `Noble1D_entropy2` (`func_rho2.c`) require `neos == 1` and `Gamma_th == Gamma_ppoly[0]`. Finalization writes `speed_limited` before the pressure check; `n_iter` and `which_routine` are written only after pressure passes.
- **Tabulated Noble2D (`Tabulated/Noble2D/tabulated_Noble2D.c`):** reuses `ghl_general_newton_raphson` with tabulated init/residual/finalization, rejects nonfinite inputs before bounding, corrects log-coordinate endpoint roundoff, returns recoverable EOS errors, and recomputes `w` with `W_final` when speed-limited.

Evidence (descriptive behavior): `GRHayL/Con2Prim/Hybrid/Noble/finalize_Noble.c`; corroboration: `GRHayL/Con2Prim/Tabulated/Noble2D/tabulated_Noble2D.c`.

### Font1D, Newman1D, Diagnostics Scope, And Archival Cerda-Duran

- **Hybrid Font1D (`Hybrid/Font1D/`):** evaluates cold `P`/`eps` at EOS-clamped density while dividing `P` by the recovered conservation-driven density in enthalpy closure (preserving conservative closure below `rho_min`). Resets `diagnostics->n_iter` at entry and accumulates outer `W`/momentum iterations across retries (reporting `0` on its `S_i gamma^ij S_j < 1e-300` shortcut after checking positive finite `D`; inner density iterations are uncounted). Each inner density loop uses Brent's exact cycle detection with geometric checkpoints (distinct from `ghl_brent`), rejecting nonconvergent cycles with `ghl_error_c2p_max_iter`, keeping only completed inner solves as retry seeds, and applying geometric damping only where log-derivative `L < -1/3` (preserving the fixed point and undamped residual check).
- **Tabulated Newman1D (`Tabulated/Newman1D/`):** uses local iterative pressure updates (`maxsteps = 300`), writes `diagnostics->n_iter = step` inside the local helper (so a `T_min` retry overwrites the prior attempt's count), and sets `which_routine` only after success.
- **Archival Cerda-Duran (`GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc`):** intentionally omitted from `Tabulated/make.code.defn`, `Doxyfile`, selectors, public headers, and tests; uses obsolete non-compiling API identifiers and supplies no C++ build/link evidence.

Evidence (descriptive behavior): `GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D_loop.c`; corroboration: `GRHayL/Con2Prim/Tabulated/Newman1D/tabulated_Newman1D_energy.c` and `GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc`.

## Sources

- [GRHayL/Con2Prim/roots.h](../../../GRHayL/Con2Prim/roots.h) - `fparams_struct`, `roots_params`, `ghl_brent` declaration
- [GRHayL/Con2Prim/brent.c](../../../GRHayL/Con2Prim/brent.c) - `ghl_brent` implementation
- [GRHayL/Con2Prim/utils_Noble.h](../../../GRHayL/Con2Prim/utils_Noble.h) - `harm_aux_vars_struct` and Noble helpers
- [GRHayL/Con2Prim/Hybrid/Noble/initialize_Noble.c](../../../GRHayL/Con2Prim/Hybrid/Noble/initialize_Noble.c) - Noble initialization and seed rules
- [GRHayL/Con2Prim/Hybrid/Noble/general_newton_raphson.c](../../../GRHayL/Con2Prim/Hybrid/Noble/general_newton_raphson.c) - shared Newton-Raphson engine
- [GRHayL/Con2Prim/Hybrid/Noble/validate_x.c](../../../GRHayL/Con2Prim/Hybrid/Noble/validate_x.c) - `ghl_validate_1D`, `ghl_validate_1D_entropy`, `ghl_validate_2D`
- [GRHayL/Con2Prim/Hybrid/Noble/finalize_Noble.c](../../../GRHayL/Con2Prim/Hybrid/Noble/finalize_Noble.c) - speed-limited closure `w = Z / W_final^2`
- [GRHayL/Con2Prim/Tabulated/Noble2D/tabulated_Noble2D.c](../../../GRHayL/Con2Prim/Tabulated/Noble2D/tabulated_Noble2D.c) - tabulated Noble2D solver
- [GRHayL/Con2Prim/utils_Palenzuela1D.h](../../../GRHayL/Con2Prim/utils_Palenzuela1D.h) - shared Palenzuela helper declarations
- [GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D.c](../../../GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D.c) - tabulated Palenzuela Brent driver and `T_min` retry
- [GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D_loop.c](../../../GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D_loop.c) - Font1D inner density loop, cycle detection, and damping
- [GRHayL/Con2Prim/Tabulated/Newman1D/tabulated_Newman1D_energy.c](../../../GRHayL/Con2Prim/Tabulated/Newman1D/tabulated_Newman1D_energy.c) - Newman1D pressure iteration and `n_iter`
- [GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc](../../../GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc) - unbuilt archival reference

## See Also

- Depends on: [Solver Matrix](solver-matrix.md)
- Depends on: [Recovery Flow](recovery-flow.md)
- Validated by: [Con2Prim Tests And Reference Data](tests-and-reference-data.md)
