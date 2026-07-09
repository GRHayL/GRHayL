# Con2Prim Gem

Compact hub for conservative-to-primitive recovery, dispatch, helper contracts,
tests, and public docs. Keep detailed method status and equations out of this
page; source, raw Doxygen, and tests remain authority.

## Read First

- Method support and dispatch: [solver matrix](con2prim/solver-matrix.md)
- Recovery order, backups, diagnostics: [recovery flow](con2prim/recovery-flow.md)
- Limits, guesses, and conversions: [limits and conversions](con2prim/limits-and-conversions.md)
- Tabulated neural-network primitive-guess routing:
  [neural-network primitive guess](con2prim/neural-network-primitive-guess.md)
- Internal numerics and root-finding:
  [internal numerics and root-finding](con2prim/internal-numerics-and-root-finding.md)
- Tests, fixtures, and HDF5 caveats: [tests and fixtures](con2prim/tests-and-fixtures.md)
- Public docs: `docs/raw/Con2Prim.dox`
- Physics background: `docs/raw/derivation.md`

## Source Authority

Treat a Con2Prim method as supported only when local evidence agrees across
`ghl_con2prim_id_t`/name mapping, public declaration, selector dispatch, build
lists, and tests where applicable. File presence alone is not support.

## Public Surface

- Headers: `GRHayL/include/ghl_con2prim.h`, `GRHayL/include/ghl.h`
- Drivers: `ghl_con2prim_hybrid_multi_method`,
  `ghl_con2prim_tabulated_multi_method`
- Selectors: `ghl_con2prim_hybrid_select_method`,
  `ghl_con2prim_tabulated_select_method`
- Helpers: `ghl_apply_conservative_limits`,
  `ghl_undensitize_conservatives`, `ghl_guess_primitives`,
  `ghl_enforce_primitive_limits_and_compute_u0`, `ghl_compute_conservs`,
  `ghl_compute_conservs_and_Tmunu`
- Tabulated NN primitive-guess API: `GHL_NN_C2P_API_VERSION`,
  `ghl_c2p_nn_model`, `ghl_c2p_nn_guess`,
  `ghl_c2p_nn_guess_primitives`, `ghl_c2p_nn_validate_model`,
  `ghl_c2p_nn_load_hdf5`, `ghl_c2p_nn_load_from_eos_hdf5`,
  `ghl_c2p_nn_free`

## Source Inventory

- Shared recovery and helper routines: `GRHayL/Con2Prim/`
- Hybrid solvers: `GRHayL/Con2Prim/Hybrid/`
- Tabulated solvers: `GRHayL/Con2Prim/Tabulated/`
- Tabulated NN primitive guesses:
  `GRHayL/Con2Prim/Tabulated/neural_network_guess/`
- Internal Brent/root helpers, Noble/Palenzuela utilities, Font/Newman helper
  paths, and magnetic/momentum contractions route through
  [internal numerics and root-finding](con2prim/internal-numerics-and-root-finding.md).
- Build lists: `GRHayL/Con2Prim/**/make.code.defn`
- Tests and data: `Unit_Tests/unit_test_con2prim_*.c`,
  `Unit_Tests/unit_test_apply_conservative_limits.c`,
  `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`,
  `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`,
  `.github/run_tests.sh`

## Key Contracts

- Initialize `ghl_con2prim_diagnostics` with `ghl_initialize_diagnostics`.
- Track densitized versus undensitized conservative inputs before calling
  selectors or helpers.
- Multi-method drivers may reset primitive guesses and try up to three backup
  routines.
- Enforce primitive bounds and compute `u0` before using recovered primitives.
- Tabulated recovery is HDF5-gated and depends on table-backed EOS setup.
