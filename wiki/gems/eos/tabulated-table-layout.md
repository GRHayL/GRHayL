# Tabulated Table Layout

> Tabulated HDF5 EOS read sequence, structural acceptance checks, memory allocation/cleanup rules, energy shift, sound-speed adjustment, and `GHL_DISABLE_HDF5` guards. · Status: confirmed
> Up: [Equation Of State (EOS)](index.md)

## Summary

GRHayL loads HDF5 nuclear EOS tables via `ghl_eos_table_stellarcollapse`, validating grid dimensions, dataset element counts, uniform strictly increasing axes, and index representability before converting to code units, tabulating derived enthalpy (`table_logh`), and adjusting relativistic sound speed. Live tabulated `ghl_eos_parameters` objects must be freed via `NRPyEOS_free_memory` before reinitialization.

## Detail

### Table Keys, Read Sequence, And Structural Acceptance Checks

`NRPyEOS_keys` in `GRHayL/include/ghl_nrpyeos_tabulated.h` defines GRHayL's tabulated quantity slots (thermodynamics, sound speed, derivatives, chemical potentials, composition fractions, nuclear averages, adiabatic index, and derived enthalpy) along with `CODE_TO_CGS_*`/`CGS_TO_CODE_*` macros.

1. **Dispatch & read:** `NRPyEOS_read_table_set_EOS_params` (`NRPyEOS_read_table_set_EOS_params.c`) verifies `ghl_eos_tabulated`, dispatches `ghl_eos_table_stellarcollapse` to `NRPyEOS_stellarcollapse_read_table`, converts via `NRPyEOS_stellarcollapse_to_ghl`, and frees the temporary adapter struct.
2. **Loader acceptance checks:**
   - Every grid dimension is `>= 2`; total point product fits both `INT_MAX / NRPyEOS_ntablekeys` and `SIZE_MAX / sizeof(double) / NRPyEOS_ntablekeys`.
   - Every HDF5 read matches its exact element count (`1` for scalar dims, `energy_shift`, and optional `have_rel_cs2`; nonempty for axis and table arrays).
   - Stored axes and code-unit transforms are finite, strictly increasing, and uniform within `1e-8` cell-index units of the first interval (`check_axis`).
   - Derived physical endpoints, inverse spacings, and their products remain finite and positive (checked inline in `NRPyEOS_read_table_set_EOS_params`), and biased integer index expressions at the first and last nodes stay representable (`validate_biased_axis`), before publication.
3. **Post-read processing:** converts pressure, specific energy, sound speed, and derivatives to code units; populates `table_eps` by exponentiating the shifted code-unit `NRPyEOS_eps_key` entries of `table_all` (`NRPyEOS_eps_key + NRPyEOS_ntablekeys * i`, retaining `energy_shift`); computes derived enthalpy via `NRPyEOS_tabulate_enthalpy_checked` (`NRPyEOS_tabulate_enthalpy.c`, subtracting `energy_shift` and rejecting nonpositive/nonfinite `h`); runs `NRPyEOS_tabulated_adjust_sound_speed` (converting non-relativistic `cs2` using enthalpy when `have_rel_cs2` is absent/0, counting nonfinite/superluminal/negative values, and clamping only when `clean_sound_speed` is true); and stores table min/max bounds (subtracting `energy_shift` for `table_eps_min` and `table_eps_max`).

Evidence (descriptive behavior): `NRPyEOS_read_table_set_EOS_params.c::NRPyEOS_read_table_set_EOS_params`; corroboration: `NRPyEOS_tabulate_enthalpy.c::NRPyEOS_tabulate_enthalpy_checked` and `NRPyEOS_tabulated_adjust_sound_speed.c`.

Evidence (descriptive behavior): `NRPyEOS_stellarcollapse.c::check_axis`; corroboration: `NRPyEOS_read_table_set_EOS_params.c::validate_biased_axis`.

### Memory Allocation, Cleanup, And `GHL_DISABLE_HDF5`

Claim status: contested; contradiction: CONTR-0010. See [CONTR-0010](../../contradictions.md#contr-0010).

- **Fresh-destination precondition & transactional failure:** `ghl_initialize_tabulated_eos` (also reached through `ghl_initialize_tabulated_eos_functions_and_params`) requires a destination with no live table/NN allocations (bytes need not be pre-zeroed). Success transfers all six table arrays (`table_logrho`, `table_logT`, `table_Y_e`, `table_all`, `table_eps`, `table_logh`) and optional `c2p_nn`; failure frees the candidate and publishes an empty `ghl_eos_tabulated` aggregate with null pointers and zero dimensions.
- **`NRPyEOS_free_memory` (`NRPyEOS_free_memory.c`):** frees and nulls all six table arrays, cached beta-equilibrium arrays (via `NRPyEOS_tabulated_free_beq_quantities`), and `eos->c2p_nn` (via `ghl_c2p_nn_free`). Repeatable on successful, failed-empty, cleaned, or zero-initialized structs; does not free the outer struct or make uninitialized garbage safe. Call before reinitializing or switching EOS families.
- **Direct low-level reader (`NRPyEOS_read_table_set_EOS_params`):** calls that reach conversion require a zero-initialized or previously cleaned destination with no live table, beta-equilibrium, or NN allocations and with `eos_type`, `table_type`, and `clean_sound_speed` set, because conversion leaves the beta-equilibrium and NN pointers untouched and both teardown and rejected-table cleanup call `NRPyEOS_free_memory`; the caller owns the outer struct. Failures before conversion preserve pre-call fields; allocation/validation failures free and null the six table pointers while scalar fields need not roll back.
- **`GHL_DISABLE_HDF5`:** `configure --disable-hdf5` excludes `GRHayL/EOS/Tabulated/` sources; tabulated parameter initializers allocate nothing and return `ghl_error_used_disabled_hdf5` (process-wide tabulated pointers stay unassigned).

Evidence (descriptive behavior): `NRPyEOS_free_memory.c::NRPyEOS_free_memory`; corroboration: `ghl_nrpyeos_tabulated.h::NRPyEOS_read_table_set_EOS_params`.

## Sources

- [GRHayL/include/ghl_nrpyeos_tabulated.h](../../../GRHayL/include/ghl_nrpyeos_tabulated.h) - `NRPyEOS_keys`, unit macros, public tabulated declarations
- [GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c) - dispatch, validation, post-read conversions
- [GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c) - `check_axis` grid-uniformity tolerance
- [GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c) - `NRPyEOS_tabulate_enthalpy_checked`
- [GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c) - `have_rel_cs2` conversion and optional clamping
- [GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c) - repeatable cleanup of table, beta-eq, and NN allocations
- [GRHayL/EOS/Tabulated/make.code.defn](../../../GRHayL/EOS/Tabulated/make.code.defn) - compiled tabulated EOS sources

## See Also

- Depends on: [EOS Initialization And Dispatch](initialization-and-dispatch.md)
- Implements: [StellarCollapse Table Adapter](stellarcollapse-table-adapter.md)
- See also: [Tabulated Interpolation And Bounds](tabulated-interpolation-and-bounds.md)
- See also: [Contradictions](../../contradictions.md#contr-0010)
