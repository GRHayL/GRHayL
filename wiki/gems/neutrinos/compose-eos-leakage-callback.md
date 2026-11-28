# CompOSE EOS Leakage Callback

> Runtime `(rho, Y_e, T)` to `(T, n_b, Y_q)` mapping, six-output chemical-potential and composition formulas, safe callback execution, HDF5/registration boundaries, and validation checks. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

A CompOSE adapter serving `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` maps GRHayL `(rho, Y_e, T)` into CompOSE `(T, n_b, Y_q)` using an explicit mass-density reference $m_{\rm ref}^{\rho}$, builds nucleon potentials under one shared energy reference $m_{\rm ref}^{E}$, reconstructs rest-mass-inclusive `mu_e` from the table's lepton branch, extracts free `X_n` and `X_p` from particle indices `10` and `11`, and writes all six outputs only after every check succeeds.

## Detail

### Runtime State Mapping: Density And Charge Fraction

CompOSE tabulates $(T\,[\mathrm{MeV}],\,n_b\,[\mathrm{fm}^{-3}],\,Y_q)$, whereas GRHayL passes geometric `rho`, dimensionless `Y_e`, and `T` in MeV. Because CompOSE defines $n_b = \sum_i B_i n_i$ without fixing a baryon mass reference, declare $m_{\rm ref}^{\rho}$ explicitly (`rho_b = m_ref^rho * n_b`) and never infer it silently from the table's model neutron mass:

$$n_b\,[\mathrm{fm}^{-3}] = \frac{\rho_b\,[\mathrm{g\,cm}^{-3}]}{m_{\rm ref}^{\rho}\,[\mathrm{g}]} \times 10^{-39}.$$

To invert GRHayL geometric density back to cgs at the EOS callback, use the mathematical inverse of the live loader's `CGS_TO_CODE_DENSITY` (`NRPyEOS_stellarcollapse_to_ghl.c`), not the unused `CODE_TO_CGS_DENSITY` macro. Leakage C separately multiplies `rho` by `NRPyLeakage_units_geom_to_cgs_D`, which is not numerically identical to `1 / CGS_TO_CODE_DENSITY`; keep both conversions explicit unless performing a deliberate end-to-end unit renormalization.

CompOSE $Y_q$ is the charge fraction of strongly interacting particles. $Y_q = Y_e$ holds only for locally neutral matter whose sole charged massive lepton is the electron (or a baryonic-only table when the host adds electrons under electron-only neutrality). With muons, $Y_q = Y_e + Y_\mu$: the three-argument callback `(rho, Y_e, T)` cannot infer an independent evolved muon fraction, and current NRPyLeakage lacks muon charged-current channels and assumes equal heavy-lepton species.

Evidence (descriptive behavior): `NRPyEOS_stellarcollapse_to_ghl.c::CGS_TO_CODE_DENSITY`; corroboration: `ghl_nrpyleakage.h::NRPyLeakage_units_geom_to_cgs_D`.

### Building The Six Leakage Outputs

For the six-output ABI order, rest-mass-inclusive `mu_e` requirement, and free-nucleon `X_n`/`X_p` definitions, see [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md).

From decoded CompOSE conserved-charge potentials $(\mu_b, \mu_q, \mu_l)$, nucleon chemical potentials are $\mu_n^{\rm C} = \mu_b$ and $\mu_p^{\rm C} = \mu_b + \mu_q$. To match the O'Connor-Ott `EOSmaker` common-reference convention, subtract one explicit shared reference $m_{\rm ref}^{E}$ from both nucleons (e.g., Hempel SFHo `EOSmaker` sets $m_{\rm ref}^{E} =$ `amu_mev` and $m_{\rm ref}^{\rho} =$ `amu_cgs`, though other tables can choose differently and $m_{\rm ref}^{E}$ need not equal $m_{\rm ref}^{\rho}$):

$$\mu_n = \mu_b - m_{\rm ref}^{E}, \qquad \mu_p = \mu_b + \mu_q - m_{\rm ref}^{E}, \qquad \widehat\mu = \mu_n - \mu_p = -\mu_q.$$

Equivalently, from raw `eos.thermo`: $\mu_n = m_n(1 + Q_3) - m_{\rm ref}^{E}$, $\mu_p = m_n(1 + Q_3 + Q_4) - m_{\rm ref}^{E}$, $\widehat\mu = -m_n Q_4$; from generated identifiers `3`, `4`, `5` (`bshift`, `q`, `l`): $\mu_n = \texttt{bshift} + m_n - m_{\rm ref}^{E}$, $\mu_p = \texttt{bshift} + m_n + q - m_{\rm ref}^{E}$, $\widehat\mu = -q$. If an intermediate table removes each particle's own rest mass, restore $m_n$ and $m_p$ using the table's exact model masses before subtracting $m_{\rm ref}^{E}$. Subtracting $m_n$ from $\mu_n$ and $m_p$ from $\mu_p$ shifts $\widehat\mu$ by $-(m_n - m_p) \approx -1.293\,\mathrm{MeV}$; choosing the wrong common $m_{\rm ref}^{E}$ leaves $\widehat\mu$ and density-derived blocking unchanged but violates the six-output ABI for other callers. Do not inject ad hoc reaction shifts in the adapter.

In CompOSE tables (which never store neutrinos), the relativistic electron chemical potential is $\mu_e = \mu_{l_e} - \mu_q$ (never subtracting $m_{\rm ref}^{E}$ or $m_e$):

- Electron-only table with local neutrality (or electron-muon table documenting $\mu_{l_e} = \mu_{l_\mu}$): $\mu_e = \mu_l - \mu_q$ ($m_n(Q_5 - Q_4)$ raw; `l - q` generated).
- Electron-muon table with unequal lepton potentials: stored effective $\mu_l$ is a density-weighted average that cannot recover $\mu_e$; reject or supply extra data.
- No-lepton table ($I_l \ne 1$): CompOSE sets $\mu_l = 0$; supply $\mu_e$ from a thermodynamically consistent external electron-positron EOS (preserving Coulomb/phase assumptions), never interpreting zero $\mu_l$ as physical $\mu_e$.

For free nucleons ($B_i = 1$), extract $X_n = Y_n^{\rm C}$ and $X_p = Y_p^{\rm C}$ from CompOSE net-particle indices `10` (neutrons) and `11` (protons) assuming negligible antinucleons, rather than grouped nuclear fractions or $Y_q$. Map an omitted species at a grid node to zero only when format metadata establish that omission denotes zero.

Evidence (descriptive behavior): `compose_to_grhayl.py::_regularize_plane`; corroboration: `compose_to_grhayl.py::_read_raw_plane` and `unit_test_tabulated_eos_compose.c::check_six_value_interpolator`.

### Safe Callback Execution, Registration, And Grid Rules

Implement the callback in nine atomic steps using local temporaries so failures leave all six caller output pointers unchanged: (1) validate initialized EOS, non-null output pointers, positive finite `rho` and `T`, finite `Y_e`; (2) convert `rho` via `1 / CGS_TO_CODE_DENSITY` and $m_{\rm ref}^{\rho}$ to $n_b$; (3) map composition to $Y_q$; (4) reject out-of-range `T`, $n_b$, or $Y_q$; (5) interpolate thermodynamic and composition fields at one common state; (6) decode raw or generated fields once; (7) apply common-reference nucleon and lepton-branch formulas; (8) verify finite outputs and physical bounds; (9) commit all six outputs only if every step succeeds. Return truthful `ghl_error_codes_t` values (reusing matching range errors or adding a dedicated code rather than reusing an unrelated code).

Adding a native CompOSE backend inside GRHayL requires more than registering a callback in `NRPyEOS_initialize_tabulated_functions.c`: a non-HDF5 backend must also rework `GHL_DISABLE_HDF5` guards across `ghl_initialize_eos_functions`, tabulated initializers, table read/free/interpolation paths, leakage entry points, and `configure` source/test filters (GRHayLib separately requires HDF5). It must also wire an explicit `ghl_eos_table_t` enum entry, allocation/cleanup, and direct tests; select backends via `ghl_initialize_tabulated_eos` (since `ghl_initialize_tabulated_eos_functions_and_params` hard-codes `ghl_eos_table_stellarcollapse`).

Grid rules: offline StellarCollapse conversions must supply at least two strictly increasing points per axis (positive uniform log $n_b$ and $T$, uniform linear $Y_q$, matching `NRPyEOS_tabulated_helpers.h`); native backends may keep CompOSE's grids and official interpolation order when verified against `compose`; neither route may extrapolate or serve $T = 0$; axis/flattening order must be verified with asymmetric indices; and higher-order CompOSE interpolation will differ off-grid from GRHayL trilinear interpolation.

Evidence (descriptive behavior): `NRPyEOS_initialize_tabulated_functions.c::NRPyEOS_initialize_tabulated_functions`; corroboration: `NRPyEOS_tabulated_helpers.h::NRPyEOS_get_interp_spots`.

### Validation Steps And Unsupported Cases

Build an asymmetric golden-state set spanning low, middle, and high interior nodes plus off-grid points, declare tolerances up front, and use an independently produced StellarCollapse oracle only if it shares the exact raw table release, `compose` version, selection, interpolation policy, and density/chemical references:

1. Compare raw `Q3`-`Q5` decoding against generated `index_thermo` / `index_yi` decoding at identical grid nodes before interpolation.
2. Compare decoded node and off-grid values against official `compose` output under the chosen interpolation policy.
3. Verify $\mu_n^{\rm C} = \mu_b$, $\mu_p^{\rm C} = \mu_b + \mu_q$, and $\mu_n - \mu_p = -\mu_q$.
4. Where a matched oracle exists, compare all six callback outputs in ABI order plus `mu_e / T`, `muhat / T`, and `(mu_e - muhat) / T` (`mu_n`, `mu_p` as ABI outputs, not blocking inputs).
5. Verify $0 \le X_n, X_p \le 1$ and table-specific composition identities without assuming $X_n + X_p = 1$ when nuclei exist.
6. Validate electron-only and no-lepton closures separately; treat muon states as unsupported until the $Y_q$/API gap is resolved.
7. Verify lower/upper range errors, $T = 0$, nonfinite inputs, missing composition, unsupported lepton/charge mappings, and unchanged outputs on failure.
8. Check EOS-side (`CGS_TO_CODE_DENSITY`) and leakage-side (`NRPyLeakage_units_geom_to_cgs_D`) density conversions independently.
9. Compare all six opacity slots, `R_source`, `Q_source`, and all three luminosities at zero and nonzero optical depths (exercising blocking and bremsstrahlung across every EOS-dependent entry point).

Repo tests `Unit_Tests/compose/test_compose_to_grhayl.py`, `Unit_Tests/compose/test_compose_failure_coverage.py`, and `Unit_Tests/unit_test_tabulated_eos_compose.c` run on the synthetic `test-analytic` table in CI; `test_profile_facts` checks only the SRO-141 profile metadata, and a full SRO-141 replay needs an externally supplied table ([EOS Tests And Reference Data](../eos/tests-and-reference-data.md)). Unsupported without additional physics or data: cold/beta-equilibrated tables; muon-bearing states without $Y_q$ closure or with unequal lepton potentials; tables lacking free `X_n`/`X_p`; grouped nuclear substitutions; generated outputs missing $m_n, m_p, I_l$; out-of-bounds extrapolation; exotic/quark/hybrid opacity channels absent from nucleon-only leakage; unpropagated $m_{\rm ref}^{\rho}$ or $m_{\rm ref}^{E}$ changes; raw $Q_7$ treated as `eps`; and undocumented web endpoints.

Evidence (descriptive behavior): `test_compose_to_grhayl.py::test_conversion_writes_regularized_stellarcollapse_schema`; corroboration: `test_compose_failure_coverage.py::test_helper_rejections`.

## Sources

- [GRHayL/include/ghl_eos_functions_declaration.h](../../../GRHayL/include/ghl_eos_functions_declaration.h) - `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`
- [GRHayL/include/ghl_nrpyleakage.h](../../../GRHayL/include/ghl_nrpyleakage.h) - `NRPyLeakage_units_geom_to_cgs_D`
- [GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c) - `NRPyEOS_initialize_tabulated_functions`
- [GRHayL/EOS/Tabulated/interpolators/NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c) - current six-field interpolator
- [GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c) - `CGS_TO_CODE_DENSITY`
- [GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h) - `NRPyEOS_get_interp_spots`
- [tools/compose/compose_to_grhayl.py](../../../tools/compose/compose_to_grhayl.py) - CompOSE raw decoding and common-reference chemical-potential mapping
- [Unit_Tests/compose/test_compose_to_grhayl.py](../../../Unit_Tests/compose/test_compose_to_grhayl.py) - offline converter verification
- [Unit_Tests/compose/test_compose_failure_coverage.py](../../../Unit_Tests/compose/test_compose_failure_coverage.py) - converter failure-path verification
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c) - C runtime integration test
- [StellarCollapse EOS Page](https://stellarcollapse.org/equationofstate.html) - target table conventions
- [EOSmaker SVN r13 Archive](https://stellarcollapse.org/EOS/EOSmaker_svn_r13.tar.gz) - `src/fix_table_units.F90` common-reference mapping
- [EOSmaker fix_table_units.F90](https://github.com/evanoconnor/EOSmaker/blob/master/src/fix_table_units.F90) - common nucleon reference subtraction
- [EOSmaker Hempel SFHo Config](https://github.com/evanoconnor/EOSmaker/blob/master/src/specific_eos/eos_table_Hempel_sc_sfho.F90) - `amu_mev` and `amu_cgs` selection

## See Also

- Depends on: [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md)
- Depends on: [CompOSE EOS Conversion And Regularization](compose-eos-conversion-and-regularization.md)
- Depends on: [Tabulated Table Layout](../eos/tabulated-table-layout.md)
- Depends on: [Tabulated Interpolation And Bounds](../eos/tabulated-interpolation-and-bounds.md)
- See also: [Build And Configure](../../architecture/build-and-configure.md)
- See also: [Physics Model And Sources](physics-model-and-sources.md)
