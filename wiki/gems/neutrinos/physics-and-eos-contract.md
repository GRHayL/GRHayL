# NRPyLeakage Physics And EOS Contract

## Purpose

This page preserves the physical meaning behind the generated NRPyLeakage C
expressions and the minimum EOS contract needed to use them. Current GRHayL
source, headers, and tests remain the implementation authority. The original
NRPy+ notebooks and the accompanying paper remain derivation evidence; they do
not guarantee byte-for-byte regeneration of the current GRHayL files.

## Physical Model

NRPyLeakage is a grey neutrino-leakage model. It computes local free emission
rates, estimates diffusion suppression from energy-averaged opacities and
optical depths, and interpolates between free-streaming and diffusion-limited
rates. It is not a Boltzmann, Monte Carlo, or moment-transport solver.

The generated rates include:

- electron and positron capture on nucleons (beta processes);
- electron-positron pair annihilation;
- transverse plasmon decay; and
- nucleon-nucleon bremsstrahlung.

The opacity model includes neutrino absorption on nucleons for electron
neutrinos and antineutrinos, plus neutrino scattering on neutrons and protons.
Heavy-lepton species have scattering opacity but no beta-process absorption in
this model. The implementation does not include neutrino momentum deposition
or absorption-driven heating in optically thin ejecta.

## Species Contract

The public structs order species as `nue`, `anue`, and `nux`:

- `nue`: electron neutrino, $\nu_e$;
- `anue`: electron antineutrino, $\bar\nu_e$;
- `nux`: one representative heavy-lepton species, any one of
  $\nu_\mu,\bar\nu_\mu,\nu_\tau,\bar\nu_\tau$.

The formulation assumes all four heavy-lepton species contribute equally.
Consequently, `nux` rates and `lum.nux` represent one species, while the matter
energy source explicitly multiplies the effective `nux` cooling rate by four.
No heavy-lepton number rate enters the electron-fraction source.

Ground truth:

- [`GRHayL/include/ghl_radiation.h`](../../../GRHayL/include/ghl_radiation.h)
- [`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
- [`NRPyLeakage_compute_neutrino_luminosities.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c)

## EOS Inputs And Meanings

The three EOS-dependent leakage routines call
`ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(eos, rho, Y_e, T, ...)`.
After the input arguments, the output-pointer order is `muhat`, `mu_e`,
`mu_p`, `mu_n`, `X_n`, `X_p`, as declared in
[`ghl_eos_functions_declaration.h`](../../../GRHayL/include/ghl_eos_functions_declaration.h).
The returned values mean:

| Value | Required meaning | Direct use in current leakage source |
| --- | --- | --- |
| `mu_e` | Electron chemical potential, including electron rest-mass energy for the current O'Connor-Ott/StellarCollapse table contract | Electron degeneracy `mu_e/T`, equilibrium neutrino degeneracy, beta rates, pair/plasmon rates |
| `mu_n` | Neutron chemical potential | Neutron degeneracy/blocking factors |
| `mu_p` | Proton chemical potential | Proton degeneracy/blocking factors |
| `muhat` | Neutron-proton chemical-potential difference; GRHayL's StellarCollapse adapter defines `muhat = mu_n - mu_p` | `muhat/T` in nucleon fraction factors and `(mu_e-muhat)/T` in equilibrium electron-neutrino degeneracy |
| `X_n` | Free-neutron mass fraction | Nucleon-nucleon bremsstrahlung composition factor |
| `X_p` | Free-proton mass fraction | Nucleon-nucleon bremsstrahlung composition factor |

`X_n` and `X_p` are free-nucleon mass fractions, not the total neutron and
proton number fractions. Current source separately sets `Y_p = Y_e` and
`Y_n = 1-Y_e` for charged-current/scattering factors. An EOS adapter must not
substitute `Y_n`/`Y_p` for `X_n`/`X_p`, particularly where nuclei or other
bound species are present.

The chemical potentials and temperature must share MeV units because current
source forms their ratios without conversion. `Y_e`, `X_n`, and `X_p` are
dimensionless. The source performs no chemical-potential zero-point or
rest-mass adjustment after the EOS call.

For the O'Connor-Ott/StellarCollapse table format read by GRHayL, the format
documentation specifies `mu_e` including electron rest-mass energy. GRHayL
reads and interpolates that field without a zero-point adjustment before the
leakage callback returns it. A replacement EOS that supplies a rest-mass-
subtracted electron chemical potential must map it to this table convention
before the leakage code forms degeneracy ratios. This is a table-format
contract, not a universal chemical-potential convention.

Ground truth:

- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h)
- [`GRHayL/include/ghl_eos_functions.h`](../../../GRHayL/include/ghl_eos_functions.h)
- [`GRHayL/EOS/Tabulated/interpolators/NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c`](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c)
- [`GRHayL/Neutrinos/NRPyLeakage/`](../../../GRHayL/Neutrinos/NRPyLeakage/)

### Chemical-Potential Convention Hazard

An alternative EOS must reproduce the combinations consumed by the code, not
only similarly named outputs. In current source the equilibrium degeneracies
are

$$
\eta_{\nu_e}^{\rm eq}=\frac{\mu_e-\widehat\mu}{T},\qquad
\eta_{\bar\nu_e}^{\rm eq}=-\eta_{\nu_e}^{\rm eq},
$$

with `muhat` consumed directly as $\widehat\mu$. The original implementation
notebook's displayed derivation also writes this in a convention with an
explicit neutron-proton rest-mass gap $Q$, while its generator code sets
`eta_hat = muhat/T`. GRHayL's StellarCollapse adapter reads the table's
`muhat` field and defines it as `mu_n - mu_p`; it does not subtract `Q` in the
leakage routine.

Therefore, do not independently add or subtract $Q$, or reconstruct
`muhat` from differently zeroed `mu_n` and `mu_p`, without first mapping the
EOS convention to the current table convention. The StellarCollapse table
documentation recommends
$\mu_{\nu_e}^{\rm eq}=\mu_e-\mu_n+\mu_p=\mu_e-\widehat\mu$, and warns that
some LS tables' stored `munu` field has an erroneous shift. The supported
table-format meaning is nevertheless `munu = mu_e - mu_n + mu_p`.

Current NRPyLeakage does not consume stored `munu`; it forms the equilibrium
combination from `mu_e` and `muhat`. Elsewhere,
[`NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c)
does consume stored `munu` values to locate a zero crossing. A correct header
comment does not alter table data, so the upstream LS stored-value caveat still
applies. Leakage adapters should supply the six named quantities above instead
of substituting the table's `munu` field.

`ghl_nrpyleakage.h` declares `NRPyLeakage_Q_npmass` and
`NRPyLeakage_ZL_Q_npmass`, but none of the five current leakage C files uses
either constant. The current generated blocks consume table `muhat` directly;
do not infer an implicit neutron-proton mass-gap correction from those
declarations or insert one into `(mu_e-muhat)/T` without changing and
revalidating the formulation. These constants are unrelated to `Q_source`.

Practical adapter validation should compare all six returned quantities and
the combinations `mu_e/T`, `mu_n/T`, `mu_p/T`, and `(mu_e-muhat)/T` against a
known StellarCollapse-table state before comparing final rates. Merely matching
pressure and internal energy does not validate this interface.

## Number And Energy Slots

`ghl_neutrino_opacities` declares two unnamed entries per species. Their
meaning is established by the generated implementation and notebooks, not by
type-level names in `ghl_radiation.h`:

- `[0]`: $j=0$, neutrino-number transport; for `nue` and `anue`, used for
  equilibrium-degeneracy interpolation, number diffusion time, and
  `R_source`. The `nux[0]` opacity is emitted and its optical depth is updated,
  but current source and luminosity routines do not consume `nux[0]`;
- `[1]`: $j=1$, neutrino-energy transport; used for energy diffusion time,
  `Q_source`, and luminosity suppression.

This mapping is a current implementation contract, not a generic promise that
every two-entry opacity container has the same meaning. Preserve `[0]`/`[1]`
at the public API boundary and document named aliases in downstream adapters.
Do not swap the slots: they use different Fermi-Dirac energy moments and
therefore generally differ numerically.

For species $\nu_i$ and slot $j$, total transport opacity is an inverse
mean free path. Schematically,

$$
\kappa_{t,j}^{\nu_e}=\kappa_{s,j}^{\nu_e n}
+\kappa_{s,j}^{\nu_e p}+\kappa_{a,j}^{\nu_e n},
$$

$$
\kappa_{t,j}^{\bar\nu_e}=\kappa_{s,j}^{\bar\nu_e n}
+\kappa_{s,j}^{\bar\nu_e p}+\kappa_{a,j}^{\bar\nu_e p},
\qquad
\kappa_{t,j}^{\nu_x}=\kappa_{s,j}^{\nu_x n}
+\kappa_{s,j}^{\nu_x p}.
$$

Optical depth is dimensionless and is the path integral
$\tau_j^{\nu_i}=\int ds\,\kappa_{t,j}^{\nu_i}$. The local GRHayL path-of-
least-resistance routine performs one six-face-neighbor update and no diagonal
integration. It does not own global convergence, outer-boundary initialization,
or AMR synchronization; callers must provide those pieces.

Ground truth:

- [`GRHayL/include/ghl_radiation.h`](../../../GRHayL/include/ghl_radiation.h)
- [`NRPyLeakage_optical_depths_PathOfLeastResistance.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c)
- [`Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`](../../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c)

## Effective Rates

Ideally, for each species and transport moment, leakage suppresses the total
free rate with its diffusion-to-loss-time ratio:

$$
R_{\nu_i}^{\rm eff}=\frac{R_{\nu_i}^{\rm free}}
 {1+t_{\nu_i,0}^{\rm diff}/t_{\nu_i,0}^{\rm loss}},\qquad
Q_{\nu_i}^{\rm eff}=\frac{Q_{\nu_i}^{\rm free}}
 {1+t_{\nu_i,1}^{\rm diff}/t_{\nu_i,1}^{\rm loss}}.
$$

The checked-in implementation has the species and unit exceptions documented
in [Current Heavy-Lepton Exception](#current-heavy-lepton-exception) and
[Current Suppression-Ratio Unit Seam](#current-suppression-ratio-unit-seam).

The ideal geometric-unit form uses
$t_{\nu_i,j}^{\rm diff}=6(\tau_j^{\nu_i})^2/\kappa_{t,j}^{\nu_i}$
in geometric units. Small optical depth recovers the free rate; large optical
depth gives diffusion-limited leakage. This is a local interpolation between
limits, not explicit neutrino propagation.

The loss time is defined through its inverse: for number leakage,
$1/t_{\nu_i,0}^{\rm loss}=R_{\nu_i}^{\rm free}/n_{\nu_i}$; for energy
leakage, $1/t_{\nu_i,1}^{\rm loss}=Q_{\nu_i}^{\rm free}/e_{\nu_i}$, where
$n_{\nu_i}$ and $e_{\nu_i}$ are the equilibrium neutrino number and energy
densities. Thus rate divided by density is an inverse loss time, not a loss
time.

### Current Heavy-Lepton Exception

The checked-in source does not apply the ideal same-species energy-loss ratio
to `nux`. Its `nux` effective-energy numerator is the free heavy-lepton rate,
and its suppression factor uses `nux` energy optical depth, `nux`
energy-transport opacity, and `nux` equilibrium energy density. However, the
inverse loss time in that factor uses the free `nue` energy rate:

$$
\left(t_{\nu_x,1}^{\rm loss}\right)^{-1}_{\rm current}
=\frac{Q_{\nu_e}^{\rm free}}{e_{\nu_x}},
$$

not $Q_{\nu_x}^{\rm free}/e_{\nu_x}$. Both the matter-source and luminosity
routines have this exception. Describe it as current observable behavior;
source provenance alone does not establish whether it is intentional.

The luminosity fixture generator exercises nonzero optical depths drawn from
1 through 1000, including this suppression branch. Its replay calls
`ghl_pert_test_fail` for all three luminosities but discards the returned
booleans, so a numerical mismatch cannot fail that executable. The
optically-thin matter-source evolution constructs all optical depths as zero
and therefore does not exercise diffusion suppression.

The electron-neutrino degeneracy is also interpolated between a transparent
value of zero and equilibrium using `exp(-tau->nue[0])`; the antineutrino uses
its own number optical depth, and heavy-lepton degeneracies are zero.

## Matter-Source Signs And Caller Meaning

Before unit conversion, current source convention is

$$
R_{\rm source}=R_{\bar\nu_e}^{\rm eff}-R_{\nu_e}^{\rm eff},
$$

$$
Q_{\rm source}=-\left(Q_{\nu_e}^{\rm eff}
+Q_{\bar\nu_e}^{\rm eff}+4Q_{\nu_x}^{\rm eff}\right).
$$

Thus positive `R_source` raises the matter electron fraction and negative
`R_source` lowers it. Positive emission/cooling rates produce negative
`Q_source`, reducing matter internal energy. The factor four sums the four
assumed-equal heavy-lepton species.

The generated C routine multiplies the number source by the atomic mass unit
and its cgs-to-geometric conversion; it applies the separate energy-source
conversion to `Q_source`. In the repo-local homogeneous, optically thin test
the caller uses

$$
\frac{dY_e}{dt}=\frac{R_{\rm source}}{\rho},\qquad
\frac{d\epsilon}{dt}=\frac{Q_{\rm source}}{\rho}.
$$

These equations are direct evidence for that test's static homogeneous-fluid
usage. A GRMHD evolution must still apply its own conservative-variable,
metric, lapse, and volume factors; do not treat the test RHS as a complete
general-relativistic update rule.

Ground truth:

- [`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
- [`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`](../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c)

## Units And Porting Boundary

Current GRHayL leakage constants assume geometrized units with
$G=c=M_\odot=1$. The public routines expect:

- `rho` in GRHayL geometric density units; generated source converts it to
  g/cm$^3$ internally;
- `T` and all chemical potentials in MeV;
- `Y_e`, `X_n`, `X_p`, and optical depths dimensionless;
- returned opacities in inverse GRHayL geometric length; and
- returned matter sources in GRHayL geometric units, consistent with the test
  divisions by `rho` above.

The fixed conversion constants live in `ghl_nrpyleakage.h`. Copying only the C
files into a code with cgs units or a different geometric mass scale requires
replacing or adapting those conversions. The HDF5 dependency comes from
GRHayL's EOS call path, not from the analytic leakage formulas themselves.

NRPyLeakage owns a separate geometric-unit conversion set in
`ghl_nrpyleakage.h`. Its macros are not numerically identical to the
`CODE_TO_CGS_*` set in `ghl_nrpyeos_tabulated.h`; replace or validate the
leakage density, length, number-rate, and energy-rate conversions coherently
rather than borrowing EOS conversion macros piecemeal.

### Current Suppression-Ratio Unit Seam

In the generated source-term and luminosity effective-rate expressions, the
current code multiplies the code-time diffusion factor
$6\tau^2/\kappa_{\rm geom}$ by an inverse loss rate obtained as a cgs rate
divided by a cgs density, hence expressed in $s^{-1}$. It does not convert that
inverse loss rate to inverse code time before forming the nominally
dimensionless suppression ratio. Holding the other quantities fixed, the
implemented ratio is larger than the unit-coherent ratio by
`NRPyLeakage_units_cgs_to_geom_T`, approximately $2.03\times10^5$.

The source does not state whether exact historical equivalence or a
unit-coherent physical ratio is intended. A port seeking exact GRHayL behavior
must preserve the current factor. A port seeking the unit-coherent formulation
must make and validate an explicit correction, for example by converting the
cgs inverse loss rate with `NRPyLeakage_units_geom_to_cgs_T`, or by expressing
both diffusion and loss times consistently in cgs. This is tracked as a
maintainer decision in [Current Contradictions](../../contradictions.md); no
single repair is prescribed here.

## Assumptions And Limitations

- Grey, energy-averaged rates: no evolved neutrino spectrum or angular
  distribution.
- Local leakage source: no nonlocal absorption/heating or neutrino momentum
  deposition in ejecta.
- Equal treatment of four heavy-lepton species, represented by one `nux`
  value.
- Free processes limited to beta capture, pair annihilation, transverse
  plasmon decay, and nucleon-nucleon bremsstrahlung.
- Opacity limited to the absorption/scattering channels described above.
- Positive finite temperature and physically consistent EOS composition are
  caller/EOS obligations; the leakage functions do not independently validate
  those physics preconditions.
- `EnsureFinite` fallbacks in generated expressions are numerical guards, not
  a substitute for a convention-correct EOS or a physical validity check.
- Path-of-least-resistance update considers six face neighbors. Full optical-
  depth initialization/convergence and transparent outer-boundary treatment
  remain caller or framework responsibilities.

## Ground Truth References

Repo-local implementation authority:

- [`GRHayL/Neutrinos/NRPyLeakage/`](../../../GRHayL/Neutrinos/NRPyLeakage/)
- [`GRHayL/include/ghl_nrpyleakage.h`](../../../GRHayL/include/ghl_nrpyleakage.h)
- [`GRHayL/include/ghl_radiation.h`](../../../GRHayL/include/ghl_radiation.h)
- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h)
- [`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`](../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c)
- [`Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`](../../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c)
- [`Unit_Tests/unit_test_nrpyleakage_luminosities.c`](../../../Unit_Tests/unit_test_nrpyleakage_luminosities.c)

Original generator, derivation, and validation material:

- https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD
- https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Formulation.ipynb
- https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Leakage_Scheme-Implementation.ipynb
- https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Numerical_implementation.ipynb
- https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Optically_thin_semi_analytic_calculation.ipynb

Physics and table-convention references:

- https://awsteiner.org/code/o2scl/html/class/eos_sn_oo.html
- https://arxiv.org/abs/2208.14487
- https://doi.org/10.1103/PhysRevD.107.044037
- https://ui.adsabs.harvard.edu/abs/1996A%26A...311..532R/abstract
- https://stellarcollapse.org/equationofstate.html
