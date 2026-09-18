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
- [`NRPyLeakage_nucleon_blocking.h`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h)
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
| `mu_n` | Neutron chemical potential | Returned by the callback for ABI compatibility; current leakage rates do not consume it directly |
| `mu_p` | Proton chemical potential | Returned by the callback for ABI compatibility; current leakage rates do not consume it directly |
| `muhat` | Neutron-proton chemical-potential difference; GRHayL's StellarCollapse adapter defines `muhat = mu_n - mu_p` | `(mu_e-muhat)/T` in the retained grey equilibrium-neutrino degeneracies and spectral moments |
| `X_n` | Free-neutron mass fraction | Neutron scattering population, charged-current transition overlap, and nucleon-nucleon bremsstrahlung composition factor |
| `X_p` | Free-proton mass fraction | Proton scattering population, charged-current transition overlap, and nucleon-nucleon bremsstrahlung composition factor |

`X_n` and `X_p` are free-nucleon mass fractions, not total neutron and proton
fractions. Current blocking uses these values so bound nucleons are not counted
as free scattering or capture targets. An EOS adapter must not substitute
`1-Y_e` and `Y_e`, particularly where nuclei or other bound species are
present.

The chemical potentials and temperature must share MeV units because current
source forms `mu_e/T` and `(mu_e-muhat)/T` without conversion. `Y_e`, `X_n`,
and `X_p` are dimensionless. Blocking does not use the absolute values of
`mu_n` or `mu_p`.

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

### Chemical Potentials And Density-Derived Blocking

An alternative EOS must reproduce the combinations consumed by the code, not
only similarly named outputs. In current source the equilibrium degeneracies
are

$$
\eta_{\nu_e}^{\rm eq}=\frac{\mu_e-\widehat\mu}{T},\qquad
\eta_{\bar\nu_e}^{\rm eq}=-\eta_{\nu_e}^{\rm eq},
$$

with `muhat` consumed directly as $\widehat\mu$. GRHayL's StellarCollapse
adapter reads the table's `muhat` field and defines it as `mu_n - mu_p`, and
the CompOSE converter writes the same difference. That stored difference
retains the neutron-proton rest-mass gap, so it is the correct argument here.

Blocking requires kinetic occupations. Full thermodynamic nucleon chemical
potentials also contain rest-energy, interaction, and producer-reference
contributions. The former implementation used `mu_n/T` and `mu_p/T` directly
for scattering and used full `muhat/T` in a same-spectrum transition quotient.
That made scattering change under a common chemical-energy-zero shift and made
the transition expression singular near equal populations. Rest-gap
subtraction alone would not remove interaction, effective-mass, or
available-volume ambiguity.

Current source instead reconstructs effective kinetic degeneracies from the
free populations. For species $N$,

$$
n_N=N_A\rho X_N=C(T)F_{1/2}(\eta_N),\qquad
C(T)=\frac{4\pi(2m_NT)^{3/2}}{(hc)^3},
$$

using the common bare mass $m_N=938.91872\,\mathrm{MeV}$. This choice follows
the density-derived free-gas blocking construction in ILEAS Appendix B,
Eqs. (69)--(71). It removes the arbitrary chemical-energy reference and uses
only the EOS population identified as free. It remains an ideal-gas
approximation: it does not recover EOS mean fields, effective masses,
available-volume corrections, or microscopic spectra.

The scattering populations are

$$
B_N=\frac{X_N}{1+\frac{2}{3}\max(\eta_N,0)}.
$$

For ordered free fractions $X_h\geq X_l$ and
$a=\eta_h-\eta_l\geq0$, the stable same-energy transition populations are

$$
Y_{h\rightarrow l}=\frac{X_h-X_l}{-\operatorname{expm1}(-a)},\qquad
Y_{l\rightarrow h}=e^{-a}Y_{h\rightarrow l}.
$$

The common normalization converts the difference of Fermi integrals into
`X_h-X_l`; `expm1` retains digits when the degeneracies are close. Equal
degeneracies use the analytic limit
$X F'_{1/2}(\eta)/F_{1/2}(\eta)$. As the lower population tends to zero,
the occupied-to-empty overlap tends to the occupied fraction and the reverse
overlap tends to zero. The degeneracy difference diverges, but each overlap-
times-shifted-beta-moment product tends to zero exponentially. At an exact
endpoint, the implementation applies that product limit without constructing
the divergent reaction shift; the both-zero state remains a neutral success.
This is the continuous extension of the installed density-derived closure,
not a general claim about interacting-EOS endpoint rates. The implementation
checks finite inputs and results and enforces
$0\leq B_N,Y_{N\rightarrow N'}\leq X_N$ for evaluated two-species states.
These identities avoid numerical quadrature, root iteration, and a new lookup
table in the leakage hot path.

The scalar $F_{-1/2}$ and inverse-$F_{1/2}$ rational fits come from Scott
Maddox's [FDINT implementation](https://github.com/scott-maddox/fdint/blob/master/fdint/_fdint.pyx),
which implements Fukushima's minimax approximations
([inverse integral](https://doi.org/10.1016/j.amc.2015.03.015),
[half-odd integral](https://doi.org/10.1016/j.amc.2015.03.009)). GRHayL's
density normalization and overlap identities are local adaptations; FDINT and
ILEAS do not supply that combined evaluator. The close-population subtraction
uses Sterbenz's exact-subtraction result (P. H. Sterbenz,
*Floating-Point Computation*, Prentice-Hall, 1974, Sec. 4.3).

Do not independently add or subtract $Q$ elsewhere, or reconstruct
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
`NRPyLeakage_ZL_Q_npmass`, and no current leakage C file uses either. Do not
insert a gap into `(mu_e-muhat)/T` or into the transition argument without a
paired reaction-model change. These constants are unrelated to
`Q_source`. In particular, `muhat(Ye=.5)` is not a rest-mass measurement:
even a classical unequal-mass gas has a temperature-dependent kinetic term.

The actual legacy SLy4 SNA/NSE and CompOSE SRO-141 energy conventions are
established in the [table adapter](../eos/stellarcollapse-table-adapter.md#producer-energy-conventions).
They remain relevant to `muhat` and equilibrium moments, but absolute `mu_n`
and `mu_p` no longer control blocking.

The overlap equations obey same-energy detailed balance within the selected
common-spectrum gas because
$Y_{l\rightarrow h}/Y_{h\rightarrow l}=e^{-a}$. Charged-current emission and
ordinary absorption use a shared reaction shift

$$
q=\widehat\mu-T(\eta_n-\eta_p),
$$

paired threshold orientation, and algebraic shifted Fermi moments. Their
spectral parent satisfies Kirchhoff balance with the EOS equilibrium chemical
potential. Production moments neglect charged-lepton mass in the phase-space
factor and use representative-energy final-state blocking. This follows the
algebraic shifted-moment structure in ILEAS Appendix B and the pairing in
Appendix C, Eqs. (100)--(109); it does not reproduce a finite-mass spectral
rate exactly. No physical clamp is applied to `q` or `q/T`: the current EOS
API exposes no authoritative mean-field bound. Large finite ratios therefore
remain a conditioning limitation. Nonfinite final results are repaired to
documented finite fallbacks and reported with
`ghl_error_nrpyleakage_nonfinite_output`; this status does not establish
accuracy for large finite ratios.

Normalized absorption can divide two representable subnormal Fermi moments.
Supported builds therefore require gradual underflow; flushing either operand
to zero changes a finite physical ratio before final-output sanitization can
detect the loss.

Independent qualification against the beta kernels in
[BNS_NURATES](https://github.com/RelNucAs/bns_nurates) found `2.0--4.0%`
differences in the four dilute emission moments and an absolute `Ye`
difference of `0.00559` after a `0.5 s` thin-gas evolution. At the published
DD2 merger points, agreement improves toward lower density, but the densest
point fails badly: common-bare-mass inversion gives `q=117.10 MeV`, while the
DD2 mean-field/effective-mass reference gives `20.22 MeV`; charged-current
rate and opacity differences range from factors of about `58` to `356` for
electron neutrinos and `11` to `26` for electron antineutrinos. The reported
128-node reference changed by at most `6.4e-6` relative from 96 nodes, so this
is physical-model disagreement rather than integration noise.

The density-derived model is numerically qualified. After reviewing the named
rate and thin-evolution comparisons, their uncertainties, and the dense DD2
failure, the owner accepted the results for this approximate leakage model and
authorized them as the golden baseline. This outcome-specific decision does not
define a transferable per-kernel tolerance or claim dense interacting-EOS
microscopic accuracy. EOS-consistent effective masses or mean-field shifts
remain an optional future accuracy improvement. See
[blocking qualification](tests-and-fixtures.md#blocking-correction-qualification).

Practical adapter validation should compare all six returned quantities,
`mu_e/T`, `(mu_e-muhat)/T`, and the derived blocking populations against known
states before comparing final rates. Raw `mu_n` and `mu_p` still belong to the
six-output ABI, but matching their ratios no longer validates active blocking.
Merely matching pressure and internal energy does not validate this interface.

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

The checked-in implementation applies this form to every species, including
`nux`; see [Heavy-Lepton Timescale](#heavy-lepton-timescale).

The implementation evaluates the ratio entirely in cgs, where
$t_{\nu_i,j}^{\rm diff}=6(\tau_j^{\nu_i})^2/(c\,\kappa_{t,j}^{\nu_i})$;
see [Suppression-Ratio Units](#suppression-ratio-units). Small optical depth
recovers the free rate; large optical depth gives diffusion-limited leakage.
This is a local interpolation between limits, not explicit neutrino
propagation.

The loss time is defined through its inverse: for number leakage,
$1/t_{\nu_i,0}^{\rm loss}=R_{\nu_i}^{\rm free}/n_{\nu_i}$; for energy
leakage, $1/t_{\nu_i,1}^{\rm loss}=Q_{\nu_i}^{\rm free}/e_{\nu_i}$, where
$n_{\nu_i}$ and $e_{\nu_i}$ are the equilibrium neutrino number and energy
densities. Thus rate divided by density is an inverse loss time, not a loss
time.

### Heavy-Lepton Timescale

Both the matter-source and luminosity routines compute the free energy rate of
one heavy-lepton species once, as the local `Q_free_nux`, and use it both as
the `nux` effective-energy numerator and in its own inverse loss time:

$$
\left(t_{\nu_x,1}^{\rm loss}\right)^{-1}
=\frac{Q_{\nu_x}^{\rm free}}{e_{\nu_x}}.
$$

The suppression factor uses `nux` energy optical depth, `nux`
energy-transport opacity, and `nux` equilibrium energy density. `lum->nux`
remains one heavy species; the matter-source cooling term multiplies the
single-species effective rate by four.

The luminosity fixture generator exercises nonzero optical depths drawn from
1 through 1000, including this suppression branch, and its replay fails on a
numerical mismatch. The optically-thin matter-source evolution constructs all
optical depths as zero and therefore does not exercise diffusion suppression.

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

### Suppression-Ratio Units

The generated source-term and luminosity effective-rate expressions build the
suppression ratio entirely in cgs. Their shared diffusion prefactor is
`6.0/NRPyLeakage_c_light`, so the diffusion time is
$6\tau^2/(c\,\kappa_{\rm cgs})$ and multiplies an inverse loss rate obtained as
a cgs rate divided by a cgs density, in $s^{-1}$. The ratio is therefore
dimensionless without any further conversion. Opacity is converted to
geometric inverse length only on writeback to `kappa`, after the ratio is
formed.

A port must not additionally apply `NRPyLeakage_units_geom_to_cgs_T` or
`NRPyLeakage_units_geom_to_cgs_L` to this prefactor: `6/L_unit` and
`6*T_unit/L_unit` differ by `NRPyLeakage_units_cgs_to_geom_T`, approximately
$2.03\times10^5$. Either express both times in cgs, as here, or convert both
consistently to code time.

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
