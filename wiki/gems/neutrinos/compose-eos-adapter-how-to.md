# CompOSE EOS Adapter How-To

## Purpose And Authority

This page describes how to supply the six EOS values required by GRHayL's
NRPyLeakage routines from a CompOSE equation of state. It covers table
selection, state and unit mapping, chemical-potential normalization,
free-nucleon fractions, callback behavior, and validation. It does not define
a new leakage model or a generic CompOSE reader.

Use these authorities in this order:

1. Current GRHayL source and headers define the callback ABI, units at the
   leakage boundary, failure behavior, and quantities consumed by checked-in
   NRPyLeakage C.
2. The CompOSE reference manual and selected table's data sheet define the
   CompOSE files, axes, fields, particle content, and thermodynamic
   conventions.
3. The O'Connor-Ott EOSmaker source defines the common-reference convention
   for EOSmaker-produced StellarCollapse target tables.

The required GRHayL callback is
`ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`. Its output-pointer
order is

```text
muhat, mu_e, mu_p, mu_n, X_n, X_p
```

as declared in
[`ghl_eos_functions_declaration.h`](../../../GRHayL/include/ghl_eos_functions_declaration.h).
See [Physics And EOS Contract](physics-and-eos-contract.md) for how current
leakage source consumes those values.

## Choose An Integration Route

| Route | Use when | Work required | Main constraint |
| --- | --- | --- | --- |
| Convert offline to StellarCollapse HDF5 | Existing GRHayL tabulated-EOS reader should remain unchanged | Generate a compatible regular grid, convert units and chemical-potential zero points, derive or aggregate fields that CompOSE does not provide in the target convention, and write every dataset required by `NRPyEOS_stellarcollapse_read_table` | Current reader requires the complete StellarCollapse dataset family, not only the six leakage fields |
| Add a native CompOSE-backed GRHayL adapter | GRHayL should own table loading and expose the existing leakage callback | Add table storage/loading, state conversion, interpolation, callback registration, and GRHayL error mapping | Current three leakage entry points retain an HDF5 compile guard even though the physical callback need not use HDF5 |
| Port NRPyLeakage into another host | Host already owns GRHD, EOS, table interpolation, and errors | Implement the equivalent six-output EOS adapter and translate GRHayL types, errors, units, and the HDF5 guard | Exact equivalence requires preserving current source formulas and all boundary conventions, including the documented suppression-ratio seam |

The first route gives closest reuse of current GRHayL. The second avoids
forcing CompOSE data into an unrelated file schema. The third is the narrowest
choice for a host that does not otherwise use GRHayL. All routes target the
same six meanings. Interpolation, range behavior, and numerical results can
differ and require validation. Here “exact equivalence” means reproducing
current observable GRHayL behavior, including existing convention seams; it
does not mean silently making those seams physically unit-coherent.

CompOSE's generated `eoscompose.h5` is not loadable by the current GRHayL
StellarCollapse reader: the schemas differ. The offline route must write the
complete schema owned by
[`NRPyEOS_stellarcollapse.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c),
not rename six CompOSE fields.

That full schema contains 19 three-dimensional quantities, three grid arrays,
three grid-size scalars, and `energy_shift`. Beyond the six leakage outputs, a
converter must handle pressure, entropy, sound speed, internal energy and its
shift, thermodynamic derivatives, the target adiabatic-index convention,
`munu`, and nuclear composition fields such as `Abar`, `Zbar`, `Xh`, and `Xa`.
Some are unit conversions; others require a target-format convention,
derivative, individual-species extraction, or table-specific nuclear-group
aggregation. If these cannot be constructed consistently, use a native
CompOSE-backed adapter instead of writing placeholders. Grouped nuclear
quantities required by the full file do not replace the individual
free-neutron and free-proton fractions required by leakage.

### Implemented Fixed Offline Profile

The repository now implements only one production instance of the offline
route:
[`sro-sly4-sna-141-regularized-v1`](../../../tools/compose/README.md).
[`compose_to_grhayl.py`](../../../tools/compose/compose_to_grhayl.py) consumes
the fixed official generated-HDF5 schema and controls for CompOSE table 141,
then emits the existing StellarCollapse schema. It adds no native CompOSE
runtime backend, mapping language, resampling, or public C API.

This profile is intentionally a physics-changing regularized surrogate. Strict
temperature-ray regression, entropy reconstruction, joint causal derivative
projection, chemical canonicalization, and constrained composition projection
replace invalid raw fields. Do not identify its output with official table 141
or use it as an oracle for unchanged CompOSE values. The embedded manifest and
tool report disclose the fixed policy and measured changes.

For current leakage, node composition closure and the six returned values are
the hard interface. Independent trilinear interpolation of `Abar`, `Zbar`, and
`Xh` does not preserve nonlinear charge closure away from nodes; this is a
recorded diagnostic, not a supported equilibrium-composition claim. Current
NRPyLeakage consumes `Xn` and `Xp`, not those three heavy-nucleus fields. The
fixed Python suite and
[`unit_test_tabulated_eos_compose.c`](../../../Unit_Tests/unit_test_tabulated_eos_compose.c)
qualify the serialization, unchanged runtime interpolation and inversion, and
the six-output callback order and range failures with runtime sound-speed
cleaning disabled.

## Prerequisites

Requirements depend on the selected route:

- Every route needs the selected table page and data sheet, its parameter and
  thermodynamic files, and individual neutron/proton composition data.
- A raw-file adapter needs a reader for the exact CompOSE file-format version.
- Generating output or using official `compose` as a reference needs a matching
  release of the official program. HDF5 output additionally needs that
  release's HDF5 writer, its documented build toolchain, and a compatible
  HDF5 library. Build switches and defaults can differ between the served
  manual and current code, so follow the instructions shipped with the chosen
  release rather than copying a switch from another version.
- The three current EOS-dependent GRHayL leakage entry points require an
  HDF5-enabled GRHayL build even if a new physical callback does not use HDF5.
  A direct port must instead resolve the equivalent compile-guard boundary
  explicitly. See [Build And CI](../../build-and-ci.md).

## Preflight The Selected Table

Do this before writing a converter or runtime callback. CompOSE properties are
table-specific even when two tables represent the same named nuclear model.

| Check | Required evidence | Stop condition |
| --- | --- | --- |
| Parameter space | Table page, data sheet, `eos.t`, `eos.nb`, and `eos.yq` cover the intended positive `T`, `n_b`, and `Y_q` states | Cold, beta-equilibrated, one-dimensional, or restricted table cannot represent evolved states |
| Thermodynamic metadata | First row of raw `eos.thermo` supplies model `m_n`, model `m_p`, and lepton flag `I_l` | Metadata missing or separated from generated data with no preserved provenance |
| Lepton model | `I_l`, particle list, and data sheet say whether electrons, muons, or neither are included and state any electron-muon chemical-potential relation | Electron chemical potential cannot be reconstructed from available fields |
| Composition | `eos.compo` exists and provides individual net neutron and proton fractions, standard indices `10` and `11`; antinucleons are negligible | Only grouped nuclei, total nucleon content, no free-nucleon entries, or appreciable antinucleons are present |
| Charge coordinate | Hydrodynamic composition can be mapped to CompOSE `Y_q` under the table's lepton and neutrality assumptions | Evolved `Y_e` does not determine `Y_q` |
| Units and source layer | Reader knows whether it consumes raw `Q3`-`Q5` or generated physical outputs | Raw scaling would be applied twice, or required neutron-mass shift cannot be restored |
| Interpolation | Exact grid arrays, axis ordering, range, and interpolation policy are known | Code would extrapolate or assume grid spacing not present in the downloaded table |

Keep the selected CompOSE table page and data sheet with the integration
documentation. Meanings of extra thermodynamic quantities, phase indices, and
groups of nuclei are defined per table.

## Acquire And Inspect CompOSE Data

CompOSE documents three supported workflows:

1. Download original ASCII files and use a project-owned reader.
2. Download the original files and use the official `compose` routines to
   read, interpolate, transform, and optionally write HDF5.
3. Use authenticated web tools to generate and visualize a customized table.

The raw parameter grids are `eos.t`, `eos.nb`, and `eos.yq`. Thermodynamic
data are in required `eos.thermo`; `eos.compo` and `eos.micro` are optional.
The official program generates `eos.quantities` and `eos.parameters`, then
produces `eos.table` or `eoscompose.h5`.

For generated output, task 1 is a required selection step, not merely a file
generation detail. Select thermodynamic identifiers `3`, `4`, and `5` and
composition-pair identifiers `10` and `11` for the six leakage outputs. Select
every additional quantity needed by the chosen route; the offline
StellarCollapse route still has to construct its complete target schema.
Unselected identifiers do not appear in generated output.

No public REST API contract is documented by CompOSE. Build durable tooling
around the published files or official code, not guessed website endpoints.

## Map The Runtime State

CompOSE's general-purpose coordinates are

$$
(T, n_b, Y_q),
$$

where $T$ is in MeV, $n_b$ is in ${\rm fm}^{-3}$, and $Y_q$ is dimensionless.
The GRHayL callback receives `(rho, Y_e, T)`, with `rho` in GRHayL geometric
units and `T` in MeV.

### Density

CompOSE defines $n_b=\sum_i B_i n_i$. It does not mandate a mass-density
reference. Define the density reference explicitly:

$$
\rho_b=m_{\rm ref}^{\rho} n_b.
$$

Treat $m_{\rm ref}^{\rho}$ as explicit adapter metadata. Convert the callback's
geometric density to cgs with the EOS-side conversion, then convert number
density from $\mathrm{cm}^{-3}$ to $\mathrm{fm}^{-3}$:

$$
n_b[\mathrm{fm}^{-3}]
=\frac{\rho_b[\mathrm{g\,cm}^{-3}]}{m_{\rm ref}^{\rho}[\mathrm{g}]}10^{-39}.
$$

Do not infer $m_{\rm ref}^{\rho}$ from CompOSE's table-model neutron mass. Use
that mass only when the target hydrodynamic or adapter convention explicitly
selects it; otherwise preserve the independently declared density reference.
Current GRHayL also has two distinct numerical conversion sets. The live EOS
table loader converts imported cgs density to code units with
`CGS_TO_CODE_DENSITY`; an adapter mapping callback density back to cgs should
use the mathematical inverse of that live conversion, not infer behavior from
the source-unused `CODE_TO_CGS_DENSITY` definition. Leakage C independently
uses `NRPyLeakage_units_geom_to_cgs_D`, which is not numerically identical to
that inverse. An exact-current route preserves and tests both observed
conversions. A
unit-coherent port may replace them only as part of a deliberate, end-to-end
density, opacity, rate, and source normalization change. See the
[units and porting boundary](physics-and-eos-contract.md#units-and-porting-boundary).

### Charge Fraction

CompOSE defines $Y_q$ from the charge of strongly interacting particles. It
equals the net electron fraction only for a locally neutral model whose only
charged massive lepton is the electron:

$$
Y_q=Y_e.
$$

With electrons and muons under local charge neutrality,

$$
Y_q=Y_e+Y_\mu.
$$

For a baryonic-only table, $Y_q$ remains defined while $Y_e$ is not part of
that EOS. It is valid to query such a table with `Y_q = Y_e` only when the host
adds electrons separately and explicitly assumes electron-only local charge
neutrality. Do not alias the coordinates when muons or another charged
component are present unless the host supplies the missing fraction and
closure.

A muon-bearing state with $Y_q\ne Y_e$ is not a drop-in current NRPyLeakage
case even when an EOS can reconstruct `mu_e`. The callback receives only
`(rho, Y_e, T)`, and checked-in leakage C independently uses $Y_p=Y_e$ and
$Y_n=1-Y_e$ in reaction and scattering factors. Supporting such states needs
an interface and leakage-formulation extension that preserves both $Y_e$ and
$Y_q$, followed by end-to-end revalidation.

## Decode CompOSE Thermodynamics

Choose exactly one input layer. Raw CompOSE storage and generated CompOSE
output use different normalizations.

### Raw `eos.thermo`

The first row supplies the selected model's neutron mass $m_n$, proton mass
$m_p$ in MeV, and lepton flag $I_l$. The raw required fields include

$$
Q_3=\frac{\mu_b}{m_n}-1,\qquad
Q_4=\frac{\mu_q}{m_n},\qquad
Q_5=\frac{\mu_l}{m_n}.
$$

Recover the full CompOSE conserved-charge potentials before applying the
StellarCollapse mapping:

$$
\mu_b=m_n(1+Q_3),\qquad
\mu_q=m_nQ_4,\qquad
\mu_l=m_nQ_5.
$$

CompOSE chemical potentials are relativistic and include particle rest-mass
energy after this reconstruction.

### Generated `eos.table` Or `eoscompose.h5`

In ASCII `eos.table`, each row begins with `T`, `n_b`, and `Y_q`, followed by
selected quantities in the order encoded by `eos.quantities`; parse against
that selection, not fixed column numbers. In the current `eoscompose.h5`
writer, axes live at `Parameters/nb`, `Parameters/t`, and `Parameters/yq`.
Inspect that group's `pointsnb`, `pointst`, `pointsyq`, and
`tabulation_scheme` attributes from the matching writer release.

CompOSE thermodynamic selection identifier `3` is $\mu_b-m_n$ in MeV,
identifier `4` is $\mu_q$ in MeV, and identifier `5` is $\mu_l$ in MeV.
These identifiers are not fixed HDF5 array offsets. In `eoscompose.h5`, locate
each identifier through `Thermo_qty/index_thermo`, then read the corresponding
slot in `Thermo_qty/thermo`. The current official writer creates the
`Thermo_qty` group only when at least one thermodynamic quantity was selected;
its quantity count is the group attribute `pointsqty`. If `bshift` denotes the
value selected by identifier `3`, recover

$$
\mu_b={\tt bshift}+m_n.
$$

Do not read `thermo[...,3]` by assumption, multiply selected identifiers
`3`-`5` by $m_n$, or reapply the raw `Q` normalization. For composition,
locate particle identifiers `10` and `11` through
`Composition_pairs/index_yi`, then read the matching slots in
`Composition_pairs/yi`; do not assume `yi[...,10]` or `yi[...,11]`. The current
official writer creates this group only for a nonzero pair selection and stores
the count in its `pointspairs` attribute. Stop if either required group is
missing or its index dataset lacks a required identifier; regenerate with the
correct task-1 selection rather than substituting another field. These paths
come from the
[official HDF5 writer](https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/raw/master/hdf5compose.f90),
because the manual tables list dataset names but not the implemented group
hierarchy. Inspect output from the matching CompOSE release before coding
against it. Preserve raw-header $m_n$, $m_p$, and $I_l$ alongside generated
output: the CompOSE manual's documented generated-HDF5 dataset list does not
provide a replacement for all three raw-header values.

The same warning applies to energy. Raw
$Q_7=e/(n_bm_n)-1$ reconstructs total energy per baryon as
$E=m_n(1+Q_7)$, including rest mass. CompOSE explicitly warns that
$E/m_n-1$ is not hydrodynamic specific internal energy because $m_nn_b$ is
not the total mass density. Map energy using the host's declared density
reference and energy convention; do not relabel raw `Q7` as GRHayL `eps`.

## Build The Six Leakage Outputs

### Nucleon Chemical Potentials And `muhat`

CompOSE's particle-potential relation gives

$$
\mu_n^{\rm C}=\mu_b,\qquad
\mu_p^{\rm C}=\mu_b+\mu_q.
$$

For an EOSmaker-produced target table, the O'Connor-Ott convention restores
each particle's rest mass, then subtracts one common reference energy
$m_{\rm ref}^{E}$ from both nucleon potentials. To reproduce that target
convention, the callback outputs are

$$
\mu_n=\mu_b-m_{\rm ref}^{E},
$$

$$
\mu_p=\mu_b+\mu_q-m_{\rm ref}^{E},
$$

$$
\widehat\mu=\mu_n-\mu_p=-\mu_q.
$$

$m_{\rm ref}^{E}$ is explicit target-table or adapter metadata. The same value
must be subtracted from both potentials. The official Hempel SFHo EOSmaker
configuration sets $m_{\rm ref}^{E}$ to `amu_mev` and its density conversion
to `amu_cgs`; other EOSmaker readers can choose different references. Do not
infer a universal GRHayL nucleon zero point or assume that an independently
defined density reference $m_{\rm ref}^{\rho}$ equals $m_{\rm ref}^{E}$.
Subtracting $m_n$ from the neutron and $m_p$ from the proton would change their
difference and is not the Hempel EOSmaker target contract.

For that target convention, the direct formulas from raw CompOSE fields are

$$
\mu_n=m_n(1+Q_3)-m_{\rm ref}^{E},
$$

$$
\mu_p=m_n(1+Q_3+Q_4)-m_{\rm ref}^{E},
$$

$$
\widehat\mu=-m_nQ_4.
$$

For generated values `bshift`, `q`, and `l` selected by CompOSE identifiers
`3`, `4`, and `5`, the equivalent nucleon mapping is

$$
\mu_n={\tt bshift}+m_n-m_{\rm ref}^{E},
$$

$$
\mu_p={\tt bshift}+m_n+q-m_{\rm ref}^{E},
\qquad \widehat\mu=-q.
$$

If an intermediate source instead stores each potential with its own particle
rest mass removed, first restore that mass and then subtract the common
reference:

$$
(\mu_n^{\rm C}-m_n)+(m_n-m_{\rm ref}^{E})
=\mu_n^{\rm C}-m_{\rm ref}^{E},
$$

$$
(\mu_p^{\rm C}-m_p)+(m_p-m_{\rm ref}^{E})
=\mu_p^{\rm C}-m_{\rm ref}^{E}.
$$

Retain the selected table's $m_n$ and $m_p$ for these transformations. Do not
substitute generic particle masses.

Test both common zero-point failure modes. Separate neutron/proton rest-mass
subtraction shifts $\widehat\mu$ by approximately $-(m_n-m_p)$, about
$-1.293$ MeV. Choosing the wrong common $m_{\rm ref}^{E}$ leaves
$\widehat\mu$ unchanged but shifts both individual potentials, corrupting the
`mu_n/T` and `mu_p/T` blocking factors.

### Electron Chemical Potential

CompOSE gives the relativistic electron potential

$$
\mu_e=\mu_{l_e}-\mu_q.
$$

Choose the applicable electron-chemical-potential reconstruction branch; the
charge-coordinate and current-leakage API gate above still applies:

- Electron-only table with local charge neutrality: the effective stored
  $\mu_l$ equals $\mu_{l_e}$, so return $\mu_e=\mu_l-\mu_q$. This is
  $m_n(Q_5-Q_4)$ for raw fields and `l-q` for generated values.
- Electron-muon table whose data sheet documents
  $\mu_{l_e}=\mu_{l_\mu}$: the same formula is valid.
- Electron-muon table with unequal lepton potentials: CompOSE's stored
  effective $\mu_l$ is a density-weighted average and cannot by itself recover
  $\mu_e$. Require additional table data or reject the table.
- No-lepton table, indicated by $I_l\ne1$: CompOSE sets $\mu_l$ to zero. Obtain
  $\mu_e$ from a consistent external electron-positron EOS at the same `T` and
  net electron density. Do not interpret zero $\mu_l$ as the physical electron
  potential. The external sector must share the simulation's thermodynamic
  closure; naive post-addition can change Coulomb corrections and phase
  behavior assumed by the selected table.

Return `mu_e` including electron rest-mass energy. Do not subtract
$m_{\rm ref}^{E}$ or electron rest mass from it. CompOSE tables never include
neutrinos, so a neutrino chemical potential is not an alternative input.

### Free Neutrons And Protons

CompOSE defines particle fraction $Y_i^{\rm C}=n_i/n_b$ and mass-number
fraction $X_i=B_iY_i^{\rm C}$. The superscript distinguishes these net
particle fractions from the current leakage routine's separate local
approximations named `Y_n` and `Y_p`. For individual free neutrons and protons,
$B_i=1$, hence

$$
X_n=Y_n^{\rm C},\qquad X_p=Y_p^{\rm C}.
$$

Use net-particle composition indices `10` for neutrons and `11` for protons.
These must represent individual unbound nucleons. Do not substitute a grouped
nuclear fraction, average nucleus, total neutron content, `Y_q`, `Y_e`, or
the leakage routine's separate approximations `1-Y_e` and `Y_e`.

CompOSE fermion fractions are particle-minus-antiparticle net fractions. The
mapping above assumes antinucleons are absent or negligible. If they are
appreciable, net $Y_n^{\rm C}$ and $Y_p^{\rm C}$ are not the total
scattering-target abundances; the data mapping and current leakage physics
both require review.

Composition data are optional in CompOSE and zero or irrelevant species may
be omitted at individual grid points. A reader must map an omitted supported
particle to zero only when the file format and table metadata establish that
omission means zero. Absence of the composition file or absence of the
free-nucleon fields from the table is not permission to invent them.

## Implement The Callback Safely

Match this sequence:

1. Validate initialized EOS state, all six non-null output pointers, positive
   finite `rho` and `T`, and finite `Y_e`.
2. Convert geometric `rho` to cgs with the inverse of the live EOS-side
   `CGS_TO_CODE_DENSITY` conversion, then to CompOSE $n_b$ with
   $m_{\rm ref}^{\rho}$. Keep the separate leakage-side density conversion
   decision explicit.
3. Derive $Y_q$ from the host composition and table assumptions.
4. Reject states outside the selected table's `T`, $n_b$, or $Y_q$ ranges.
5. Interpolate the required thermodynamic and composition fields at one
   common state.
6. Decode raw or generated fields exactly once.
7. Apply the common-reference chemical-potential formulas and the selected
   lepton branch.
8. Validate finite outputs and physical composition bounds.
9. Commit all six output pointers only after every step succeeds.

Use temporaries so an error leaves all caller outputs unchanged, matching the
current GRHayL interpolation and leakage failure boundary. Return a
`ghl_error_codes_t`; use existing range errors where their meanings match.
Reject unsupported table structure during initialization when possible. If a
new runtime failure has no truthful existing code, add an explicit error
contract instead of returning an unrelated value.

When adding a native GRHayL backend, register the callback where tabulated EOS
function pointers are initialized; current registration lives in
[`NRPyEOS_initialize_tabulated_functions.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c).
The three EOS-dependent leakage C files return
`ghl_error_used_disabled_hdf5` before calling the callback in no-HDF5 builds.
Those three guards are only the first boundary. A non-HDF5 native backend must
also replace or restructure the current tabulated-EOS gates in
`ghl_initialize_eos_functions`, both tabulated initializers, function-pointer
registration, table read/free/interpolation paths, and `configure` source and
test filtering; GRHayLib separately requires HDF5. Audit every
`GHL_DISABLE_HDF5` site through the
[tabulated table contract](../eos/tabulated-table-contract.md#hdf5-build-gate)
and [Build And CI](../../build-and-ci.md). Editing the three leakage guards or
replacing only the callback is insufficient.

A native GRHayL backend also needs an explicit table-type/dispatch entry,
owned allocation and cleanup, loader and interpolation state, build-manifest
and HDF5 decisions, function-pointer initialization, truthful error mapping,
and direct tests. Adding only the six-value callback body leaves lifecycle and
dispatch unwired. Select an implemented table type through
`ghl_initialize_tabulated_eos`, which accepts `ghl_eos_table_t` directly. The
convenience wrapper `ghl_initialize_tabulated_eos_functions_and_params`
currently hard-codes `ghl_eos_table_stellarcollapse`; update that wrapper only
if the new backend must be selected through it. See the
[tabulated table contract](../eos/tabulated-table-contract.md).

## Grid And Interpolation Contract

The official `compose` program supports direct interpolation of individual
quantities with selectable orders and can generate a customized grid. Current
StellarCollapse grid requirements are owned by
[tabulated interpolation and bounds](../eos/tabulated-interpolation-and-bounds.md#current-stellarcollapse-grid-contract).
For CompOSE integration:

- an offline StellarCollapse conversion must generate at least two strictly
  increasing points on every axis: positive, uniformly log-spaced $n_b$ and
  $T$, plus uniformly linearly spaced $Y_q$;
- a native CompOSE adapter may retain CompOSE's exact grids and official
  interpolation policy, but must be validated against official `compose`
  output;
- neither route may silently extrapolate beyond the source table;
- $T=0$ is incompatible with current leakage ratios and GRHayL's logarithmic
  temperature interpolation, even when a CompOSE table contains that point;
- axis order and flattened array order must be verified with asymmetric test
  points, not inferred from equal-size grids.

Do not claim numerical equivalence between CompOSE's higher-order options and
GRHayL trilinear interpolation. Exact table nodes should agree; off-grid
differences must be measured against the selected reference policy.

## Validation Gates

Build a golden-state set covering low, middle, and high interior points plus
off-grid points. Include non-equal axis indices and avoid using only symmetric
states. Declare absolute/relative tolerances before comparison and compare only
states derived from the same underlying EOS, source table, and reference
convention.

Use an independently produced oracle only when it comes from the same raw table
release, CompOSE program version, quantity selection, interpolation policy,
and target density and chemical-potential conventions. A table from the same
named EOS family but a different source/version is useful as a harness input,
not as a numerical oracle for this mapping.

For each state:

1. At the same source grid node, decode raw `Q3`-`Q5`, then read generated
   identifiers through `index_thermo` and composition identifiers through
   `index_yi`. Compare $\mu_b-m_n$, $\mu_q$, $\mu_l$, $Y_n^{\rm C}$, and
   $Y_p^{\rm C}$ before interpolation. This catches double scaling and
   identifier/slot confusion.
2. Compare decoded values at exact nodes and off-grid states with official
   `compose` output under the selected interpolation policy.
3. Check
   $\mu_n^{\rm C}=\mu_b$,
   $\mu_p^{\rm C}=\mu_b+\mu_q$, and
   $\mu_n-\mu_p=-\mu_q$.
4. If an independently converted StellarCollapse oracle exists for the same
   raw table release, CompOSE version, grid and interpolation policy, and
   density and chemical-potential conventions, compare all callback outputs in
   ABI order and the leakage combinations `mu_e/T`, `mu_n/T`, `mu_p/T`,
   `muhat/T`, and `(mu_e-muhat)/T`.
5. Check $0\le X_n\le1$, $0\le X_p\le1$, and table-specific composition
   identities. Do not require $X_n+X_p=1$ where nuclei or other baryons exist.
6. Test electron-only and no-lepton closures independently. Treat muon-bearing
   states as unsupported until the charge-coordinate/API gate above is solved.
7. Confirm lower/upper range errors, zero temperature, nonfinite inputs/data,
   missing composition, unsupported lepton/charge mappings, and unchanged
   outputs on every error path.
8. Compare the current EOS-side and leakage-side density conversions
   separately. Do not infer their coherence from opacity agreement.
9. After callback validation, compare all six opacity slots, `R_source`,
   `Q_source`, and all three luminosities at zero and nonzero optical depth.
   Pressure or energy agreement alone does not test this interface. The
   opacity-only routine reads but does not consume `X_n` or `X_p`; their live
   bremsstrahlung use is in the combined source-term and luminosity routines,
   so validate those paths to exercise free-nucleon composition.

For a direct NRPyLeakage port, validate density, opacity, optical-depth,
source, and luminosity conversions together. Separately decide whether exact
compatibility requires preserving the current
[suppression-ratio unit seam](physics-and-eos-contract.md#current-suppression-ratio-unit-seam).

## Unsupported Without More Physics Or Data

Do not present these cases as mechanical conversions:

- cold or beta-equilibrated tables lacking independent positive `T` and $Y_q$
  coverage for evolved states;
- a muon-bearing state where evolved `Y_e` does not determine $Y_q$;
- unequal electron and muon lepton potentials when only effective $\mu_l$ is
  available;
- tables without individual free-neutron and free-proton composition;
- grouped nuclear fractions used in place of free nucleons;
- generated output whose model masses, lepton flag, or raw/generated field
  identity was not preserved;
- extrapolation outside any source axis;
- exotic, quark, or hybrid matter whose important opacity channels are not
  represented by current nucleon-only leakage formulas;
- changing either density or chemical-potential reference from the selected
  target convention without propagating and revalidating that change, or
  assuming the two independently defined references are equal without table
  evidence;
- treating raw $Q_7$ as hydrodynamic specific internal energy;
- treating a website endpoint as a stable CompOSE API without a published
  contract.

## Repo-Local Ground Truth

- Callback declaration:
  [`GRHayL/include/ghl_eos_functions_declaration.h`](../../../GRHayL/include/ghl_eos_functions_declaration.h)
- Callback registration:
  [`GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c)
- Current six-field interpolation:
  [`NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c`](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c)
- StellarCollapse reader and field names:
  [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c)
- StellarCollapse-to-GRHayL density conversion:
  [`NRPyEOS_stellarcollapse_to_ghl.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c)
- Current interpolation grid assumptions:
  [`NRPyEOS_tabulated_helpers.h`](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h)
- Fixed offline converter and policy:
  [`tools/compose/compose_to_grhayl.py`](../../../tools/compose/compose_to_grhayl.py) and
  [`tools/compose/README.md`](../../../tools/compose/README.md)
- Converter and runtime integration tests:
  [`Unit_Tests/compose/`](../../../Unit_Tests/compose/) and
  [`Unit_Tests/unit_test_tabulated_eos_compose.c`](../../../Unit_Tests/unit_test_tabulated_eos_compose.c)
- Leakage constants and unit conversions:
  [`GRHayL/include/ghl_nrpyleakage.h`](../../../GRHayL/include/ghl_nrpyleakage.h)
- EOS and leakage behavior:
  [Physics And EOS Contract](physics-and-eos-contract.md) and
  [API And Data](api-and-data.md)

## Ground Truth References

- [CompOSE Reference Manual](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf) (Version 3.01, served under the `manual_v3.00.pdf` path): sections 3.1, 3.3-3.5, 4.1, 4.2, 7.3, 7.5, tables 3.2, 3.3, 7.1, 7.7, 7.8, and appendices A-C.
- [CompOSE Quick Guide For Users](https://compose.obspm.fr/download/pdf/CompOSE_Quick_Guide_for_Users.pdf)
- [CompOSE table catalog](https://compose.obspm.fr/table/)
- [CompOSE software page](https://compose.obspm.fr/software/)
- [Official CompOSE code repository](https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/tree/master)
- [Official CompOSE HDF5 writer](https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/raw/master/hdf5compose.f90)
- [O'Connor-Ott StellarCollapse EOS page](https://stellarcollapse.org/equationofstate.html)
- [O'Connor-Ott EOSmaker archive](https://stellarcollapse.org/EOS/EOSmaker_svn_r13.tar.gz): `src/fix_table_units.F90` defines the common-reference nucleon-potential mapping.
- [EOSmaker common-reference conversion](https://github.com/evanoconnor/EOSmaker/blob/master/src/fix_table_units.F90)
- [EOSmaker Hempel SFHo configuration](https://github.com/evanoconnor/EOSmaker/blob/master/src/specific_eos/eos_table_Hempel_sc_sfho.F90): selects `amu_mev` and `amu_cgs` for that converter.
