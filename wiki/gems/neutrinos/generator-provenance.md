# NRPyLeakage Generator Provenance

## Purpose

This page routes questions about where NRPyLeakage's generated C came from,
why its formula blocks contain `tmp_*` variables, and how ancestral generator
evidence relates to current GRHayL. For callable behavior, use the current
[implementation flow](implementation-flow.md), public headers, build manifest,
and tests. The external tutorial repository is provenance and derivation
evidence, not current GRHayL source authority.

External claims below describe the repository and saved notebook contents
reviewed for this page. The linked `master` branch can move. The durable
physics and generation facts needed by GRHayL readers are therefore summarized
locally; recheck upstream contents before attempting regeneration.

## Authority Boundary

The Python/NRPy+ generator is not present under
[`GRHayL/Neutrinos/`](../../../GRHayL/Neutrinos/). It lives in Leo Werneck's
`Tabulated_EOS_IllinoisGRMHD` repository, alongside generated standalone and
Einstein Toolkit implementations. That repository explains the symbolic
expressions and generation process behind GRHayL's code, but it does not define
GRHayL's present API, error handling, build inclusion, units header, or tests.

Use this authority order:

1. Current GRHayL C, headers, manifests, and tests are ground truth for code
   that builds and runs here.
2. The external implementation notebook is direct evidence for the symbolic
   construction and original generation flow.
3. The external formulation and validation notebooks are derivation and
   cross-check evidence. They do not prove current behavior.
4. The external checked-in C is ancestral comparison evidence. Do not replace
   current GRHayL files with it or assume rerunning the notebook reproduces
   current files byte-for-byte.

## Four Notebook Roles

### Formulation

`Tutorial-Neutrino_leakage-Formulation.ipynb` is the physics derivation. It
develops the leakage model from GRMHD source terms through:

- beta processes, electron-positron pair annihilation, transverse plasmon
  decay, and nucleon-nucleon bremsstrahlung;
- neutrino absorption and scattering opacities, optical depths, and diffusion
  times;
- free, diffusion-limited, and effective number/energy rates; and
- the final lepton-number source `R` and cooling source `Q`.

It also discusses neutrino pressure, which is broader than the five-file
current GRHayL NRPyLeakage build boundary. Use it to recover the meaning and
derivation of expressions, not to infer an additional GRHayL API.

### Symbolic Implementation And C Generation

`Tutorial-Leakage_Scheme-Implementation.ipynb` is the generator source. It:

- imports the vendored `nrpy_core` modules, SymPy, and Astropy constants;
- declares interaction switches, physical constants, and unit conversions;
- builds Fermi-Dirac approximations of orders zero through five as symbolic
  expressions;
- constructs free emission/cooling rates for the four processes, neutrino
  number and energy densities, scattering and absorption opacities, diffusion
  times, effective rates, and `R_source`/`Q_source`;
- obtains C expressions through `outputC.outputC`, writes C functions through
  `outputC.outCfunction`, and uses `optical_depth_helpers` for the
  path-of-least-resistance update; and
- emits generated C and a header; in Einstein Toolkit mode it places
  `make.code.defn` beside them under `src/` and writes the CCL files at the
  thorn root.

The notebook's geometric-unit construction uses
`M = M_sun`, `L = G M / c^2`, `time = L / c`, `D = M / L^3`,
`R = D / time`, and `Q = c^2 D / time`, with `Q` converted to
`MeV / (cm^3 s)`. Current numeric constants and conversions remain owned by
[`ghl_nrpyleakage.h`](../../../GRHayL/include/ghl_nrpyleakage.h); do not copy
numbers from the notebook into current code without a source-level review.

### ZelmaniLeak Numerical Walkthrough

`Tutorial-Neutrino_leakage-Numerical_implementation.ipynb` is an annotated
walkthrough of the older `ZelmaniLeak` thorn. It describes its Fortran/C++
drivers, leakage and optical-depth routines, interpolation helpers,
deleptonization profile support, neutrino pressure, scheduling, and parameter
checks. It is useful comparative evidence for numerical choices, but it is not
the NRPy+ generator and is not a description of GRHayL's public interface.

### Optically Thin Semi-Analytic Check

`Tutorial-Neutrino_leakage-Optically_thin_semi_analytic_calculation.ipynb` is
an independent validation path. It interpolates a tabulated EOS, evaluates the
optically thin reaction rates, forms `R` and `Q`, and evolves the resulting
homogeneous-fluid electron fraction and energy/temperature system. It helps
interpret the intent of
[`unit_test_nrpyleakage_optically_thin_gas.c`](../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c),
but it neither generates GRHayL C nor establishes current test pass/fail
behavior.

## Preserved Generation Flow

The generator notebook has a `generate_ET_thorn` switch. Its saved code derives
the module name and output directory from that switch:

- false: `standalone/NRPyLeakage/`;
- true: `standalone/NRPyLeakageET/src/` for generated C, its header, and
  `make.code.defn`; `param.ccl`, `schedule.ccl`, and `interface.ccl` are written
  one directory above `src/`.

At initialization, the notebook recursively removes and recreates its selected
`standalone/<module>` output directory. Run it only in a disposable clone and
review the selected mode before execution.

The reviewed notebook's final driver invokes generators in this order:

1. Fermi-Dirac integral function;
2. low-level combined source-term/opacity functions for NRPy and HARM constant
   sets;
3. the combined driver that selects a constant set;
4. a grid-level optical-depth function;
5. low-level opacity functions for both constant sets;
6. the opacity driver;
7. the generated header; and
8. in Einstein Toolkit mode, parameter, schedule, interface, and build files.

The repository contents reviewed for this page included two checked-in output
families:

- `standalone/NRPyLeakage/` contains the standalone header, Fermi function,
  combined source/opacity driver and NRPy/HARM kernels, opacity driver and
  NRPy/HARM kernels, and optical-depth routine.
- `NRPyLeakageET/src/` contains the integrated Einstein Toolkit thorn sources,
  including analogous generated kernels, path-of-least-resistance code, and
  luminosity code.

These checked-in directories matter because the reviewed notebook's saved
output path and checked-in outputs did not form a single turnkey regeneration
contract. The reviewed contents supplied no pinned Python environment or
lockfile for the four notebooks. Regeneration can therefore vary with SymPy,
Astropy, NRPy+, and notebook state.

## Why Generated C Uses `tmp_*`

`nrpy_core/outputC.py` enables common-subexpression elimination by default,
uses the variable prefix `tmp`, and calls SymPy `cse` with numbered symbols
beginning at `tmp_0`. The implementation notebook passes the complete output
set for a kernel to `outputC.outputC`, so repeated algebra is hoisted into
ordered `const` temporaries before final assignments.

The same defaults are visible repo-locally in
[`GRHayL/Flux_Source/nrpy/outputC.py`](../../../GRHayL/Flux_Source/nrpy/outputC.py).
That file is a sibling NRPy+ copy used by another GRHayL module, not the
identical ancestral `nrpy_core/outputC.py`; it corroborates the CSE defaults
without replacing the external generator as provenance.

The `tmp_*` names therefore have no physics meaning. They are compiler-style
common subexpressions. Their numbering and grouping depend on symbolic input,
output ordering, SymPy/NRPy+ processing, and generator version. Read the named
symbolic expressions in the implementation notebook or the named outputs at
the end of a current C block; do not document a stable meaning for an
individual `tmp_N`.

## Relationship To Current GRHayL

The current build manifest
[`GRHayL/Neutrinos/NRPyLeakage/make.code.defn`](../../../GRHayL/Neutrinos/NRPyLeakage/make.code.defn)
lists five adapted C files:

- [`NRPyLeakage_Fermi_Dirac_integrals.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c)
  retains the generated approximation expressions, but GRHayL returns
  `ghl_error_codes_t` and writes through an output pointer instead of exiting
  on an unsupported key.
- [`NRPyLeakage_compute_neutrino_opacities.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c)
  and
  [`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
  preserve recognizable `tmp_*` sequences from the external standalone NRPy
  constant kernels, around a GRHayL EOS/API/error-handling adaptation.
- [`NRPyLeakage_compute_neutrino_luminosities.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c)
  preserves a generated formula block visible in the external
  `NRPyLeakageET_compute_neutrino_luminosities_nrpy_constants.c`. No
  corresponding luminosity generator was found among the notebooks reviewed
  for this page, so that external C file is the surviving ancestral evidence
  inspected for this block.
- [`NRPyLeakage_optical_depths_PathOfLeastResistance.c`](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c)
  extracts the pointwise face-average, neighbor-candidate, and minimum-path
  calculation seen in the external Einstein Toolkit routine into a
  GRHayL struct-based call.

GRHayL adaptations include different public types and argument shapes,
tabulated-EOS dispatch, HDF5-disabled returns, propagated Fermi/EOS errors,
struct writeback, and finite-value handling. Current constants have one
active NRPy family in
[`ghl_nrpyleakage.h`](../../../GRHayL/include/ghl_nrpyleakage.h). The same
header also retains dormant `NRPyLeakage_ZL_*` macros, but current leakage C
does not reference them. By contrast, the reviewed upstream standalone output
contains paired NRPy/HARM kernel families and a runtime selector. GRHayL has no
current runtime NRPy/HARM selector. Review exact source before assuming an
ancestral formula or failure mode survived unchanged.

## Regeneration Rule

No current repo-local command regenerates all five GRHayL files. If a formula
must change:

1. edit or reconstruct the named symbolic expression using the external
   implementation notebook and vendored `nrpy_core` as provenance;
2. generate into a disposable directory;
3. compare formulas against the current GRHayL kernel;
4. reapply current GRHayL API, EOS, error, unit, and finite-handling contracts;
   and
5. verify the built sources and tests routed by
   [Implementation Flow](implementation-flow.md) and
   [Tests And Fixtures](tests-and-fixtures.md).

Do not treat generated C replacement as a mechanical copy operation.

## Ground Truth References

- Original tutorial and implementation repository:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD
- Physics formulation notebook:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Formulation.ipynb
- Symbolic implementation and C-generation notebook:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Leakage_Scheme-Implementation.ipynb
- ZelmaniLeak numerical implementation walkthrough:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Numerical_implementation.ipynb
- Optically thin semi-analytic notebook:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Optically_thin_semi_analytic_calculation.ipynb
- Checked-in standalone generated outputs:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/standalone/NRPyLeakage
- Checked-in Einstein Toolkit implementation and generated kernels:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/NRPyLeakageET/src
- NRPy+ support used by the generator, including `outputC.py` and
  `optical_depth_helpers.py`:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/nrpy_core
