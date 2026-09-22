# Flux Source Generated NRPy Boundary

## Routing Purpose

Use this page when deciding whether a Flux_Source file is a generator script,
NRPy support module, checked-in generated C kernel, or build-list entry. It
records repo-local evidence only and does not claim upstream generator
provenance beyond checked-in files.

## Generator Scripts

Flux_Source keeps Python generator/source scripts beside generated C kernels:

- [GRHayL/Flux_Source/GRHayL_rhs.py](../../../GRHayL/Flux_Source/GRHayL_rhs.py)
  imports local NRPy support and calls source-term C generation.
- [GRHayL/Flux_Source/GRMHD_equations_new_version.py](../../../GRHayL/Flux_Source/GRMHD_equations_new_version.py)
  contains symbolic GRMHD equation support used by Flux_Source generators.
- [GRHayL/Flux_Source/IGM_All_Source_Terms.py](../../../GRHayL/Flux_Source/IGM_All_Source_Terms.py)
  contains the `ghl_calculate_source_terms` generation path.
- [GRHayL/Flux_Source/IGM_All_fluxes.py](../../../GRHayL/Flux_Source/IGM_All_fluxes.py)
  contains the naming/output path for HLLE flux kernels across directions and
  variants.
- [GRHayL/Flux_Source/IGM_Characteristic_Speeds.py](../../../GRHayL/Flux_Source/IGM_Characteristic_Speeds.py)
  contains the naming/output path for characteristic-speed kernels.

## Command Status

Install the pinned symbolic dependency in a disposable environment, then
generate into a new or empty staging directory from the repository root:

```sh
python3 -m pip install -r GRHayL/Flux_Source/requirements.txt
python3 GRHayL/Flux_Source/GRHayL_rhs.py /tmp/grhayl-flux-stage
```

The driver refuses destinations inside `GRHayL/Flux_Source`, never removes an
existing tree, emits the complete family, and checks exact output-set equality
against the root and variant Flux_Source build manifests. It rejects other
SymPy versions: the pin stabilizes symbolic simplification and CSE ordering.
The source-term generator separately normalizes the printer's known
`pow(tmp_35, 1.0/2.0)` form to `sqrt(tmp_35)`.

| Family | Status | Evidence |
| --- | --- | --- |
| source terms | **supported** | Pinned staged regeneration emits the checked-in interface and the explicit normalization preserves its `sqrt` form. |
| characteristic speeds | **supported** | Pinned staged regeneration emits all three current `ghl_*` status-returning interfaces. |
| four HLLE families | **supported** | Pinned staged regeneration emits all direct variants with current interfaces. |

Reviewers must syntax-compile every staged C file and compare it with the
checked-in generated result before copying any output into the source tree.

## NRPy Support Modules

[GRHayL/Flux_Source/nrpy/](../../../GRHayL/Flux_Source/nrpy/) is the local NRPy
support tree imported by the Flux_Source Python scripts. It contains modules
for BSSN/ADM conversion, finite differences, code output, indexed expressions,
parameter handling, and related support. Treat this directory as source input
for local generation, not as generated documentation output.

## Checked-In C Kernels

The build list
[GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
compiles shared Flux_Source kernels and routes variant subdirectories:

- `ghl_calculate_source_terms.c`
- `ghl_calculate_characteristic_speed_dirn0.c`
- `ghl_calculate_characteristic_speed_dirn1.c`
- `ghl_calculate_characteristic_speed_dirn2.c`
- `hybrid/`
- `hybrid_entropy/`
- `tabulated/`
- `tabulated_entropy/`

Variant-local `make.code.defn` files compile the checked-in
`ghl_calculate_HLLE_fluxes_dirn*_<variant>.c` kernels. These C files are part
of the built source tree even when they are generated or derived from Python
scripts.

`configure --disable-hdf5` retains the tabulated and tabulated-entropy flux
sources in generated build targets. The Ubuntu-Clang `c2p-failure` no-HDF5
variant checks the public symbols; table-dependent flux tests and generators
remain excluded in that mode.

## Drift Rules

Generated C and Python sources can drift. Review them together when changing:

- characteristic-speed equations or direction handling;
- HLLE flux conservative outputs or entropy outputs;
- source terms, metric derivative use, or extrinsic curvature use;
- variable naming, include lists, or generated function signatures;
- build-list entries that decide which checked-in kernels compile.

Do not copy generated formulas into KB pages. Link the source files and route
readers to [docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox) and
[docs/raw/derivation.md](../../../docs/raw/derivation.md) for read-only
equation evidence.

Checked-in C plus active build configuration define current executable
behavior. Python files define generation intent. Review both together and use
the staging command rather than hand-editing generated C alone.

## Evidence Links

- [GRHayL/Flux_Source/GRHayL_rhs.py](../../../GRHayL/Flux_Source/GRHayL_rhs.py)
- [GRHayL/Flux_Source/GRMHD_equations_new_version.py](../../../GRHayL/Flux_Source/GRMHD_equations_new_version.py)
- [GRHayL/Flux_Source/IGM_All_Source_Terms.py](../../../GRHayL/Flux_Source/IGM_All_Source_Terms.py)
- [GRHayL/Flux_Source/IGM_All_fluxes.py](../../../GRHayL/Flux_Source/IGM_All_fluxes.py)
- [GRHayL/Flux_Source/IGM_Characteristic_Speeds.py](../../../GRHayL/Flux_Source/IGM_Characteristic_Speeds.py)
- [GRHayL/Flux_Source/requirements.txt](../../../GRHayL/Flux_Source/requirements.txt)
- [GRHayL/Flux_Source/nrpy/](../../../GRHayL/Flux_Source/nrpy/)
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
- [GRHayL/Flux_Source/hybrid/make.code.defn](../../../GRHayL/Flux_Source/hybrid/make.code.defn)
- [GRHayL/Flux_Source/hybrid_entropy/make.code.defn](../../../GRHayL/Flux_Source/hybrid_entropy/make.code.defn)
- [GRHayL/Flux_Source/tabulated/make.code.defn](../../../GRHayL/Flux_Source/tabulated/make.code.defn)
- [GRHayL/Flux_Source/tabulated_entropy/make.code.defn](../../../GRHayL/Flux_Source/tabulated_entropy/make.code.defn)
- [configure](../../../configure)
