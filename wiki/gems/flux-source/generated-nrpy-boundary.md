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

Run generation only in a disposable copy. Current repo-local probe classifies
entry points as follows:

| Family | Status | Evidence |
| --- | --- | --- |
| source terms | **supported (disposable probe verified)** | From `GRHayL/Flux_Source`, `python3 GRHayL_rhs.py` exits zero and writes `./ghl_calculate_source_terms.c`. It requires Python, SymPy, and the checked-in `nrpy/` tree. Running from repo root fails the script's relative `nrpy/` import path. |
| characteristic speeds | **unknown / no supported command** | `GRHayL_rhs.py` imports the module but comments out its generation call. Its callable still spells unprefixed parameter types such as `primitive_quantities`, unlike checked-in public signatures. |
| four HLLE families | **unknown / no supported command** | Main-loop directory creation and flux calls are commented out. The callable derives both filename and C symbol suffix from `Ccodesdir` and still spells unprefixed parameter types. Do not invent a command from these comments. |

Disposable source-term regeneration produced C that passed
`cc -std=c99 -Wall -Wextra -fsyntax-only` with the public include directory,
but did not match the checked-in file byte-for-byte: one expression used
`pow(tmp_35, 1.0/2.0)` where checked-in C uses `sqrt(tmp_35)`. This is a
generator-drift signal, not permission to overwrite checked-in source. No
characteristic-speed or HLLE regeneration was run.

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

Variant-local `make.code.defn` files compile the 12 checked-in
`ghl_calculate_HLLE_fluxes_dirn*_<variant>.c` kernels. These C files are part
of the built source tree even when they are generated or derived from Python
scripts.

`configure --disable-hdf5` filters the six tabulated and tabulated-entropy
sources from generated build targets. Their files and public declarations
remain present; source presence is not link evidence in that mode.

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
behavior. Python files define generation intent. Neither side alone proves
reproducibility. Until speed/flux commands are restored and validated, do not
hand-edit generated C as though it were primary symbolic authority and do not
promise full-family regeneration.

## Evidence Links

- [GRHayL/Flux_Source/GRHayL_rhs.py](../../../GRHayL/Flux_Source/GRHayL_rhs.py)
- [GRHayL/Flux_Source/GRMHD_equations_new_version.py](../../../GRHayL/Flux_Source/GRMHD_equations_new_version.py)
- [GRHayL/Flux_Source/IGM_All_Source_Terms.py](../../../GRHayL/Flux_Source/IGM_All_Source_Terms.py)
- [GRHayL/Flux_Source/IGM_All_fluxes.py](../../../GRHayL/Flux_Source/IGM_All_fluxes.py)
- [GRHayL/Flux_Source/IGM_Characteristic_Speeds.py](../../../GRHayL/Flux_Source/IGM_Characteristic_Speeds.py)
- [GRHayL/Flux_Source/nrpy/](../../../GRHayL/Flux_Source/nrpy/)
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
- [GRHayL/Flux_Source/hybrid/make.code.defn](../../../GRHayL/Flux_Source/hybrid/make.code.defn)
- [GRHayL/Flux_Source/hybrid_entropy/make.code.defn](../../../GRHayL/Flux_Source/hybrid_entropy/make.code.defn)
- [GRHayL/Flux_Source/tabulated/make.code.defn](../../../GRHayL/Flux_Source/tabulated/make.code.defn)
- [GRHayL/Flux_Source/tabulated_entropy/make.code.defn](../../../GRHayL/Flux_Source/tabulated_entropy/make.code.defn)
- [configure](../../../configure)
