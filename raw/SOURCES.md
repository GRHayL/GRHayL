# Source Manifest

> Pointer manifest for GRHayL KB source material. Code, Doxygen source (`docs/raw/`), configuration, tests, and build inputs stay in place. Status is `living` for mutable sources, `frozen` for version-pinned specifications, and `frozen-historical` for superseded snapshots.

## Aggregate Sources

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `grhayl-core-and-headers` | Top-level library manifest, `GRHayL_Core` shared routines, and public headers under `GRHayL/make.code.defn`, `GRHayL/GRHayL_Core/**`, and `GRHayL/include/**`. | living | ingested |
| `grhayl-atmosphere-con2prim-eos` | Atmosphere, Con2Prim (including C++ Cerda-Duran and NN guess), and hybrid/tabulated EOS sources under `GRHayL/Atmosphere/**`, `GRHayL/Con2Prim/**`, and `GRHayL/EOS/**`. | living | ingested |
| `grhayl-flux-induction-reconstruction` | Flux/source, vector-potential induction, and PLM/PPM/WENO-z reconstruction sources under `GRHayL/Flux_Source/**`, `GRHayL/Induction/**`, and `GRHayL/Reconstruction/**`. | living | ingested |
| `grhayl-neutrinos-and-compose` | NRPyLeakage radiation sources under `GRHayL/Neutrinos/**` and offline CompOSE EOS conversion scripts under `tools/compose/**`. | living | ingested |
| `grhaylib-cactus-implementation` | Direct-compile Cactus/Einstein Toolkit thorn CCL files, symlink layout, and initialization sources under `implementations/GRHayLib/**`. | living | ingested |
| `unit-tests-and-sample-tables` | C unit-test drivers, `data_gen/` reference generators, Python CompOSE tests, and sample HDF5 tables under `Unit_Tests/**`. | living | ingested |
| `doxygen-raw-docs` | Hand-authored Doxygen pages, GRMHD derivation, stylesheets, and BibTeX inputs under `docs/raw/**`. | living | ingested |
| `scripts-and-ci-automation` | Manifest parser helpers under `scripts/**` and GitHub Actions workflows, composite actions, and test runner under `.github/**`. | living | ingested |

## Cited Code And Config Sources

Add an exact local path here only when no aggregate row above covers it.

| Source | Status |
| --- | --- |
| `AGENTS.md` | living |
| `README.md` | living |
| `AUTHORS` | living |
| `LICENSE` | living |
| `THIRD_PARTY_NOTICES` | living |
| `configure` | living |
| `generate_makefile.sh` | living |
| `Doxyfile` | living |
| `DoxygenLayout.xml` | living |
| `codecov.yml` | living |
| `.clang-format` | living |
| `.editorconfig` | living |
| `.gitattributes` | living |
| `.gitignore` | living |
| `tools/kb_lint.py` | living |

## External Sources

| Source | Tier | Provenance | Status | Ingest | Notes |
| --- | --- | --- | --- | --- | --- |
| `https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md` | background | Andrej Karpathy `LLM Wiki` gist note. | frozen | partial | Source, compiled-page, and format-rule separation. |
| `https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#onpushpull_requestpull_request_targetpathspaths-ignore` | external-spec | GitHub Actions workflow syntax reference for `paths-ignore`. | living | partial | Ignored-only push/PR filtering behavior for workflows. |
| `https://raw.githubusercontent.com/GRHayL/TestData/main/EOS/simple_table.h5` | external-spec | GRHayL/TestData sample tabulated EOS HDF5 file. | living | partial | Downloaded fallback table for HDF5 unit tests. |
| `https://stellarcollapse.org/EOS/` | background | StellarCollapse tabulated EOS repository index. | living | partial | Upstream catalog of nuclear EOS HDF5 tables. |
| `https://stellarcollapse.org/EOS/EOSmaker_svn_r13.tar.gz` | external-spec | StellarCollapse EOSmaker SVN r13 tarball. | frozen | partial | Legacy StellarCollapse HDF5 table conversion reference. |
| `https://stellarcollapse.org/EOS/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2` | external-spec | StellarCollapse LS220 tabulated EOS archive. | frozen | partial | Sample LS220 HDF5 table referenced in test setup. |
| `https://stellarcollapse.org/EOS/SLy4_3335_rho391_temp163_ye66.h5.bz2` | external-spec | StellarCollapse SLy4 tabulated EOS archive. | frozen | partial | Sample SLy4 HDF5 table used by NRPyLeakage tests. |
| `https://stellarcollapse.org/SROEOS.html` | background | Schneider-Roberts-Ott SRO EOS overview page. | living | partial | SRO Skyrme nuclear EOS table provenance and naming. |
| `https://stellarcollapse.org/equationofstate.html` | external-spec | StellarCollapse tabulated EOS format documentation. | living | partial | HDF5 dataset names, units, and chemical-potential signs. |
| `https://compose.obspm.fr/download/3D/SRO/SLy4/eos.pdf` | external-spec | CompOSE SRO SLy4 table datasheet PDF. | frozen | partial | SRO SLy4 table grid bounds and chemical-potential conventions. |
| `https://compose.obspm.fr/download/pdf/CompOSE_Quick_Guide_for_Users.pdf` | background | CompOSE Quick Guide for Users PDF. | living | partial | Table download and `compose` workflow overview. |
| `https://compose.obspm.fr/download/pdf/manual_v3.00.pdf` | external-spec | CompOSE Manual v3.00 PDF. | frozen | partial | CompOSE ASCII/HDF5 columns, `Q_i`, `mu_b`/`mu_q`/`mu_l`, and units. |
| `https://compose.obspm.fr/eos/141` | background | CompOSE EOS table 141 (SRO SLy4) catalog page. | living | partial | Catalog entry for SRO SLy4 nuclear EOS table. |
| `https://compose.obspm.fr/software/` | background | CompOSE official software downloads page. | living | partial | Official CompOSE Fortran utilities and `hdf5compose`. |
| `https://compose.obspm.fr/table/` | background | CompOSE EOS table browser. | living | partial | Online CompOSE table selection and metadata portal. |
| `https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/raw/master/hdf5compose.f90` | external-spec | Upstream CompOSE `hdf5compose.f90` HDF5 writer. | living | partial | Exact CompOSE HDF5 dataset layout and variable indices. |
| `https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/tree/master` | background | Upstream CompOSE Fortran source repository. | living | partial | Reference CompOSE interpolation and table generation tools. |
| `https://github.com/evanoconnor/EOSmaker/blob/master/src/fix_table_units.F90` | external-spec | EOSmaker unit-conversion routine `fix_table_units.F90`. | living | partial | StellarCollapse energy-shift and unit-conversion conventions. |
| `https://github.com/evanoconnor/EOSmaker/blob/master/src/specific_eos/eos_table_Hempel_sc_sfho.F90` | external-spec | EOSmaker Hempel SFHo table adapter source. | living | partial | Reference CompOSE/Hempel-to-StellarCollapse chemical potentials. |
| `https://awsteiner.org/code/o2scl/html/class/eos_sn_oo.html` | background | O2scl `eos_sn_oo` class documentation. | living | partial | Supernova EOS chemical-potential and LS `munu` sign notes. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD` | background | Upstream `Tabulated_EOS_IllinoisGRMHD` repository. | living | partial | Ancestral NRPyLeakage standalone code and notebooks. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Leakage_Scheme-Implementation.ipynb` | external-spec | Ancestral NRPyLeakage implementation notebook. | living | partial | Historical NRPy+ generator notebook and `tmp_*` CSE provenance. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Formulation.ipynb` | external-spec | Ancestral NRPyLeakage formulation notebook. | living | partial | Historical derivation of emission/absorption rates and opacities. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Numerical_implementation.ipynb` | background | Annotated walkthrough of the older Fortran/C++ `ZelmaniLeak` thorn. | living | partial | Comparative `ZelmaniLeak` walkthrough; not the NRPy+ generator. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Optically_thin_semi_analytic_calculation.ipynb` | background | Ancestral NRPyLeakage optically thin validation notebook. | living | partial | Semi-analytic optically thin leakage comparison provenance. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/NRPyLeakageET/src` | background | Ancestral `NRPyLeakageET` Einstein Toolkit source tree. | living | partial | Historical ET thorn layout prior to GRHayL integration. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/nrpy_core` | background | Ancestral `nrpy_core` helper tree. | living | partial | Historical NRPy 1 code-generation helpers for NRPyLeakage. |
| `https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/standalone/NRPyLeakage` | background | Ancestral standalone `NRPyLeakage` C tree. | living | partial | Historical standalone C layout prior to GRHayL port. |
| `https://github.com/scott-maddox/fdint` | background | `fdint` Fermi-Dirac integral package repository. | living | partial | Rational-approximation Fermi-Dirac integral origin. |
| `https://github.com/scott-maddox/fdint/blob/master/fdint/_fdint.pyx` | external-spec | `fdint` Cython implementation `_fdint.pyx`. | living | partial | Reference Fukushima Fermi-Dirac polynomial coefficients. |
| `https://github.com/RelNucAs/bns_nurates` | background | `bns_nurates` neutrino microphysics repository. | living | partial | Dense-matter nucleon blocking comparison reference. |
| `https://github.com/GRHayL/GRHayLET/blob/main/IllinoisGRMHD/schedule.ccl` | background | Upstream `GRHayLET` `IllinoisGRMHD/schedule.ccl`. | living | partial | Downstream Cactus schedule context for GRHayLib callers. |
| `https://github.com/GRHayL/GRHayLET/blob/main/IllinoisGRMHD/src/backward_compatible_initialize.c` | background | Upstream `GRHayLET` initialization helper. | living | partial | Downstream parameter compatibility context for GRHayLib. |
| `https://github.com/GRHayL/GRHayLET/blob/main/IllinoisGRMHD/src/evaluate_phitilde_and_A_gauge_rhs.c` | background | Upstream `GRHayLET` Lorenz gauge RHS C driver. | living | partial | Downstream staggered `A_i`/`phitilde` gauge call pattern. |
| `https://github.com/GRHayL/GRHayLET/blob/scupp/GRHayLMHDX/IllinoisGRMHDX/src/evaluate_phitilde_and_A_gauge_rhs.cxx` | background | Upstream `GRHayLET` CarpetX Lorenz gauge RHS driver. | living | partial | Downstream CarpetX staggered gauge call pattern. |
| `https://github.com/GRHayL/GRHayLET/tree/main/GRHayLHD` | background | Upstream `GRHayLET` `GRHayLHD` thorn tree. | living | partial | Downstream Carpet hydro thorn context. |
| `https://github.com/GRHayL/GRHayLET/tree/main/GRHayLHDX` | background | Upstream `GRHayLET` `GRHayLHDX` thorn tree. | living | partial | Downstream CarpetX hydro thorn context. |
| `https://bitbucket.org/zach_etienne/wvuthorns/src/master/IllinoisGRMHD/src/Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs.C` | background | Legacy `WVUThorns` `IllinoisGRMHD` Lorenz gauge RHS source. | living | partial | Historical staggered Lorenz gauge finite-difference convention. |
| `https://arxiv.org/abs/1711.00868` | background | Siegel and Metzger (2018) arXiv landing page. | frozen | partial | Neutrino leakage formulation background. |
| `https://arxiv.org/abs/2208.14487` | external-spec | Werneck et al. (2023) NRPyLeakage arXiv landing page. | frozen | partial | NRPyLeakage emission, absorption, and opacity formulation. |
| `https://arxiv.org/html/1808.00006v2` | external-spec | ILEAS (Ardevol-Pulpillo et al., 2019) HTML preprint. | frozen | partial | Detailed-balance nucleon blocking and chemical-potential conventions. |
| `https://doi.org/10.1016/j.amc.2015.03.009` | external-spec | Fukushima (2015) AMC half-odd Fermi-Dirac paper DOI. | frozen | partial | Piecewise rational approximations for half-odd Fermi-Dirac integrals. |
| `https://doi.org/10.1016/j.amc.2015.03.015` | external-spec | Fukushima (2015) AMC inverse Fermi-Dirac paper DOI. | frozen | partial | Piecewise rational approximations for inverse Fermi-Dirac integrals. |
| `https://doi.org/10.1093/mnras/stz613` | external-spec | Ardevol-Pulpillo et al. (2019) ILEAS MNRAS 485 DOI. | frozen | partial | Published ILEAS nucleon blocking and shifted-moment pairing reference. |
| `https://doi.org/10.1103/PhysRevD.107.044037` | external-spec | Werneck et al., Phys. Rev. D 107, 044037 (2023). | frozen | partial | Published NRPyLeakage equations and test verification. |
| `https://doi.org/10.1103/PhysRevD.111.063053` | external-spec | Chiesa et al., Phys. Rev. D 111, 063053 (2025). | frozen | partial | High-density nucleon blocking qualification (`bns_nurates`). |
| `https://ui.adsabs.harvard.edu/abs/1996A%26A...311..532R/abstract` | external-spec | Ruffert, Janka, and Schaefer (1996) A&A ADS record. | frozen | partial | Classical neutrino leakage rates and blocking factors. |
| `https://zenodo.org/records/14811397` | background | Zenodo record 14811397 (EOS/CompOSE table archive). | frozen | partial | External EOS dataset provenance cited by StellarCollapse adapter. |

## Exclusions

- `Makefile`, `build/`, `test/`, `lib/`, `docs/html/`, `docs/latex/`, and `nrpy/`.
- Compiled object files (`*.o`), executables, generated binary/ASCII reference outputs (`*.bin`, `*.asc`), and downloaded `*.h5` tables.
- Tool caches (`__pycache__/`, `.mypy_cache/`, `.pytest_cache/`), scratch logs, transcripts, audit/scan/token/latest reports, and uncommissioned planning or ranking files.
