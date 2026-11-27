# Nucleon Blocking And EOS Conventions

> Density-derived free-nucleon blocking in `NRPyLeakage_nucleon_blocking.h`, Kirchhoff reaction-energy pairing, dense-matter limitations, and BNS_NURATES/DD2 qualification. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

`GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h` reconstructs free-neutron and free-proton kinetic degeneracies $\eta_n, \eta_p$ from $(\rho_{\rm cgs}, T, X_n, X_p)$ using a common bare nucleon mass and evaluates scattering and charged-current transition overlaps with stable `expm1` algebra and analytic endpoint limits. This eliminates dependence on arbitrary chemical-potential energy zeros and removes the old equal-population pole without hot-path quadrature or root solves, while preserving the EOS thermodynamic $\widehat\mu$ in equilibrium degeneracies and the paired reaction shift $q = \widehat\mu - T(\eta_n - \eta_p)$. External qualification against BNS_NURATES establishes `2.0--4.0%` dilute emission agreement and exposes large dense-matter discrepancies on DD2 in a comparison that includes nucleon effective masses, mean fields, weak magnetism/recoil, and decay corrections, none of which GRHayL includes.

## Detail

### Installed Density-Derived Blocking Model

Every EOS-dependent NRPyLeakage entry point includes the source-private header `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h`. For each free species $N \in \{n, p\}$, it inverts the nonrelativistic ideal-gas density relation with common bare mass $m_N = 938.91872\,\mathrm{MeV}$ (following ILEAS Appendix B, Eqs. (69)-(71)):

$$n_N = N_A \rho X_N = C(T) F_{1/2}(\eta_N), \qquad C(T) = \frac{4\pi(2 m_N T)^{3/2}}{(hc)^3}, \qquad B_N = \frac{X_N}{1 + \frac{2}{3}\max(\eta_N, 0)}.$$

For ordered free fractions $X_h \ge X_l$ and $a = \eta_h - \eta_l \ge 0$, the same-energy transition overlaps are:

$$Y_{h \to l} = \frac{X_h - X_l}{-\mathrm{expm1}(-a)}, \qquad Y_{l \to h} = e^{-a} Y_{h \to l},$$

using `expm1` and Sterbenz exact subtraction for close populations, plus the analytic equal-population limit $X F'_{1/2}(\eta) / F_{1/2}(\eta)$. At an exactly-one-zero endpoint, the occupied-to-empty overlap tends to the occupied fraction, the reverse overlap tends to zero, and every overlap-times-shifted-beta-moment product tends to zero exponentially; the implementation applies that zero-product limit directly without forming the divergent degeneracy difference or shift (this is the continuous extension of the ideal-gas closure, not a general interacting-EOS claim). A both-zero state returns a neutral result. Evaluated two-species states enforce $0 \le B_N, Y_{N \to N'} \le X_N$.

Scalar $F_{-1/2}$ and inverse-$F_{1/2}$ minimax rational fits are ported from FDINT / Fukushima under BSD-3 attribution, avoiding runtime quadrature, root iteration, or extra EOS lookups. GRHayL's density normalization, overlap algebra, endpoint handling, and reaction pairing are local adaptations.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h::nrpyl_compute_population`; corroboration: `Unit_Tests/unit_test_nrpyleakage_physics.c::main`.

### Equilibrium Pairing And Why Raw EOS `mu_n`/`mu_p` Fail

The EOS thermodynamic difference $\widehat\mu = \mu_n - \mu_p$ still defines $\eta_{\nu_e}^{\rm eq} = (\mu_e - \widehat\mu) / T$, preserving the EOS equilibrium target. Charged-current emission and ordinary absorption share one reaction-energy shift:

$$q = \widehat\mu - T(\eta_n - \eta_p),$$

with opposite threshold orientations and algebraic shifted Fermi moments derived from a common spectral parent satisfying Kirchhoff detailed balance (ILEAS Appendix B and Appendix C, Eqs. (100)-(109), distinguishing ordinary and stimulated absorption). Production uses a grey reduction that neglects charged-lepton mass in the phase-space polynomial and evaluates final-state lepton blocking at a representative energy.

In a nonrelativistic quasiparticle EOS, $\mu_i^{\rm table} = \mu_i^{\rm full} - E_{\rm ref}$ and $\eta_i^{\rm kin} = (\mu_i^{\rm table} + E_{\rm ref} - m_i - U_i) / T$, where $E_{\rm ref}$ is a producer zero-point shift and $U_i$ is the interaction mean field. Because the six-output callback exposes neither $E_{\rm ref}$, $U_i$, nor effective masses (and the table's `energy_shift` is strictly for logarithmic internal-energy storage), absolute `mu_n / T` and `mu_p / T` cannot serve as kinetic degeneracies. The legacy implementation fed raw `mu_n`, `mu_p` to scattering and full `muhat / T` to a transition quotient, creating zero-point dependence and a pole near equal populations; subtracting a vacuum mass gap alone cannot fix that quotient because even a classical unequal-mass gas adds a temperature-dependent kinetic term.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h::nrpyl_compute_beta_emission_moments`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h::nrpyl_compute_beta_absorption_moments_with_fermi`.

### Free Populations And Dense-Matter Effective-Mass Limitation

Using free fractions `X_n`, `X_p` instead of total `1 - Y_e`, `Y_e` prevents counting bound nucleons as free targets. In legacy SLy4 at $\rho = 10^{14}\,\mathrm{g\,cm^{-3}}, T = 1\,\mathrm{MeV}$: at $Y_e = 0.1$, $X_n \approx 0.33682542, X_p \approx 3.1828 \times 10^{-9}$; at $Y_e = 0.5$, $X_n \approx 1.7645 \times 10^{-10}, X_p \approx 2.1357 \times 10^{-8}$. Regularized CompOSE SRO-141 converter fractions are surrogate inputs rather than microscopic truth: at $\rho = 10^{14}\,\mathrm{g\,cm^{-3}}, T = 1\,\mathrm{MeV}, Y_e = 0.1$, raw $X_p \approx 3.1734 \times 10^{-9}$ vs. regularized $3.1254 \times 10^{-7}$; at $Y_e = 0.5$, raw $(1.7662 \times 10^{-10}, 2.1389 \times 10^{-8})$ vs. regularized $(0, 3.8423 \times 10^{-6})$.

Common-bare-mass inversion does not reproduce interacting dense matter. At raw SRO-141 $\rho = 10^{14}\,\mathrm{g\,cm^{-3}}, T = 10\,\mathrm{MeV}, Y_e = 0.1$, correcting the chemical reference and subtracting microscopic mean fields yields kinetic degeneracies $(3.17720, -0.60888)$, matching effective-mass density inversion $(3.17579, -0.60776)$, whereas common bare mass shifts the minority transition overlap by `28.47%` (and `10.27%` occupancy shift at $Y_e = 0.5$).

### Independent BNS_NURATES Qualification And Runtime Cost

External qualification compared GRHayL's blocking evaluator against BNS_NURATES beta kernels and six published DD2 merger states (Chiesa et al., *Phys. Rev. D* **111**, 063053, 2025). Scratch drivers, DD2 tables, and modified BNS_NURATES quadrature are external and absent from this repository.

At the dilute SLy4 state (with BNS_NURATES mean-field, effective-mass, weak-magnetism, and decay corrections disabled), the four $\nu_e$ and $\bar\nu_e$ emission moments differed by `2.0--4.0%`, and a `0.5 s` optically thin evolution (RK4 refinement from 2,000 to 4,000 steps changed final $Y_e$ by `3.4e-12`) yielded:

| Quantity | GRHayL | BNS_NURATES | Difference |
| --- | ---: | ---: | ---: |
| `Ye` | `0.6081532252` | `0.6025608998` | `+0.0055923` abs |
| specific energy | `0.6585731569` | `0.6584636022` | `+0.0166%` |
| temperature [MeV] | `0.9951385806` | `0.9950932673` | `+0.00455%` |

With full BNS_NURATES DD2 physics enabled (mean fields, effective masses, weak magnetism/recoil, decay), candidate-minus-reference differences for `(number emission, energy emission, number opacity, energy opacity)` were:

| Point | `rho` [g cm^-3] | `T` [MeV] | `Ye` | electron neutrino [%] | electron antineutrino [%] |
| --- | ---: | ---: | ---: | --- | --- |
| A | `6.92e14` | `12.41` | `0.0716` | `+35482, +17680, +19254, +5668` | `+1010, +2535, +1013, +2532` |
| B | `9.95e13` | `16.89` | `0.0576` | `+13.23, +12.74, +31.93, +14.37` | `+102.57, +157.66, +100.49, +155.69` |
| C | `9.89e12` | `8.92` | `0.0569` | `+0.48, +0.65, +18.58, +1.49` | `+45.84, +54.94, +42.61, +53.40` |
| D | `1.01e12` | `6.70` | `0.1153` | `-0.33, -0.28, +12.14, +2.98` | `+32.19, +38.71, +30.48, +38.02` |
| E | `9.99e10` | `3.59` | `0.1350` | `+0.25, +0.15, +9.71, +2.90` | `+17.71, +20.77, +16.18, +20.15` |
| F | `1.02e10` | `2.14` | `0.1874` | `+0.89, +0.77, +6.62, +2.52` | `+11.54, +13.10, +10.02, +12.46` |

At Point A, the full-model comparison differs by factors of ~`58--356` ($\nu_e$) and `11--26` ($\bar\nu_e$). Common-bare-mass inversion gives `q = 117.10 MeV` vs. `20.22 MeV` from DD2 mean fields and effective masses (DD2 effective masses are about `280 MeV`); this ~`96.9 MeV` gap in $q$ is a likely contributor that the comparison does not isolate. Quadrature refinement from 96 to 128 nodes shifted reference channels by at most `6.4e-6` relative, so the disagreement is not reference quadrature error. Published M1 Eddington factors classify A-D as diffusion-like and E-F as transitioning without providing leakage optical depths, so optical-depth-aware evolution tests remain needed. Neither project policy, ILEAS, nor BNS_NURATES defines a formal per-kernel acceptance threshold; regenerated replay binaries serve as implementation-derived regression baselines for this approximate model.

Whole-BNS runtime overhead is unmeasured: an alternating 500-invocation whole-process test of the optically thin generator vs. an `LD_PRELOAD` stub (`3.315 s` vs. `3.205 s`, `3.43%`) included loading an ~`879 MB` EOS table and cannot be scaled by a `35%` leakage share. Replacing common bare mass or adding coherent nuclear scattering requires an extended EOS interface exposing effective masses, mean-field shifts, or bound nuclear composition.

Evidence (descriptive behavior): `Unit_Tests/unit_test_nrpyleakage_physics.c::main`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h::nrpyl_compute_population`.

## Sources

- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h) - `nrpyl_compute_population`, `nrpyl_compute_beta_emission_moments`, `nrpyl_compute_beta_absorption_moments_with_fermi`, FDINT BSD-3 provenance
- [Unit_Tests/unit_test_nrpyleakage_physics.c](../../../Unit_Tests/unit_test_nrpyleakage_physics.c) - table-free blocking, endpoint, and Kirchhoff reference checks
- [ILEAS (Ardevol-Pulpillo et al., arXiv:1808.00006v2)](https://arxiv.org/html/1808.00006v2) - Appendices B and C density-derived degeneracies and reaction pairing
- [BNS_NURATES Repository](https://github.com/RelNucAs/bns_nurates) - independent beta-process reference kernels
- [Chiesa et al., PRD 111, 063053 (2025)](https://doi.org/10.1103/PhysRevD.111.063053) - six DD2 BNS merger comparison states
- [CompOSE Reference Manual v3.00](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf) - sections 3.5 and 4.2.2 thermodynamic and microscopic fields
- [FDINT Repository](https://github.com/scott-maddox/fdint) - scalar Fermi-Dirac minimax implementation
- [FDINT _fdint.pyx](https://github.com/scott-maddox/fdint/blob/master/fdint/_fdint.pyx) - ported `fdm1h` and `ifd1h` rational fits
- [Fukushima (2015) Half-Odd Fermi Integrals](https://doi.org/10.1016/j.amc.2015.03.009) - $F_{-1/2}$ minimax fit
- [Fukushima (2015) Inverse Fermi Integrals](https://doi.org/10.1016/j.amc.2015.03.015) - inverse $F_{1/2}$ minimax fit

## See Also

- Depends on: [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md)
- Depends on: [StellarCollapse Table Adapter](../eos/stellarcollapse-table-adapter.md)
- Implements: [Implementation Flow](implementation-flow.md)
- Validated by: [Neutrinos Tests And Reference Data](tests-and-reference-data.md)
- See also: [Generator Provenance](generator-provenance.md)
