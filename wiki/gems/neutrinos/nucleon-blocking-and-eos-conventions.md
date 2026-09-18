# Nucleon Blocking And EOS Conventions

## Purpose

Nucleon blocking estimates how occupied neutron and proton states suppress
scattering and charged-current reactions. Three separate questions matter:
which nucleons participate, what kinetic spectrum represents them, and how
the reaction energy is paired between emission and absorption. The installed
model answers these questions with a fast ideal-gas approximation. Its density
inversion and population-overlap algebra are numerically stable and reference-
invariant; the reaction-energy shift can remain poorly conditioned at large
finite `|q/T|`. The owner accepted its measured qualification results for
leakage use while retaining the dense interacting-matter limitation.

See the [physics contract](physics-and-eos-contract.md) for species, units, and
equilibrium conventions, the [table adapter](../eos/stellarcollapse-table-adapter.md#producer-energy-conventions)
for producer-specific EOS meanings.

## Installed Blocking Model

The [private blocking helper](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h)
is shared by standalone opacities, opacities with GRMHD sources, and
luminosities. It uses the EOS callback's free-neutron and free-proton fractions
`Xn` and `Xp`. For each free population it inverts the nonrelativistic ideal-gas
density relation with a common bare nucleon mass to obtain a kinetic degeneracy
`eta_n` or `eta_p`.

For two occupations with the same kinetic spectrum, the transition overlap is

$$
Y_{np}=\frac{X_n-X_p}{1-\exp[-(\eta_n-\eta_p)]},\qquad
Y_{pn}=e^{-(\eta_n-\eta_p)}Y_{np}.
$$

The implementation evaluates this identity with stable `expm1` algebra and an
analytic equal-population limit. At an exactly-one-zero endpoint, the occupied-
to-empty overlap tends to the occupied fraction, the reverse overlap tends to
zero, and every overlap-times-shifted-beta-moment product tends to zero. The
code applies that analytic product limit directly instead of constructing the
divergent kinetic-degeneracy difference or reaction-energy shift. The both-zero
state returns a neutral result. For evaluated two-species states, overlaps are
projected to their initial free populations. This removes the former pole near
equal populations, keeps the result bounded, and makes blocking independent
of a common chemical-energy reference.

This endpoint is the continuous extension of GRHayL's installed density-
derived closure. It is not a claim that charged-current rates at a
single-species composition are generally zero for a richer interacting EOS.

The inverse half-order Fermi integral uses scalar minimax fits from FDINT and
Fukushima. No numerical quadrature, root iteration, or additional EOS lookup
runs in the leakage hot path. The source carries the FDINT license and cites
the original fits. GRHayL's density normalization, overlap identity, endpoint
handling, and reaction pairing are local adaptations.

## Equilibrium And Reaction-Energy Pairing

The EOS thermodynamic difference

$$
\widehat\mu=\mu_n-\mu_p
$$

still defines the electron-neutrino equilibrium degeneracy through
`(mu_e-muhat)/T`. It is not replaced by the auxiliary gas's kinetic chemical
difference. This preserves the EOS equilibrium target.

Charged-current emission and ordinary absorption use the common shift

$$
q=\widehat\mu-T(\eta_n-\eta_p).
$$

Their threshold orientations and algebraic shifted moments come from one
spectral parent, so the parent kernels satisfy Kirchhoff detailed balance.
Production uses a grey reduction: it neglects charged-lepton mass in the
phase-space polynomial and evaluates final-state lepton blocking at a
representative energy. This structure follows ILEAS Appendices B and C,
including its distinction between ordinary and stimulated absorption. It is
an approximation to finite-mass spectral rates, not an exact weak-interaction
calculation.

## Why EOS Chemical Potentials Cannot Directly Supply Blocking

For a nonrelativistic quasiparticle model,

$$
\mu_i^{\rm table}=\mu_i^{\rm full}-E_{\rm ref},\qquad
\eta_i^{\rm kin}=\frac{\mu_i^{\rm table}+E_{\rm ref}-m_i-U_i}{T}.
$$

The common reference `E_ref` is bookkeeping. The mean field `U_i` and effective
mass describe physical interactions. The existing callback does not return
those quantities. Absolute `mu_n/T` and `mu_p/T` therefore depend on an
unknown energy zero and cannot safely serve as kinetic degeneracies. The EOS
table's `energy_shift` instead supports logarithmic storage of specific
internal energy; it is not this chemical reference.

The earlier implementation fed absolute nucleon potentials to scattering and
used `muhat/T` in a population quotient. A common chemical-reference shift
changed scattering, and the quotient could be singular or unbounded near
equal populations. Merely subtracting a vacuum neutron-proton mass difference
does not repair that quotient: even an unequal-mass classical gas has a
temperature- and composition-dependent kinetic contribution to `muhat`.

## Free Populations And Producer Conventions

Blocking uses `Xn` and `Xp`, not total fractions `1-Ye` and `Ye`. A cold
legacy SLy4 state at `rho=1e14 g/cm^3`, `T=1 MeV`, and `Ye=.1` has approximately
`Xn=.33682542` and `Xp=3.1828e-9`; treating `.9` and `.1` as available targets
would count bound nucleons as free. At `Ye=.5`, the same table gives free
fractions near `1.7645e-10` and `2.1357e-8`, so equal total fractions need not
represent appreciable free populations.

Regularized converter fractions are still model inputs rather than recovered
microscopic truth. At the cold SRO-141 state `rho=1e14 g/cm^3`, `T=1 MeV`, and
`Ye=.1`, raw `Xp` is about `3.1734e-9`, while the regularized callback returns
about `3.1254e-7`. At `Ye=.5`, raw free fractions are approximately
`1.7662e-10,2.1389e-8`; the callback returns `0,3.8423e-6`. Closure repair can
therefore affect strongly blocked channels. Nonuniform matter may also require
an explicit available-volume convention.

## Effective-Mass Limitation

A common bare mass is sufficient to define the installed stable approximation,
but it does not reproduce interacting dense matter. A matched raw SRO-141 state
at `rho=1e14 g/cm^3`, `T=10 MeV`, and `Ye=.1` illustrates the missing input.
Correcting the table's chemical reference and subtracting the microscopic mean
fields gives degeneracies about `3.17720,-.60888`; independently inverting the
free densities with the table's effective masses gives `3.17579,-.60776`.
Using a common bare mass changes the minority transition overlap by about
`28.47%` in this nearly dissociated state. At `Ye=.5`, the corresponding
occupancy difference is about `10.27%`.

These occupancy differences do not directly equal leakage-rate or luminosity
errors. They show why reference invariance and bounded algebra cannot establish
dense-matter accuracy without effective masses or mean fields.

## Independent Qualification

Physical qualification used the beta kernels in
[BNS_NURATES](https://github.com/RelNucAs/bns_nurates) and the six published DD2
merger states from Chiesa et al., *Physical Review D* **111**, 063053 (2025),
[doi:10.1103/PhysRevD.111.063053](https://doi.org/10.1103/PhysRevD.111.063053).
The comparison separated two questions.

### Dilute matched model

At the existing dilute SLy4 state, with BNS_NURATES mean-field,
effective-mass, weak-magnetism, and decay corrections disabled, the four
electron-neutrino and electron-antineutrino emission moments differed by
`2.0--4.0%`. An independent `0.5 s` thin-gas evolution gave:

| Quantity | GRHayL candidate | BNS_NURATES | Difference |
| --- | ---: | ---: | ---: |
| `Ye` | `0.6081532252` | `0.6025608998` | `+0.0055923` absolute |
| specific energy | `0.6585731569` | `0.6584636022` | `+0.0166%` |
| temperature [MeV] | `0.9951385806` | `0.9950932673` | `+0.00455%` |

GRHayL RK4 refinement from 2,000 to 4,000 steps changed final `Ye` by
`3.4e-12`. This is positive dilute, optically thin evidence. No requester,
project, or cited source supplies a physical discrepancy limit that converts
it into formal application acceptance.

### Full DD2 model

The full comparison enabled BNS_NURATES DD2 mean fields, effective masses,
weak magnetism/recoil, and decay corrections. Candidate-minus-reference
differences for number emission, energy emission, number opacity, and energy
opacity were:

| Point | `rho` [g cm^-3] | `T` [MeV] | `Ye` | electron neutrino [%] | electron antineutrino [%] |
| --- | ---: | ---: | ---: | --- | --- |
| A | `6.92e14` | `12.41` | `0.0716` | `+35482, +17680, +19254, +5668` | `+1010, +2535, +1013, +2532` |
| B | `9.95e13` | `16.89` | `0.0576` | `+13.23, +12.74, +31.93, +14.37` | `+102.57, +157.66, +100.49, +155.69` |
| C | `9.89e12` | `8.92` | `0.0569` | `+0.48, +0.65, +18.58, +1.49` | `+45.84, +54.94, +42.61, +53.40` |
| D | `1.01e12` | `6.70` | `0.1153` | `-0.33, -0.28, +12.14, +2.98` | `+32.19, +38.71, +30.48, +38.02` |
| E | `9.99e10` | `3.59` | `0.1350` | `+0.25, +0.15, +9.71, +2.90` | `+17.71, +20.77, +16.18, +20.15` |
| F | `1.02e10` | `2.14` | `0.1874` | `+0.89, +0.77, +6.62, +2.52` | `+11.54, +13.10, +10.02, +12.46` |

At A, GRHayL's common-bare-mass inversion gives `q=117.10 MeV`; the DD2
mean-field/effective-mass reference gives `20.22 MeV`. DD2 effective masses are
about `280 MeV`, and the missing dense-matter information leaves about
`96.9 MeV` in the reconstructed shift.

The widened scratch-only BNS_NURATES integration changed the reported channels
by at most `6.4e-6` relative from 96 to 128 nodes. This numerical uncertainty
is far below the dense physical disagreement. Published M1 Eddington factors
classify A--D as diffusion-like and E--F as transitioning, but the six states
do not provide leakage optical depths. Diffusion can suppress sensitivity to
free emission while the opacity still controls diffusion and integrated
optical depth. An optical-depth-aware comparison is therefore needed before
deciding the effect on a full leakage evolution.

## Acceptance Status

- The evaluator's algebra, endpoint behavior, spectral-parent Kirchhoff
  pairing, and reported numerical comparisons pass their stated checks.
- Whole-BNS overhead is unmeasured. The earlier `3.43%` whole-process timing
  includes loading an approximately 879 MB EOS table and cannot be treated as
  leakage-local cost or scaled by the user-supplied `35%` leakage share. The
  production evaluator remains algebraic and adds no quadrature, root solve,
  or EOS lookup.
- The owner accepts the measured dilute rates, thin trajectory, and dense DD2
  comparison set for this approximate leakage model. Regenerated replay results
  are the authorized golden baseline.
- Dense DD2 rates and opacities show that the common-bare-mass model does not
  reproduce mean-field physics. This limits microscopic accuracy claims but
  does not block the accepted leakage approximation.
- The reaction shift is not physically clamped because the current EOS API
  supplies no authoritative mean-field bound. Large finite `|q/T|` remains a
  conditioning limitation; observable nonfinite-output errors do not remove
  that finite-value sensitivity.
- Scattering uses only the EOS free-neutron and free-proton fractions. Bound
  nucleons and coherent nuclear scattering are outside the accepted model.
- Replacing the common-bare-mass inversion or expanding the scattering targets
  requires an owner decision and a richer EOS/API contract. The current EOS
  boundary does not expose the effective masses, mean-field shifts, or bound
  nuclear composition needed to implement either change unambiguously.

An optional future accuracy change may expose EOS-consistent effective masses
or mean-field shifts at the leakage boundary, or define a separately justified
fallback for EOSs that lack them. An optical-depth-aware leakage comparison
would then assess the full evolution effect.
ILEAS's reported scheme-level transport agreement is not a per-kernel GRHayL
error allowance and is not reused as one.

## Ground Truth References

- [ILEAS](https://arxiv.org/html/1808.00006v2), especially Appendices B and C:
  density-derived gas degeneracies, algebraic shifted moments, and detailed
  balance pairing.
- [BNS_NURATES](https://github.com/RelNucAs/bns_nurates) and
  [Chiesa et al.](https://doi.org/10.1103/PhysRevD.111.063053): independent beta
  kernels and the DD2 BNS comparison states.
- [CompOSE manual, sections 3.5 and 4.2.2](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf):
  thermodynamic and microscopic field meanings.
- [FDINT](https://github.com/scott-maddox/fdint) and Fukushima's
  [forward](https://doi.org/10.1016/j.amc.2015.03.009) and
  [inverse](https://doi.org/10.1016/j.amc.2015.03.015) Fermi-integral fits:
  scalar evaluator provenance, not a nuclear-matter accuracy bound.
