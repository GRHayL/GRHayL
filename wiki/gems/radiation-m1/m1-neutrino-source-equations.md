# M1 neutrino source equations

This leaf gives the current grey one-body neutrino source equations and the
local E/F-plus-number update order. It uses the provider’s frozen rate bundle;
it does not define weak-interaction microphysics. The relevant implementation
is [ghl_m1_neutrino_sources.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c)
and the public declarations are in
[ghl_m1.h](../../../GRHayL/include/ghl_m1.h#L1261).

## Frozen rate fields

For one species, the single-species source path consumes the one-body fields

$$
\eta_N,\ \eta_E,\ \kappa_{a,N},\ \kappa_{a,E},\ \kappa_s,
\qquad
\kappa_{\rm tr}=\kappa_{a,E}+\kappa_s.
$$

The provider also supplies comoving grey targets `n_eq` and `J_eq`, a positive
mean-energy diagnostic, and charged-current subset fields
\(\eta_{N,cc},\kappa_{a,N,cc}\). Radiation consumes these frozen values. It
does not perform EOS lookup, weak-equilibrium calculations, unit conversion, or
rate refresh inside the local source solve; see the
[rate-provider contract](rate-provider-contract.md) and
[provider header](../../../GRHayL/include/ghl_neutrino_rate_provider.h#L14).

`eta_N`, `eta_E`, and the absorption/scattering coefficients are nonnegative.
The validator enforces the provider’s aggregate and charged-current
identities. Separated pair/plasmon/bremsstrahlung fields are not valid inputs
to this single-species equation set; they require the paired operation described
in [PAIR_SOURCE_MODEL.md](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md).

## Aggregate energy/momentum source

The comoving energy relaxation residual is

$$
Q=\eta_E-\kappa_{a,E}J.
$$

The radiation-side four-force represented by the current E/F source is

$$
G^\mu=Q u^\mu-\kappa_{\rm tr}H^\mu.
$$

Projecting with the current conventions gives the undensitized conservative
source components

$$
S_E=QW+\kappa_{\rm tr}H_n,
$$

$$
S_i=QWV_i-\kappa_{\rm tr}H_i.
$$

Here `H_n = -V_i H^i`, and `H_i` is covariant. The source implementation
uses `V_cov[i]` and `HD[i]` exactly for these lower-index terms; see
[ghl_m1_neutrino_sources.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c#L25).
Scattering contributes to \(\kappa_{\rm tr}\) and damps the comoving flux;
it does not appear in the number reaction below.

## Number source

After deriving the endpoint number-current normalization,

$$
S_N=\eta_N-\kappa_{a,N}\frac{N}{\Gamma_N}.
$$

This is an undensitized local source. In an explicit conservative RHS it is
multiplied by \(\alpha\sqrt{\gamma}\), just like the E/F interaction
projections. The number source uses the current `Gamma_N`, not the fluid `W`.

## Endpoint backward-Euler number update

The public number update uses `dt_alpha = alpha * dt` and an endpoint current
normalization:

$$
N^{n+1}=\frac{N^{\rm base}+\Delta t_\alpha\eta_N}
 {1+\Delta t_\alpha\kappa_{a,N}/\Gamma_N^{n+1}}.
$$

`N_base` is the source base. In the normal dispatcher this is the
transport-predicted state, not the pre-transport `state_input`. In the default
implicit path, the endpoint \(\Gamma_N^{n+1}\) is derived from the converged
E/F endpoint before complete-state repair. The number update uses that
provisional endpoint normalization; the candidate is then repaired and
checked. Current is derived again from the published repaired endpoint for
bounds, diagnostics, and exchange, without repeating the number update. The
denominator is validated as finite and positive. Scattering-only number
evolution leaves `N` unchanged because \(\kappa_{a,N}=0\) and the corresponding
emissivity identity is enforced by the rate bundle.

The formula does not apply `N_floor` internally. Number-floor repair is a
separate explicit operation, and an undesired floor mutation is rejected by
the conserving paired source path. See
[ghl_m1_neutrino_sources.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c#L237)
and [number-current definitions](m1-number-current-and-transport.md).

## Implicit E/F residual

The local homogeneous implicit solve treats the densitized E/F vector as its
four unknowns:

$$
\mathbf U=(\widetilde E,\widetilde F_x,\widetilde F_y,\widetilde F_z)^T,
\qquad
\widetilde E=\sqrt{\gamma}E,\quad
\widetilde F_i=\sqrt{\gamma}F_i.
$$

For a substep base \(\mathbf U_{\rm base}\), the residual is

$$
R_0(\mathbf U)=\widetilde E-\widetilde E_{\rm base}
-\Delta t_\alpha\sqrt{\gamma}\,S_E(\mathbf U),
$$

$$
R_i(\mathbf U)=\widetilde F_i-\widetilde F_{i,\rm base}
-\Delta t_\alpha\sqrt{\gamma}\,S_i(\mathbf U).
$$

The trial E/F state is undensitized by `sqrt_detgamma` before closure and
source evaluation. Metric, fluid primitives, and all rate fields are frozen;
the trial closure, comoving moments, and E/F source are reevaluated. `N` and
the number source are intentionally absent from this four-variable residual.
Repair is not applied inside the residual; admissibility is enforced by trial
checks and safeguarded Newton steps. See
[ghl_m1_neutrino_implicit_residual.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit_residual.c#L6).

The Jacobian is a deterministic 4x4 finite difference of this residual. A
forward perturbation is used normally; a backward one-sided perturbation is
allowed when the forward trial fails only the admissibility check. Other hard
errors propagate. This is documented in
[ghl_m1_neutrino_implicit_jacobian.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit_jacobian.c#L5).

## Local update order and compatibility branches

The established default source policy is the local implicit route:

1. solve E/F with the frozen-rate residual;
2. derive the current from the converged E/F endpoint before complete-state
   repair;
3. apply the endpoint number backward-Euler formula;
4. repair/check the complete endpoint;
5. rederive the published endpoint current and assemble diagnostics/exchange.

The opt-in branched compatibility policy can select an explicit thin E/F step
when `dt_alpha*kappa_a_E < 1` and `dt_alpha*kappa_s < 1`. Its number update is
still the endpoint-Gamma backward-Euler formula. The thick and
scattering-dominated compatibility predictors use a comoving backward-Euler
relaxation,

$$
J^*=\frac{J+\Delta\tau\eta_E}{1+\Delta\tau\kappa_{a,E}},
\qquad
H_i^*=\frac{H_i}{1+\Delta\tau\kappa_{\rm tr}},
\qquad
\Delta\tau=\frac{\Delta t_\alpha}{W},
$$

then use the \(\chi=1/3\) isotropic boost predictor before the same endpoint
number/current/exchange checks. These are internal source branches, not
alternative transport or closure methods. The dispatch and transaction rules
are in [ghl_m1_neutrino_source_update.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c#L612).

## Photon and obsolete-method boundary

The projection structure is shared M1 mathematics. The photon whitepaper’s
substitution `J_eq = a_R T^4` and photon opacity formulas are not neutrino
rules; current neutrino equilibrium targets and emissivities come from the
provider. The older whitepaper’s HLL path and reduced number-current formula
must not be substituted for the current Rusanov/full-current implementation.

## Focused evidence

- [Neutrino source-update tests](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
  covers endpoint number updates, source branches, pair boundaries, diagnostics,
  and transactional failures.
- [Finite-difference Jacobian tests](../../../Unit_Tests/unit_test_m1_fd_jacobian.c)
  covers residual/Jacobian construction, admissibility, Newton safeguarding,
  and rollback paths.
- [Seeded source/invariant tests](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  covers instantaneous source projections and conservation-oriented invariants.
