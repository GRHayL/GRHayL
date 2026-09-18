# Grey electron-flavor pair collision model

`ghl_m1_solve_neutrino_pair_source_update` receives the electron neutrino and
antineutrino states together, in `{nue, anue}` order. It remains a pointwise library
operation: metric, matter primitives, and provider coefficients are frozen for
the stage. The host supplies the transport-predicted states and owns matter
recovery and final publication.

## Provider coefficients

Electron-flavor `eta_N`, `eta_E`, `kappa_a_N`, `kappa_a_E`, and `kappa_s`
describe the independent charged-current and scattering operator. Pair,
plasmon, and bremsstrahlung emission is supplied separately in
`eta_N_pair[c]` and `eta_E_pair[c]`. The number emissivity for a given process
must be the same for both members of the pair. The energy emissivities may
differ. Existing `n_eq` and `J_eq` are the comoving equilibrium targets.
Callers constructing bundles themselves must initialize the complete struct,
including zero process arrays when those channels are absent. Existing callers
using aggregate electron pair rates must migrate; this changes the public
rate struct's size and requires rebuilding consumers.

The production provider retains the raw number and energy emissivities. It
does not force the emitted mean energy `eta_E_pair/eta_N_pair` to equal the
equilibrium mean energy `J_eq/n_eq`. These are different spectral averages.
The table-free reference backend supplies synthetic, symmetric pair-number
emission and remains a test model, not a weak-interaction prescription.

The lumped heavy-flavor species continues to use its existing aggregate
approximation. Its state is already summed over the four heavy flavors and
has zero electron-lepton weight; callers must not apply another multiplicity.

## Shared number reaction

For each process let `q_c` be its common number emissivity. At the beginning
of a pair substep, compute each species' number-current normalization
`Gamma_s` from its E/F state. During that substep use `n_s = N_s/Gamma_s`.
The model's inverse number reaction is proportional to both occupancies:

```
R_c = q_c * (1 - n_e*n_a/(n_eq_e*n_eq_a)).
```

With `h = alpha*dt_sub`, backward Euler reduces to a scalar quadratic for a
shared undensitized Eulerian number increment `d`:

```
d = h * sum_c(q_c) *
    (1 - (N_e+d)*(N_a+d)/(Gamma_e*Gamma_a*n_eq_e*n_eq_a)).
```

Use the admissible root continuous with `d=0` as `h` tends to zero. Both
species receive that same `d`; neither number is independently clipped.
Consequently the pair operator preserves `N_e-N_a` up to arithmetic
roundoff. If either configured number floor would be violated, the candidate
is rejected. Pair processes produce no matter electron-fraction increment.

Numerically, solve the equivalent quadratic for the smaller surviving
population and recover the larger one from the original number difference.
This avoids subtracting a nearly complete annihilation increment from a much
larger initial population.

## Energy and momentum

After solving the number reaction, freeze the partner's updated comoving
number density within the E/F solve. The energy-weighted absorption for
species `s`, with partner `p`, is

```
K_s = sum_c(eta_E_pair_s[c]/J_eq_s) * (n_p_new/n_eq_p).
Q_s = sum_c(eta_E_pair_s[c]) - K_s*J_s.
S_E_s = W*Q_s + K_s*Hn_s.
S_i_s = W*Q_s*V_i - K_s*HD_i_s.
```

These use the same comoving moments and Eulerian source projections as the
existing M1 collision operator. Both species' E/F candidates use the same
pair substep schedule. At comoving equilibrium (`n=n_eq`, `J=J_eq`, `H=0`),
the number and energy/momentum sources vanish.

This is an isotropic grey emission and absorption approximation with separate
number and energy weighting. It freezes the number-current normalizations
and partner occupancies during each respective subsolve and is first-order
in time. It is not an energy-resolved or angle-resolved annihilation kernel.
The need for joint neutrino/antineutrino information and the additional
spectral/angular modelling in annihilation are discussed in
[Foucart's moment-transport formulation](https://academic.oup.com/mnras/article/475/3/4186/4810555).
The occupancy-product and time-splitting equations above specify this
library's approximation; they are not a claim of exact equivalence to that
paper's transport scheme.

## Opacity supplied to transport

The host-prepared four-point face operation still receives a total transport
opacity. For an electron species, the provider's scalar `kappa_tr` now contains
only the independent absorption and scattering contribution. Add the pair
inverse-energy contribution when preparing the opacity:

```
kappa_transport_s = rates_s.kappa_tr
                  + sum_c(eta_E_pair_s[c]/J_eq_s) * (n_partner/n_eq_partner).
```

Evaluate the partner density from the transport-stage state using
`n_partner = N_partner/(W-Hn_partner/J_partner)` and the same closure/comoving
moments as the number current. The host owns the existing cell-to-face opacity
construction. Passing scalar `kappa_tr` alone when pair channels are active
would omit their absorption from the thick-limit flux suppression. The lumped
heavy-flavor scalar already contains its aggregate pair contribution. Supply
the total as the face-operation operand; do not overwrite `kappa_tr` in the
validated rate bundle.

## Composition, exchange, and failure

The canonical paired update first advances independent charged-current and
scattering sources into temporary states, then advances the pair operator.
Charged-current lepton exchange is accumulated from the first operation.
Energy and momentum exchange are computed from the final published radiation
increments and have exactly opposite densitized matter counterparts.

`n_b_cons` is the undensitized Eulerian baryon number density, `W*rho/m_b`,
compatible with the undensitized charged-current number increment. Thus
`dYe = -dL_rad_cc/n_b_cons`. A host using a densitized baryon variable must
convert its normalization before calling this interface.

Both species' outputs are initialized to their source bases with zero exchange
packets. Validation, source solves, bounds checks, and exchange assembly must
all succeed before either candidate is published. Number-floor injections
are not part of this conserving operation; a host that repairs input states
must track that repair separately. The host applies a common admissible
limiter to both radiation and all matter increments before matter recovery.

Single-species source operations cannot evaluate electron-flavor pair
collisions. They reject separated pair rates and old aggregate electron
number rates containing a non-charged-current contribution. Callers must use
the paired source operation for these rates.
