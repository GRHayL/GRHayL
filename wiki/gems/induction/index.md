# Induction

> Router for vector-potential magnetic HLL fluxes, ADM/BSSN staggered interpolators, Lorenz scalar-potential gauge RHS, tests, and verification commands.

| Page | Go here when... |
| --- | --- |
| [HLL Flux](hll-flux.md) | You need `ghl_HLL_vars` cross-product packing, `ghl_HLL_flux_with_B` vs. `ghl_HLL_flux_with_Btilde`, sample offsets, or denominator hazards. |
| [Interpolation And Staggering](interpolation-and-staggering.md) | You need cell-centered ADM/BSSN and vertex-centered ADM `2x2x2`/`3x3x3` stencil layouts, `ghl_induction_interp_vars`, or downstream GRHayLET alignment. |
| [Gauge RHS](gauge-rhs.md) | You need `ghl_calculate_phitilde_rhs`, shift upwinding, `sqrtg_Ai_stencil[3][2]`, `Lorenz_damping_factor`, or caller-assembled `A_i` gauge terms. |
| [Induction Tests And Reference Data](tests-and-reference-data.md) | You need Induction HLL, `ccc_ADM`/`ccc_BSSN`/`vvv_ADM`, and ET_Legacy gauge RHS tests, private helpers, generators, and coverage gaps. |
| [Induction Verification Workflows](verification-workflows.md) | You need targeted local Induction build/run commands, CI matrix mapping, and `generate_makefile.sh` caveats. |

Back to [GRHayL Gems](../index.md).
