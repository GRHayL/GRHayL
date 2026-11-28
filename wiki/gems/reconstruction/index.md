# Reconstruction

> Router for GRHayL face reconstruction: left-face stencil conventions, PLM limiters, PPM flattening/steepening, WENO-z, and Reconstruction tests.

| Page | Go here when... |
| --- | --- |
| [Face And Stencil Conventions](face-and-stencil-conventions.md) | You need left-face `Ur`/`Ul` orientation at `i-1/2`, 4-/5-/6-point stencils, `enum reconstruction_stencil`, or wrapper-vs-helper semantics. |
| [PLM Limiters](plm-limiters.md) | You need `ghl_minmod_reconstruction`, `ghl_mc_reconstruction`, `ghl_superbee_reconstruction`, and shared `ghl_minmod`/`ghl_maxmod` helpers. |
| [PPM Flow](ppm-flow.md) | You need `ghl_ppm_reconstruction`, `ghl_ppm_reconstruction_with_steepening`, `ghl_compute_ftilde`, flattening, steepening, or `ppm_*` parameters. |
| [WENO-z Method](wenoz-method.md) | You need `ghl_wenoz_reconstruction`, `ghl_wenoz_reconstruction_right_left_faces`, fixed weights/regularization, or positivity/overflow caveats. |
| [Reconstruction Tests And Reference Data](tests-and-reference-data.md) | You need PLM, WENO-z, and ET_Legacy PPM tests, generators, CI vs. `run_tests.sh` differences, and coverage bounds. |

Back to [GRHayL Gems](../index.md).
