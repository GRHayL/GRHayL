# Con2Prim

> Router for conservative-to-primitive solvers, multi-method recovery, primitive/conservative limits, neural-network guesses, root-finding numerics, and Con2Prim tests.

| Page | Go here when... |
| --- | --- |
| [Solver Matrix](solver-matrix.md) | You need `ghl_con2prim_id_t` support status across hybrid/tabulated EOS, selector/build seams, `ghl_con2prim_diagnostics`, or archival Cerda-Duran status. |
| [Recovery Flow](recovery-flow.md) | You need the ordered 7-step recovery sequence, densitized vs. undensitized boundaries, backup retries, NN retry hook, or `GHL_DISABLE_HDF5` stubs. |
| [Limits And Conversions](limits-and-conversions.md) | You need `ghl_apply_conservative_limits`, `ghl_guess_primitives`, `ghl_enforce_primitive_limits_and_compute_u0`, `ghl_limit_utilde_and_compute_v`, or `ghl_compute_conservs*`. |
| [Neural-Network Primitive Guess](neural-network-primitive-guess.md) | You need `GHL_NN_C2P_API_VERSION`, `ghl_c2p_nn_model` HDF5 layout (`grhayl_nn_c2p`), overflow-safe `x` inference, fallback behavior, or `enable_backup_nn_primitive_guess`. |
| [Internal Numerics And Root-Finding](internal-numerics-and-root-finding.md) | You need `ghl_brent`, Noble Newton validation/closure, Palenzuela `evolve_temp` callbacks, Font1D cycle detection, Newman iteration counters, or Cerda-Duran archival notes. |
| [Con2Prim Tests And Reference Data](tests-and-reference-data.md) | You need hybrid/tabulated/NN/ET_Legacy Con2Prim tests, Newman history sanitizer requirements, `pyghl append` setup, and targeted check commands. |

Back to [GRHayL Gems](../index.md).
