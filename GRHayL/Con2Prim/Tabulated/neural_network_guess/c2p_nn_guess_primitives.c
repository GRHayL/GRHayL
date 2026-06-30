#include "ghl_c2p_nn.h"

void ghl_c2p_nn_guess_primitives(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims) {

  ghl_tabulated_primitive_guess_aux aux;
  ghl_tabulated_compute_primitive_guess_auxiliaries(metric_adm, cons_undens, prims, &aux);

  // Guess x and W using the neural network
  const ghl_nn_c2p_input_t nn_input = { aux.q, aux.r, aux.s, aux.t };
  const ghl_nn_c2p_guess_t nn_guess = ghl_c2p_nn_guess(eos->c2p_nn, nn_input);

  // Enforce physical bounds on x
  const double x_lo = 1.0 + aux.q - aux.s;
  const double x_hi = 2.0 + 2.0 * aux.q - aux.s;
  const double x = ghl_clamp(nn_guess.x, x_lo, x_hi);

  // Complete the primitive guess using the neural-network x and W guesses. If
  // the W guess is unavailable, the helper falls back to Eq. (42) of 1712.07538.
  ghl_tabulated_primitive_guess_from_x_and_W(params, eos, metric_adm,
                                             cons_undens, &aux,
                                             x, nn_guess.W, prims);
}
