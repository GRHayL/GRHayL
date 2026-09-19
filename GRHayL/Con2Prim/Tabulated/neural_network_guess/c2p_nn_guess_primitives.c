#include "ghl_c2p_nn.h"

/**
 * @brief Build a tabulated primitive initial guess from a neural network.
 *
 * A present model must already be validated. With usable finite contractions
 * and candidate it takes the existing NN path; a missing model or unusable
 * numerical state/candidate returns the atmosphere initial guess. The fallback
 * preserves initialized finite magnetic components and supplies zero Eulerian
 * velocity; it is not a successful primitive recovery. Arguments must be
 * non-NULL, and parameters, EOS atmosphere values, and metric lapse/shift must
 * already satisfy their public contracts. Arbitrarily corrupted EOS or model
 * objects are outside this helper's contract.
 */
void ghl_c2p_nn_guess_primitives(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims) {

  if(eos->c2p_nn != NULL && isfinite(cons_undens->rho) && cons_undens->rho > 0.0) {
    ghl_tabulated_primitive_guess_aux aux;
    ghl_tabulated_compute_primitive_guess_auxiliaries(metric_adm, cons_undens, prims, &aux);

    if(isfinite(aux.q) && isfinite(aux.r)
          && isfinite(aux.s) && isfinite(aux.t)) {
      // Guess x using the neural network
      const ghl_nn_c2p_input_t nn_input = { aux.q, aux.r, aux.s, aux.t };
      const ghl_nn_c2p_guess_t nn_guess = ghl_c2p_nn_guess(eos->c2p_nn, nn_input);

      // Enforce the NN training bracket before completing the primitive guess.
      const double x_lo = 1.0 + aux.q - aux.s;
      const double x_hi = 2.0 + 2.0 * aux.q - aux.s;
      if(isfinite(x_lo) && isfinite(x_hi) && x_lo <= x_hi && isfinite(nn_guess.x)) {
        const double x = ghl_clamp(nn_guess.x, x_lo, x_hi);
        if(isfinite(x)) {
          ghl_tabulated_primitive_guess_from_x(params, eos, metric_adm,
                                               cons_undens, &aux, x, prims);
          return;
        }
      }
    }
  }

  // Reuse the shared completion helper's atmosphere fallback.
  const ghl_tabulated_primitive_guess_aux aux = { 0 };
  ghl_tabulated_primitive_guess_from_x(
        params, eos, metric_adm, cons_undens, &aux, NAN, prims);
}
