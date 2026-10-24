#include "../../utils_Palenzuela1D.h"

/**
 * @ingroup c2p_internal
 * @brief Evaluate tabulated thermodynamics and W for the entropy residual.
 * @author Leo Werneck
 *
 * When temperature is evolved, uses entropy to invert the EOS for
 * temperature. Otherwise evaluates the EOS at the supplied temperature.
 * The incoming primitive electron fraction is bounded, not reconstructed.
 * Records the EOS evaluation status in fparams->eos_error.
 *
 * @param[in] x Trial value of hW
 * @param[in] params Recovery parameters supplying the Lorentz-factor cap
 * @param[in] eos Tabulated EOS parameters
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] fparams Normalized contractions and callback state from roots.h
 * @param[in,out] prims Primitive state receiving rho, pressure, and internal energy
 * @param[out] W_ptr Limited Lorentz factor
 *
 * See Palenzuela et al. (2015), arXiv:1505.01607, and Siegel et al. (2018),
 * arXiv:1712.07538. Entropy formulation: Werneck et al. (2023), arXiv:2208.14487.
 */
static void
compute_rho_P_eps_T_W_entropy(
      const double x,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims,
      double *restrict W_ptr) {

  // Step 1: First compute rho and W
  double W;
  compute_rho_W_from_x_and_conservatives(x, params, cons_undens, fparams, &prims->rho, &W);

  // Step 2: Now compute P and eps
  ghl_error_codes_t error;
  if(fparams->evolve_T) {
    prims->entropy = cons_undens->entropy/W;
    ghl_tabulated_enforce_bounds_rho_Ye_S(eos, &prims->rho, &prims->Y_e, &prims->entropy);
    error = ghl_tabulated_compute_P_eps_T_from_S(eos, prims->rho, prims->Y_e, prims->entropy, &prims->press, &prims->eps, &prims->temperature);
  }
  else {
    // If the temperature is not evolved, use the input guess to determine
    // the remaining primitives. Note that in this case one must provide
    // the appropriate temperature instead of the default guess T = T_min.
    ghl_tabulated_enforce_bounds_rho_Ye_T(eos, &prims->rho, &prims->Y_e, &prims->temperature);
    error = ghl_tabulated_compute_P_eps_from_T(eos, prims->rho, prims->Y_e, prims->temperature, &prims->press, &prims->eps);
  }

  // Step 3: Store the EOS status of this evaluation for the caller
  fparams->eos_error = error;

  // Step 4: Set the output
  *W_ptr   = W;
}

/* ghl_tabulated_Palenzuela1D_entropy selects the entropy callback for
 * ghl_tabulated_Palenzuela1D. See ghl_con2prim.h for the public calling contract.
 */
ghl_error_codes_t ghl_tabulated_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  const ghl_error_codes_t error = ghl_tabulated_Palenzuela1D(
        compute_rho_P_eps_T_W_entropy, params, eos, metric_adm, cons_undens, prims,
        diagnostics);
  if(error == ghl_success) {
    diagnostics->which_routine = ghl_con2prim_id_Palenzuela1D_entropy;
  }
  return error;
}
