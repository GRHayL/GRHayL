#include "../../../utils_Noble.h"

/**
 * @ingroup Con2Prim
 * @brief Solves for primitive variables with the density-based entropy
 *        variant of the Noble1D method.
 *
 * @details Solves the momentum equation directly for rest-mass density,
 * using the entropy conservative variable to close the hybrid EOS. The
 * energy equation is not used.
 *
 * @param[in] params pointer to ghl_parameters struct
 * @param[in] eos pointer to ghl_eos_parameters struct
 * @param[in] metric_adm pointer to ghl_metric_quantities struct with ADM metric
 * @param[in] metric_aux pointer to ghl_ADM_aux_quantities struct
 * @param[in] cons_undens pointer to undensitized conservative variables
 * @param[in,out] prims initial guess on input and recovered primitives on output
 * @param[out] diagnostics Con2Prim diagnostics, written on successful recovery
 *
 * @returns ghl_success or a specific Con2Prim failure code
 */

ghl_error_codes_t ghl_hybrid_Noble1D_entropy2(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double gnr_out[1];

  harm_aux_vars_struct harm_aux;

  double rho0, Z_last;
  ghl_error_codes_t error = ghl_initialize_Noble_entropy(
        params, eos, metric_adm, metric_aux, cons_undens, prims, &harm_aux, &rho0,
        &Z_last);
  if(error) {
    return error;
  }
  if(rho0 <= 0.0) {
    return ghl_error_neg_rho;
  }

  gnr_out[0] = rho0;

  const ghl_error_codes_t retval = ghl_general_newton_raphson(
        eos, &harm_aux, 1, Z_last, gnr_out, ghl_validate_1D_entropy, ghl_func_rho2);

  rho0 = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if(retval != ghl_success) {
    return retval;
  }
  // Calculate v^2:

  const double rel_err = fabs((harm_aux.D - rho0) / harm_aux.D);
  const double utsq = (rel_err > 1e-15)
                            ? (harm_aux.D - rho0) * (harm_aux.D + rho0) / (rho0 * rho0)
                            : 0.0;
  if(utsq < 0.0) {
    return ghl_error_neg_vsq;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  const double Wsq = 1.0 + utsq;
  const double W = sqrt(Wsq);

  const double Gamma_ppoly
        = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos, rho0)];
  const double p_final = cons_undens->entropy * pow(rho0, Gamma_ppoly) / harm_aux.D;
  const double eps_final = ghl_hybrid_compute_epsilon(eos, rho0, p_final);
  const double w = rho0 * (1.0 + eps_final) + p_final;
  const double Z = w * Wsq;

  prims->rho = rho0;

  diagnostics->speed_limited |= ghl_finalize_Noble_entropy(
        params, eos, metric_adm, metric_aux, cons_undens, &harm_aux, Z, W, prims);
  if(prims->press <= 0.0) {
    return ghl_error_neg_pressure;
  }

  /* Done! */
  diagnostics->n_iter = harm_aux.n_iter;
  diagnostics->which_routine = ghl_con2prim_id_Noble1D_entropy2;
  return ghl_success;
}
