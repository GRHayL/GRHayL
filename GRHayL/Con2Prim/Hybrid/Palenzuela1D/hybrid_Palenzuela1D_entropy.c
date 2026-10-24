#include "../../utils_Palenzuela1D.h"

/**
 * @ingroup c2p_internal
 * @brief Evaluate hybrid thermodynamics and W for the entropy residual.
 * @author Leo Werneck
 *
 * Derives pressure from conservative entropy, floors it at the cold pressure,
 * and evaluates specific internal energy from the hybrid EOS.
 *
 * @param[in] x Trial value of hW
 * @param[in] params Recovery parameters supplying the Lorentz-factor cap
 * @param[in] eos Hybrid EOS parameters
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] fparams Normalized contractions and callback state from roots.h
 * @param[in,out] prims Primitive state receiving rho, pressure, and internal energy
 * @param[out] W_ptr Limited Lorentz factor
 *
 * See Palenzuela et al. (2015), arXiv:1505.01607, and Siegel et al. (2018),
 * arXiv:1712.07538. Entropy formulation: Werneck et al. (2023), arXiv:2208.14487.
 */
static void
compute_rho_P_eps_W_entropy(
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

  // Step 2: Compute P (Eq. 19 of arXiv:0808.3140)
  prims->entropy     = cons_undens->entropy / W;
  const double Gamma = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos, prims->rho)];
  prims->press       = prims->entropy * pow(prims->rho, Gamma - 1.0);

  // Step 3: Enforce the cold-pressure floor and compute specific internal energy
  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &P_cold, &eps_cold);
  // Don't let P be less than P_cold
  if(prims->press < P_cold) prims->press = P_cold;
  // Compute eps
  prims->eps = eps_cold + (prims->press - P_cold) / ( (eos->Gamma_th - 1.0) * prims->rho);

  // Step 4: Set the output
  *W_ptr   = W;
}

/* ghl_hybrid_Palenzuela1D_entropy selects the entropy callback for
 * ghl_hybrid_Palenzuela1D. See ghl_con2prim.h for the public calling contract.
 */
ghl_error_codes_t ghl_hybrid_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  const ghl_error_codes_t error = ghl_hybrid_Palenzuela1D(
        compute_rho_P_eps_W_entropy, params, eos, metric_adm, cons_undens, prims,
        diagnostics);
  if(error == ghl_success) {
    diagnostics->which_routine = ghl_con2prim_id_Palenzuela1D_entropy;
  }
  return error;
}
