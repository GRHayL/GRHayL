#include "con2prim.h"

/* Function    : ghl_enforce_primitive_limits_and_compute_u0()
 * Description : Applies limits to rho_b, pressure, and v^i, then
                 recomputes epsilon and (if needed) entropy

 * Inputs      : params         - ghl_parameters struct with parameters
 *                                for the simulation
 *             : eos            - ghl_eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - ghl_metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - ghl_primitive_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : prims          - returns primitives within floors and ceilings
 *             : speed_limited  - tracks whether velocity was speed-limited
 *
 */

int ghl_enforce_primitive_limits_and_compute_u0(
    const ghl_parameters *restrict params,
    const ghl_eos_parameters *restrict eos,
    const ghl_metric_quantities *restrict ADM_metric,
    ghl_primitive_quantities *restrict prims) {

  // The density floor and ceiling is always applied
  prims->rho = MIN(MAX(prims->rho,eos->rho_min),eos->rho_max);

  // Hybrid EOS specific floors and ceilings
  switch(eos->eos_type) {
    case ghl_eos_hybrid: {
      // Pressure and epsilon must be recomputed
      // Compute P and eps
      double P_cold;
      double eps_cold;
      ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &P_cold, &eps_cold);

      // Set P_min and P_max
      const double P_min = P_cold;
      // Adjust P_max based on Psi6
      const double P_max = (ADM_metric->sqrt_detgamma > params->psi6threshold) ? 1e5*P_cold : 100.0*P_cold;

      // Now apply floors and ceilings to P
      if(prims->press < P_min) prims->press = P_min;
      else if((prims->rho < 100.0*eos->rho_atm || ADM_metric->sqrt_detgamma > params->psi6threshold) && prims->press > P_max)
        prims->press = P_max;

      // Now recompute eps and, if needed, entropy
      prims->eps = eps_cold + (prims->press-P_cold)/(eos->Gamma_th-1.0)/prims->rho;
      if(params->evolve_entropy)
        prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);
      break;}

    // Tabulated EOS specific floors and ceilings
    case ghl_eos_tabulated:
      // Apply floors and ceilings to Y_e and T
      prims->Y_e         = MIN(MAX(prims->Y_e,eos->Y_e_min), eos->Y_e_max);
      prims->temperature = MIN(MAX(prims->temperature, eos->T_min),eos->T_max);

      // Additional variables used for the EOS call
      if(params->evolve_entropy) {
        ghl_tabulated_compute_P_eps_S_from_T(eos,
                                             prims->rho, prims->Y_e, prims->temperature,
                                             &prims->press, &prims->eps, &prims->entropy);
      }
      else {
        ghl_tabulated_compute_P_eps_from_T(eos,
                                           prims->rho, prims->Y_e, prims->temperature,
                                           &prims->press, &prims->eps);
      }
      break;

    case ghl_eos_simple:
      // Apply floors and ceilings to P
      if(prims->press < eos->press_min) prims->press = eos->press_min;
      else if(prims->press > eos->press_max) prims->press = eos->press_max;

      // Now recompute eps and, if needed, entropy
      prims->eps = prims->press/(prims->rho * (eos->Gamma_th-1.0));
      if(params->evolve_entropy)
        prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);
      break;
    default:
      ghl_error("Unknown EOS type %d\n", eos->eos_type);
  }

  // Finally, apply speed limit to v and compute u^0
  return ghl_limit_v_and_compute_u0(params, ADM_metric, prims);
}
