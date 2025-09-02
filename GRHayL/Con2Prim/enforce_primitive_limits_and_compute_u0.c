#include "ghl_con2prim.h"

/* Function     : ghl_enforce_primitive_limits_and_compute_u0()
 * Description  : Applies limits to rho_b, pressure, and v^i, then
                  recomputes epsilon and (if needed) entropy
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_enforce_primitive_limits_and_compute_u0
*/

ghl_error_codes_t ghl_enforce_primitive_limits_and_compute_u0(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims,
      bool *restrict speed_limited) {

  // Hybrid EOS specific floors and ceilings
  switch(eos->eos_type) {
    case ghl_eos_hybrid:
    {
      // Apply floors and ceilings to rho
      prims->rho = MIN(MAX(prims->rho, eos->rho_min), eos->rho_max);

      // Pressure and epsilon must be recomputed
      // Compute P and eps
      double P_cold;
      double eps_cold;
      ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &P_cold, &eps_cold);

      // Set P_min and P_max
      const double P_min = P_cold;

      const bool inhorizon = ADM_metric->sqrt_detgamma > params->psi6threshold;
      // Adjust P_max based on Psi6
      const double P_max = inhorizon ? 1e5 * P_cold : 100.0 * P_cold;

      // Now apply floors and ceilings to P
      if(prims->press < P_min) {
        prims->press = P_min;
      }
      else if((prims->rho < 100.0 * eos->rho_atm || inhorizon) && prims->press > P_max) {
        prims->press = P_max;
      }

      // Now recompute eps and, if needed, entropy
      prims->eps = eps_cold + (prims->press - P_cold) / (eos->Gamma_th - 1.0) / prims->rho;
      if(params->evolve_entropy) {
        prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);
      }
      break;
    }

    // Tabulated EOS specific floors and ceilings
    case ghl_eos_tabulated:
      // Apply floors and ceilings to rho, Y_e and T
      ghl_tabulated_enforce_bounds_rho_Ye_T(
            eos, &prims->rho, &prims->Y_e, &prims->temperature);

      // Additional variables used for the EOS call
      if(params->evolve_entropy) {
        ghl_tabulated_compute_P_eps_S_from_T(
              eos, prims->rho, prims->Y_e, prims->temperature, &prims->press,
              &prims->eps, &prims->entropy);
      }
      else {
        ghl_tabulated_compute_P_eps_from_T(
              eos, prims->rho, prims->Y_e, prims->temperature, &prims->press,
              &prims->eps);
      }
      break;

    case ghl_eos_simple:
      // Apply floors and ceilings to rho
      prims->rho = MIN(MAX(prims->rho, eos->rho_min), eos->rho_max);

      // Apply floors and ceilings to P
      prims->press = MIN(MAX(prims->press, eos->press_min), eos->press_max);

      // Now recompute eps and, if needed, entropy
      prims->eps = prims->press / (prims->rho * (eos->Gamma_th - 1.0));
      if(params->evolve_entropy) {
        prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);
      }
      break;
    default:
      return ghl_error_unknown_eos_type;
  }

  // Finally, apply speed limit to v and compute u^0
  return ghl_limit_v_and_compute_u0(params, ADM_metric, prims, speed_limited);
}
