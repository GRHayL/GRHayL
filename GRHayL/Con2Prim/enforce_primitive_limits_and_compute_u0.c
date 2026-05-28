#include "ghl_con2prim.h"

/**
 * @ingroup Con2Prim
 * @brief Applies limits to density, pressure, and velocity, then
          recomputes epsilon and (if needed) entropy.
 *
 * @details
 * This function applies limits to the primitive variables, generally used
 * after solving for them with a conservative-to-primitive method. It
 * then speed-limits the velocity and solves for the 0th component of the
 * 4-velocity, which is needed for computing the conservatives and stress-
 * energy tensor.
 *
 * @param[in] params:     pointer to ghl_parameters struct
 *
 * @param[in] eos:        pointer to ghl_eos_parameters struct
 *
 * @param[in] ADM_metric: pointer to ghl_metric_quantities struct with ADM metric
 *
 * @param[in,out] prims: pointer to ghl_primitive_quantities struct
 *
 * @returns GRHayL error code
 */
ghl_error_codes_t ghl_enforce_primitive_limits_and_compute_u0(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims,
      bool *restrict speed_limited) {

  /**
   * # Step-by-step Procedure
   * First, the code splits based on EOS type.
   */
  switch(eos->eos_type) {
    case ghl_eos_hybrid:
    {
      /**
       * ## Hybrid EOS
       * The code
       * 1. Applies floors and ceilings to the density
       */
      prims->rho = ghl_clamp(prims->rho, eos->rho_min, eos->rho_max);

      /**
       * 2. Recomputes the cold pressure and epsilon
       */
      double P_cold;
      double eps_cold;
      ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &P_cold, &eps_cold);

      /**
       * 3. Sets pressure limits:
       *    \f$ P_{\mathrm{min}} = P_{\mathrm{cold}} \f$
       *
       *    \f$
       *    P_{\mathrm{max}} = \begin{cases}
       *    10^5 P_{\mathrm{cold}} & \psi^6 > \psi^6_\mathrm{threshold} \\
       *    100P_{\mathrm{cold}} & \mathrm{otherwise}
       *    \end{cases}
       *    \f$
       *    For the maximum pressure, the estimate is dependent on whether the point
       *    is approximately within the horizon.
       */
      const double P_min = P_cold;

      const bool inhorizon = ADM_metric->sqrt_detgamma > params->psi6threshold;
      const double P_max = inhorizon ? 1e5 * P_cold : 100.0 * P_cold;

      /**
       * 4. Apply pressure limits. Again, the upper limit is more complicated. It
       *    is only applied if either the point is in the horizon or the density
       *    is small (\f$ < 100 \rho_\mathrm{atm} \f$).
       */
      if(prims->press < P_min) {
        prims->press = P_min;
      }
      else if((prims->rho < 100.0 * eos->rho_atm || inhorizon) && prims->press > P_max) {
        prims->press = P_max;
      }

      /**
       * 4. For consistency, recompute eps and, if needed, entropy.
       */ 
      prims->eps = eps_cold + (prims->press - P_cold) / (eos->Gamma_th - 1.0) / prims->rho;
      if(params->evolve_entropy) {
        prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);
      }
      break;
    }

    case ghl_eos_tabulated:
      /**
       * ## Tabulated EOS
       * The code
       * 1. Applies floors and ceilings to density, \f$ Y_e \f$, and temperature.
       */
      ghl_tabulated_enforce_bounds_rho_Ye_T(
            eos, &prims->rho, &prims->Y_e, &prims->temperature);

      /**
       * 2. Compute pressure, epsilon, and optionally entropy.
       */ 
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
      /**
       * ## Simple EOS
       * The code
       * 1. Applies floors and ceilings to density and pressure.
       */
      prims->rho   = ghl_clamp(prims->rho,   eos->rho_min,   eos->rho_max);
      prims->press = ghl_clamp(prims->press, eos->press_min, eos->press_max);

      /**
       * 2. For consistency, recompute eps and, if needed, entropy.
       */ 
      prims->eps = prims->press / (prims->rho * (eos->Gamma_th - 1.0));
      if(params->evolve_entropy) {
        prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);
      }
      break;
    default:
      return ghl_error_unknown_eos_type;
  }

  /** Finally, we speed limit the velocity and compute \f$ u^0 \f$ */
  return ghl_limit_v_and_compute_u0(params, ADM_metric, prims, speed_limited);
}
