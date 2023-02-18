#include "con2prim.h"

/* Function    : enforce_primitive_limits_and_compute_u0()
 * Description : Applies limits to rho_b, pressure, and v^i, then
                 recomputes epsilon and (if needed) entropy

 * Inputs      : params         - GRHayL_parameters struct with parameters
 *                                for the simulation
 *             : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - primitive_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : prims          - returns primitives within floors and ceilings
 *             : diagnostics    - tracks whether velocity was speed-limited
 *
 */

void enforce_primitive_limits_and_compute_u0(const GRHayL_parameters *restrict params,
                                            const eos_parameters *restrict eos,
                                            const metric_quantities *restrict metric,
                                            primitive_quantities *restrict prims,
                                            con2prim_diagnostics *restrict diagnostics) {

  // The density floor and ceiling is always applied
  prims->rho = MIN(MAX(prims->rho,eos->rho_min),eos->rho_max);

  // Hybrid EOS specific floors and ceilings
  if( eos->eos_type == grhayl_eos_hybrid ) {
    // Pressure and epsilon must be recomputed
    // Compute P and eps
    double P_cold = 0.0;
    double eps_cold = 0.0;
    eos->hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &P_cold, &eps_cold);

    // Set P_min and P_max
    const double P_min = 0.9*P_cold;
    // Adjust P_max based on Psi6
    const double P_max = (metric->psi6 > params->psi6threshold) ? 1e5*P_cold : 100.0*P_cold;

    // Now apply floors and ceilings to P
    if(prims->press<P_min) prims->press = P_min;
    // Finally, perform the last check
    if((prims->rho < 100.0*eos->rho_atm || metric->psi6 > params->psi6threshold) && prims->press>P_max) {
      prims->press = P_max;
    }
    // Now recompute eps and, if needed, entropy
    prims->eps = eps_cold + (prims->press-P_cold)/(eos->Gamma_th-1.0)/prims->rho;
    if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims->rho, prims->press, &prims->entropy);

  // Tabulated EOS specific floors and ceilings
  } else if( eos->eos_type == grhayl_eos_tabulated ) {
    grhayl_warn("No tabulated EOS support yet! Sorry!");
  //  // Apply floors and ceilings to Y_e and T
  //  const double xye   = MIN(MAX(prims->Y_e,eos->Ye_min),eos->Ye_atm);
  //  const double xtemperature = MIN(MAX(prims->temperature,eos->T_atm ),eos->T_max );

  //  // Additional variables used for the EOS call
  //  const double xrho  = prims->rho;
  //  double xprs        = 0.0;
  //  double xeps        = 0.0;
  //  double xeps        = 0.0;
  //  double xent        = 0.0;
  //  WVU_EOS_P_eps_and_S_from_rho_Ye_T(xrho,xye,xtemp, &xprs,&xeps,&xent);
  //  // Now update the primitives (rho has already been set)
  //  prims->Y_e = xye;
  //  prims->temperature = xtemp;
  //  prims->press = xprs;
  //  prims->eps = xeps;
  //  prims->entropy = xent;
  }

  // Finally, apply speed limit to v and compute u^0
  limit_v_and_compute_u0(eos, metric, prims, diagnostics);
}
