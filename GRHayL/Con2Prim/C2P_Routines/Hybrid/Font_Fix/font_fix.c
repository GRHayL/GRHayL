#include "con2prim.h"

/**********************************
 * Piecewise Polytropic EOS Patch *
 *    Font fix: function call     *
 **********************************/
int font_fix(
      const GRHayL_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics) {


  // Use the new Font fix subroutine
  double u_x=1e100, u_y=1e100, u_z=1e100; // Set to insane values to ensure they are overwritten.
  /************************
   * New Font fix routine *
   ************************/
  int check = font_fix_hybrid_EOS(eos, metric, cons, prims, &u_x, &u_y, &u_z);
  diagnostics->failure_checker+=10000;
  diagnostics->font_fix=1;

  //Translate to HARM primitive now:
  double utcon1 = metric->adm_gupxx*u_x + metric->adm_gupxy*u_y + metric->adm_gupxz*u_z;
  double utcon2 = metric->adm_gupxy*u_x + metric->adm_gupyy*u_y + metric->adm_gupyz*u_z;
  double utcon3 = metric->adm_gupxz*u_x + metric->adm_gupyz*u_y + metric->adm_gupzz*u_z;

  //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
  limit_utilde_and_compute_v(eos, metric, &utcon1, &utcon2, &utcon3, prims, &diagnostics->speed_limited);
  prims->rho = cons->rho/(metric->lapse*prims->u0*metric->psi6);

  double K_ppoly, Gamma_ppoly;
  eos->hybrid_get_K_and_Gamma(eos, prims->rho, &K_ppoly, &Gamma_ppoly);

  // After that, we set P = P_cold
  eos->hybrid_compute_P_cold(eos, prims->rho, &prims->press);

  // and compute epsilon from rho and pressure
  prims->eps = prims->press/(prims->rho*(Gamma_ppoly-1.0));
  if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims->rho, prims->press, &prims->entropy);
  return check;
}
