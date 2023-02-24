#include "con2prim.h"
#include "../harm_u2p_util.h"

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

  if( eos->eos_type != grhayl_eos_hybrid ) {
    grhayl_warn("Font fix is only implemented for Hybrid EOS! The eos_type is set to %d.", eos->eos_type); //TODO: make a better error message
    return 1;
  }

  // Use the new Font fix subroutine
  double u_xl=1e100, u_yl=1e100, u_zl=1e100; // Set to insane values to ensure they are overwritten.
  /************************
   * New Font fix routine *
   ************************/
  int check = font_fix_hybrid_EOS(eos, metric, cons, prims, &u_xl, &u_yl, &u_zl);
  diagnostics->failure_checker+=10000;
  diagnostics->font_fixes=1;

  //Translate to HARM primitive now:
  double utcon1 = metric->adm_gupxx*u_xl + metric->adm_gupxy*u_yl + metric->adm_gupxz*u_zl;
  double utcon2 = metric->adm_gupxy*u_xl + metric->adm_gupyy*u_yl + metric->adm_gupyz*u_zl;
  double utcon3 = metric->adm_gupxz*u_xl + metric->adm_gupyz*u_yl + metric->adm_gupzz*u_zl;

  //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
  limit_utilde_and_compute_v(eos, metric, &utcon1, &utcon2, &utcon3, prims, diagnostics);
  prims->rho = cons->rho/(metric->lapse*prims->u0*metric->psi6);

  double K_ppoly, Gamma_ppoly;
  eos->hybrid_get_K_and_Gamma(eos, prims->rho, &K_ppoly, &Gamma_ppoly);

  // After that, we set P = P_cold
  eos->hybrid_compute_P_cold(eos, prims_guess->rho, &prims_guess->press);

  // and compute epsilon from rho and pressure
  prims_guess->eps = prims_guess->press/(prims_guess->rho*(Gamma_ppoly-1.0));
  if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims_guess->rho, prims_guess->press, &prims_guess->entropy);
  return check;
}
