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
      const primitive_quantities *restrict prims,
      primitive_quantities *restrict prims_guess,
      con2prim_diagnostics *restrict diagnostics) {

  if( eos->eos_type != 0 ) {
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
  diagnostics->font_fixes++;


  //Translate to HARM primitive now:
  double utcon1 = metric->adm_gupxx*u_xl + metric->adm_gupxy*u_yl + metric->adm_gupxz*u_zl;
  double utcon2 = metric->adm_gupxy*u_xl + metric->adm_gupyy*u_yl + metric->adm_gupyz*u_zl;
  double utcon3 = metric->adm_gupxz*u_xl + metric->adm_gupyz*u_yl + metric->adm_gupzz*u_zl;

  double u0;
  //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
  limit_utilde_and_compute_v(eos, metric, &u0, &utcon1, &utcon2, &utcon3, prims_guess, diagnostics);

  prims_guess->rho = cons->rho/(metric->lapse*u0*metric->psi6);

  double K_ppoly, Gamma_ppoly;
  eos->hybrid_get_K_and_Gamma(eos, prims_guess->rho, &K_ppoly, &Gamma_ppoly);

  // After that, we compute P_cold
  double P_cold = K_ppoly*pow(prims_guess->rho, Gamma_ppoly);
  double energy_u = P_cold/(Gamma_ppoly-1.0);

  prims_guess->press = pressure_rho0_u(eos, prims_guess->rho, energy_u);
  prims_guess->eps = energy_u/prims_guess->rho;
  if( params->evolve_entropy ) eos->hybrid_compute_entropy_function(eos, prims_guess->rho, prims_guess->press, &prims_guess->entropy);
  return 0;
}
