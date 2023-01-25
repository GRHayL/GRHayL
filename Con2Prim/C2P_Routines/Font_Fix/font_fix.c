#include "con2prim.h"
#include "../harm_u2p_util.h"
#include <stdio.h>

/**********************************
 * Piecewise Polytropic EOS Patch *
 *    Font fix: function call     *
 **********************************/
int font_fix(const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons,
             const primitive_quantities *restrict prims,
             primitive_quantities *restrict prims_guess,
             con2prim_diagnostics *restrict diagnostics) {

  int check=1;

  double u0;

  double utcon1;
  double utcon2;
  double utcon3;

  if( eos->eos_type == 0 ) {
    // Use the new Font fix subroutine
    double u_xl=1e100, u_yl=1e100, u_zl=1e100; // Set to insane values to ensure they are overwritten.
    /************************
     * New Font fix routine *
     ************************/
    check = font_fix_hybrid_EOS(eos, metric, cons, prims, &u_xl, &u_yl, &u_zl);

    //Translate to HARM primitive now:
    utcon1 = metric->adm_gupxx*u_xl + metric->adm_gupxy*u_yl + metric->adm_gupxz*u_zl;
    utcon2 = metric->adm_gupxy*u_xl + metric->adm_gupyy*u_yl + metric->adm_gupyz*u_zl;
    utcon3 = metric->adm_gupxz*u_xl + metric->adm_gupyz*u_yl + metric->adm_gupzz*u_zl;

    if (check==1) {
//TODO: error checking
//      CCTK_VInfo(CCTK_THORNSTRING,"Font fix failed!");
//          CCTK_VInfo(CCTK_THORNSTRING,"i,j,k = %d %d %d, stats.failure_checker = %d x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",i,j,k,stats.failure_checker,X[index],Y[index],Z[index],index,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6]);
      return check;
    }
  /*************************************************************/
  } else {
    grhayl_warn("Font fix is only implemented for Hybrid EOS! The eos_type is set to %d.", eos->eos_type); //TODO: make a better error message
    return 5;
  } //EOS if

  diagnostics->failure_checker+=10000;
  diagnostics->font_fixes++;
  //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
  limit_utilde_and_compute_v(eos, metric, &u0, &utcon1, &utcon2, &utcon3, prims_guess, diagnostics);

  prims_guess->rho = cons->rho/(metric->lapse*u0*metric->psi6);

  double K_ppoly, Gamma_ppoly;
  eos->hybrid_get_K_and_Gamma(eos, prims_guess->rho, &K_ppoly, &Gamma_ppoly);

  // After that, we compute P_cold
  double P_cold = K_ppoly*pow(prims_guess->rho, Gamma_ppoly);

  double energy_u = P_cold/(Gamma_ppoly-1.0);

  if( eos->eos_type == 0 ) {
    prims_guess->press = pressure_rho0_u(eos, prims_guess->rho, energy_u);
  }
  return 0;
}
