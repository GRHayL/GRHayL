#include "con2prim_header.h"
#include "EOS_hybrid_header.h"
#include "harm_u2p_util.h"
#include <stdio.h>

/**********************************
 * Piecewise Polytropic EOS Patch *
 *    Font fix: function call     *
 **********************************/
int font_fix(const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const conservative_quantities *restrict cons_undens,
             const primitive_quantities *restrict prims,
             primitive_quantities *restrict prims_guess,
             con2prim_diagnostics *restrict diagnostics,
             double *restrict u0L_ptr) {

  int check=1;

  double u0L = *u0L_ptr;

  double utcon1;
  double utcon2;
  double utcon3;

  if( eos->eos_type == 0 ) {
    // Use the new Font fix subroutine
    double u_xl=1e100, u_yl=1e100, u_zl=1e100; // Set to insane values to ensure they are overwritten.
    /************************
     * New Font fix routine *
     ************************/
    check = font_fix_hybrid_EOS(eos, metric, cons_undens, prims, &u_xl, &u_yl, &u_zl);

    //Translate to HARM primitive now:
    utcon1 = metric->adm_gxx*u_xl + metric->adm_gxy*u_yl + metric->adm_gxz*u_zl;
    utcon2 = metric->adm_gxy*u_xl + metric->adm_gyy*u_yl + metric->adm_gyz*u_zl;
    utcon3 = metric->adm_gxz*u_xl + metric->adm_gyz*u_yl + metric->adm_gzz*u_zl;
    if (check==1) {
//TODO: error checking
//      CCTK_VInfo(CCTK_THORNSTRING,"Font fix failed!");
//          CCTK_VInfo(CCTK_THORNSTRING,"i,j,k = %d %d %d, stats.failure_checker = %d x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",i,j,k,stats.failure_checker,X[index],Y[index],Z[index],index,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6]);
      return check;
    }
  /*************************************************************/
  } else {
    printf("Font fix is only implemented for Hybrid EOS! The eos_type is set to %d.",eos->eos_type); //TODO: make a better error message
    return 1;
  } //EOS if

  diagnostics->failure_checker+=10000;
  diagnostics->font_fixes++;
  //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
  limit_velocity_and_convert_utilde_to_v(eos, metric, &utcon1, &utcon2, &utcon3, cons_undens->rho, prims_guess, diagnostics);
  double gijuiuj = metric->adm_gxx*SQR(utcon1) +
    2.0*metric->adm_gxy*utcon1*utcon2 + 2.0*metric->adm_gxz*utcon1*utcon3 +
    metric->adm_gyy*SQR(utcon2) + 2.0*metric->adm_gyz*utcon2*utcon3 +
    metric->adm_gzz*SQR(utcon3);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  u0L = (au0m1+1.0)*metric->lapseinv;
  prims_guess->rho = cons_undens->rho/(metric->lapse*u0L);

  int polytropic_index = find_polytropic_K_and_Gamma_index(eos,prims_guess->rho);
  double Gamma_ppoly_tab = eos->Gamma_ppoly_tab[polytropic_index];

  double uu_energy = prims_guess->press/(Gamma_ppoly_tab-1.0);

  if( eos->eos_type == 0 ) {
    prims_guess->press = pressure_rho0_u(eos, prims_guess->rho, uu_energy);
  }

  return 0;
}
