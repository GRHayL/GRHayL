/* We evolve forward in time a set of functions called the
 *  "conservative variables", and any time the conserv's
 *  are updated, we must solve for the primitive variables
 *  (rho, pressure, velocities) using a Newton-Raphson
 *  technique, before reconstructing & evaluating the RHSs
 *  of the MHD equations again.
 *
 * This file contains the driver routine for this Newton-
 *  Raphson solver. Truncation errors in conservative
 *  variables can lead to no physical solutions in
 *  primitive variables. We correct for these errors here
 *  through a number of tricks described in the appendices
 *  of http://arxiv.org/pdf/1112.0568.pdf.
 *
 * This is a wrapper for the 2d solver of Noble et al. See
 *  harm_utoprim_2d.c for references and copyright notice
 *  for that solver. This wrapper was primarily written by
 *  Zachariah Etienne & Yuk Tung Liu, in 2011-2013.
 *
 * For optimal compatibility, this wrapper is licensed under
 *  the GPL v2 or any later version.
 *
 * Note that this code assumes a simple gamma law for the
 *  moment, though it would be easy to extend to a piecewise
 *  polytrope. */

// Standard #include's
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "Symmetry.h"
#include "con2prim_header.h"

void con2prim_loop_kernel(const GRMHD_parameters *restrict params, const eos_parameters *restrict eos,
                          metric_quantities *restrict metric, conservative_quantities *restrict cons,
                          primitive_quantities *restrict prims, con2prim_diagnostics *restrict diagnostics) {


  // Only attempt a primitive recovery if this is the first
  // attempt at doing so or if the previous attempt failed.
//  if( (loop_count == 0) || (diagnostics->c2p_fail_flag != 0) ) {
//TODO: if is in surrounding code rn

//  //FIXME: might slow down the code.
//  if(robust_isnan(cons->rho*cons->S_x*cons->S_y*cons->S_z*cons->tau*prims->Bx*prims->By*prims->Bz)) {
//    CCTK_VWARN(CCTK_WARN_ALERT,"NAN FOUND: index = %d, x,y,z = %e %e %e , st_i = %e %e %e, rhostar = %e, tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
//               i,j,k,x[i],y[i],z[i],index,
//               cons->S_x,cons->S_y,cons->S_z,cons->rho,cons->tau,prims->Bx,prims->By,prims->Bz,
//               metric->adm_gxx,metric->adm_gxy,metric->adm_gxz,metric->adm_gyy,metric->adm_gyy,metric->adm_gzz,metric->psi6);
//    diagnostics->nan_found++;
//  }

  // Here we save the original values of conservative variables in cons_orig for debugging purposes.
  struct conservative_quantities cons_orig;
  cons_orig = *cons;

  int check=0;
  diagnostics->which_routine  = None;
  if(cons->rho>0.0) {
    // Apply the tau floor
    if( eos->eos_type == 0 )
      apply_tau_floor(params, eos, metric, prims, cons, diagnostics);

    for(int ii=0;ii<3;ii++) {
      check = con2prim(params, eos, metric, cons, prims, diagnostics);
      if(check==0) ii=4;
      else diagnostics->failure_checker+=100000;
    }
  } else {
    diagnostics->failure_checker+=1;
    reset_prims_to_atmosphere(eos, prims);
    diagnostics->rho_star_fix_applied++;
  }

  if( check != 0 ) {
    //--------------------------------------------------
    //----------- Primitive recovery failed ------------
    //--------------------------------------------------
    // Increment the failure flag
    diagnostics->c2p_fail_flag += 1;
    if( diagnostics->c2p_fail_flag > 4 ) {
      // Sigh, reset to atmosphere
      reset_prims_to_atmosphere( eos, prims );
      diagnostics->atm_resets++;
      // Then flag this point as a "success"
      check = 0;
      diagnostics->c2p_fail_flag = 0;
//TODO: change to prinf
//      if( eos->type == 0 ) {
//        CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
//                   cons_orig.tau,cons_orig.rho,cons_orig.S_x,cons_orig.S_y,cons_orig.S_z,cons_orig.rho/metric->psi6,eos->rho_atm,prims->Bx,prims->By,prims->Bz,metric->adm_gxx,metric->adm_gxy,metric->adm_gxz,metric->adm_gyy,metric->adm_gyy,metric->adm_gzz,metric->lapse);
//      }
//      else if( eos->type == 1 ) {
//        CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: %e %e %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
//                   cons_orig.tau,cons_orig.rho,cons_orig.S_x,cons_orig.S_y,cons_orig.S_z,cons_orig.Y_e,cons_orig.entropy,cons_orig.rho/metric->psi6,eos->rho_atm,prims->Bx,prims->By,prims->Bz,metric->adm_gxx,metric->adm_gxy,metric->adm_gxz,metric->adm_gyy,metric->adm_gyy,metric->adm_gzz, metric->lapse);
//      }
    }
//    else {
//      // Increment the number of gridpoints which will need the average fix
//      num_of_conservative_averagings_needed++;
//    }
  }

  if( check == 0 ) {
    //--------------------------------------------------
    //---------- Primitive recovery succeeded ----------
    //--------------------------------------------------
    // Enforce limits on primitive variables and recompute conservatives.
    double TUPMUNU[10],TDNMUNU[10];
    enforce_limits_on_primitives_and_recompute_conservs(params, eos, metric, prims, cons, TUPMUNU, TDNMUNU, &diagnostics->failure_checker);

//    if(update_Tmunu) {
//      int ww=0;
//      eTtt[i] = TDNMUNU[ww++];
//      eTtx[i] = TDNMUNU[ww++];
//      eTty[i] = TDNMUNU[ww++];
//      eTtz[i] = TDNMUNU[ww++];
//      eTxx[i] = TDNMUNU[ww++];
//      eTxy[i] = TDNMUNU[ww++];
//      eTxz[i] = TDNMUNU[ww++];
//      eTyy[i] = TDNMUNU[ww++];
//      eTyz[i] = TDNMUNU[ww++];
//      eTzz[i] = TDNMUNU[ww  ];
//    }

    //Now we compute the difference between original & new conservatives, for diagnostic purposes:
    diagnostics->error_int_numer += fabs(cons->tau - cons_orig.tau) + fabs(cons->rho - cons_orig.rho) + fabs(cons->S_x - cons_orig.S_x)
                      + fabs(cons->S_y - cons_orig.S_y) + fabs(cons->S_z - cons_orig.S_z) + fabs(cons->Y_e - cons_orig.Y_e);
    diagnostics->error_int_denom += cons_orig.tau + cons_orig.rho + fabs(cons_orig.S_x) + fabs(cons_orig.S_y) + fabs(cons_orig.S_z)
                      + cons_orig.Y_e;

//TODO: Remove cctk function   if(stats.nan_found==1) { CCTK_VWARN(CCTK_WARN_ALERT,"Found NAN while imposing speed limit"); diagnostics->nan_found++; }
    if(check!=0) {
      diagnostics->failures++;
      if(exp(metric->bssn_phi*6.0)>params->psi6threshold) {
        diagnostics->failures_inhoriz++;
        diagnostics->pointcount_inhoriz++;
      }
    }
    diagnostics->pointcount++;
    /***************************************************************************************************************************/
  }
//  } // if( (loop_count == 0) || (c2p_fail_flag != 0) )
}
