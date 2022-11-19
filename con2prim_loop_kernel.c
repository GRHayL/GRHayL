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
                          primitive_quantities *restrict prims, con2prim_diagnostics *restrict diagnostics, stress_energy *restrict Tmunu) {

  //FIXME: might slow down the code.
  if(isnan(cons->rho*cons->S_x*cons->S_y*cons->S_z*cons->tau*prims->Bx*prims->By*prims->Bz)) {
    CCTK_VWARN(CCTK_WARN_ALERT,"NaN found at start of C2P kernel: st_i = %e %e %e, rho_* = %e, ~tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
               cons->S_x,cons->S_y,cons->S_z,cons->rho,cons->tau,prims->Bx,prims->By,prims->Bz,
               metric->adm_gxx,metric->adm_gxy,metric->adm_gxz,metric->adm_gyy,metric->adm_gyy,metric->adm_gzz,metric->psi6);
    diagnostics->nan_found++;
  }

  // Here we save the original values of conservative variables in cons_orig for debugging purposes.
  conservative_quantities cons_orig = *cons;

  int check=0;
  if(cons->rho>0.0) {
    // Apply the tau floor
    if( eos->eos_type == 0 )
      apply_tau_floor(params, eos, metric, prims, cons, diagnostics);

/*************************************************************************/
    // declare some variables for the C2P routine.
    conservative_quantities cons_undens;
    primitive_quantities prims_guess;

    // Set the conserved variables required by the con2prim routine
    undensitize_conservatives( eos, params->main_routine, metric, prims, cons, &cons_undens );
  
    /************* Conservative-to-primitive recovery ************/

    for(int which_guess=1;which_guess<3;which_guess++) {

      // Set primitive guesses
      guess_primitives( eos, params->main_routine, which_guess, metric, prims, cons, &prims_guess );
      int check = C2P_Select_Hybrid_Method(eos, params->main_routine, metric, &cons_undens, &prims_guess, diagnostics);

      if( (check != 0) && (params->backup_routine[0] != None) ) {
        // Backup 1 triggered
        diagnostics->backup[0] = 1;
        // Recompute guesses
        guess_primitives( eos,params->backup_routine[0], which_guess, metric, prims, cons, &prims_guess );
        // Backup routine #1
        check = C2P_Select_Hybrid_Method(eos, params->backup_routine[0], metric, &cons_undens, &prims_guess, diagnostics);

        if( (check != 0) && (params->backup_routine[1] != None) ) {
          // Backup 2 triggered
          diagnostics->backup[1] = 1;
          // Recompute guesses
          guess_primitives( eos,params->backup_routine[1], which_guess, metric, prims, cons, &prims_guess );
          // Backup routine #2
          check = C2P_Select_Hybrid_Method(eos, params->backup_routine[1], metric, &cons_undens, &prims_guess, diagnostics);

          if( (check != 0) && (params->backup_routine[2] != None) ) {
            // Backup 3 triggered
            diagnostics->backup[2] = 1;
            // Recompute guesses
            guess_primitives( eos,params->backup_routine[2], which_guess, metric, prims, cons, &prims_guess );
            // Backup routine #3
            check = C2P_Select_Hybrid_Method(eos, params->backup_routine[2], metric, &cons_undens, &prims_guess, diagnostics);
          }
        }
      }
      /*************************************************************/
  
      if(check!=0) {
        check = font_fix(eos, metric, cons, prims, &prims_guess, diagnostics);
//if(fabs(cons_orig.tau - 3.8712996396e-09) < 1.0e-18 ) {
//  CCTK_VINFO("Font Fix final prims: %.16e %.16e\n  %.16e %.16e %.16e", prims_guess.rho, prims_guess.press, prims_guess.vx, prims_guess.vy, prims_guess.vz);
//}
        diagnostics->font_fixes++;
      }
  
      if(check==0) {
  //       Check for NAN!
        if( isnan(prims_guess.rho*prims_guess.press*prims_guess.eps*prims_guess.vx*prims_guess.vy*prims_guess.vz) ) {
          CCTK_VINFO("***********************************************************");
          CCTK_VINFO("NAN found in function %s (file: %s)",__func__,__FILE__);
          CCTK_VINFO("Input conserved variables:");
          CCTK_VINFO("rho_*, ~tau, ~S_{i}: %e %e %e %e %e", cons->rho, cons->tau, cons->S_x, cons->S_y, cons->S_z);
          CCTK_VINFO("Undensitized conserved variables:");
          CCTK_VINFO("D, tau, S_{i}: %e %e %e %e %e", cons_undens.rho, cons_undens.tau, cons_undens.S_x, cons_undens.S_y, cons_undens.S_z);
          CCTK_VINFO("Output primitive variables:");
          CCTK_VINFO("rho, P: %e %e", prims_guess.rho, prims_guess.press);
          CCTK_VINFO("v: %e %e %e", prims_guess.vx, prims_guess.vy, prims_guess.vz);
          CCTK_VINFO("***********************************************************");
        }
  
        *prims = prims_guess;
  //CCTK_VINFO("cons: rho=%.16e, ~tau=%.16e, ~S=(%.16e, %.16e, %.16e),", cons->rho,cons->tau,cons->S_x,cons->S_y,cons->S_z);
  //CCTK_VINFO("prims: rho=%.16e, press=%.16e, vx=%.16e, vy=%.16e, vz=%.16e",prims->rho,prims->press,prims->vx,prims->vy,prims->vz);
  //CCTK_VINFO("      B=(%.16e, %.16e, %.16e)",prims->Bx,prims->By,prims->Bz);
        which_guess=3;
      } //If we didn't find a root, then try again with a different guess.
    }
    diagnostics->failure_checker+=100000;

/************************************************************************/
//    check = con2prim(params, eos, metric, cons, prims, diagnostics);
  } else {
    diagnostics->failure_checker+=1;
    reset_prims_to_atmosphere(eos, prims, diagnostics);
    diagnostics->rho_star_fix_applied++;
  }

  if( check != 0 ) {
    //--------------------------------------------------
    //----------- Primitive recovery failed ------------
    //--------------------------------------------------
CCTK_VINFO("C2P and FF failures! Reseting to atm...");
    // Sigh, reset to atmosphere
    reset_prims_to_atmosphere( eos, prims, diagnostics);
    diagnostics->atm_resets++;
    // Then flag this point as a "success"
    check = 0;
//TODO: change to prinf
    CCTK_VINFO("Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
               cons_orig.tau,cons_orig.rho,cons_orig.S_x,cons_orig.S_y,cons_orig.S_z,cons_orig.rho/metric->psi6,eos->rho_atm,
               prims->Bx,prims->By,prims->Bz,metric->adm_gxx,metric->adm_gxy,metric->adm_gxz,metric->adm_gyy,metric->adm_gyy,metric->adm_gzz,metric->lapse);
  }

  //--------------------------------------------------
  //---------- Primitive recovery succeeded ----------
  //--------------------------------------------------
  // Enforce limits on primitive variables and recompute conservatives.
  double TUPMUNU[10],TDNMUNU[10];
  enforce_limits_on_primitives_and_recompute_conservs(params, eos, metric, prims, cons, TUPMUNU, TDNMUNU, Tmunu, diagnostics);

  //Now we compute the difference between original & new conservatives, for diagnostic purposes:
//CCTK_VINFO("Cons: tau rho S %.16e %.16e %.16e %.16e %.16e", cons->tau, cons->rho, cons->S_x, cons->S_y, cons->S_z);
//CCTK_VINFO("Orig: tau rho S %.16e %.16e %.16e %.16e %.16e", cons_orig.tau, cons_orig.rho, cons_orig.S_x, cons_orig.S_y, cons_orig.S_z);
  diagnostics->error_int_numer += fabs(cons->tau - cons_orig.tau) + fabs(cons->rho - cons_orig.rho) + fabs(cons->S_x - cons_orig.S_x)
                    + fabs(cons->S_y - cons_orig.S_y) + fabs(cons->S_z - cons_orig.S_z);
  diagnostics->error_int_denom += cons_orig.tau + cons_orig.rho + fabs(cons_orig.S_x) + fabs(cons_orig.S_y) + fabs(cons_orig.S_z);

  if(check!=0) {
    diagnostics->failures++;
    if(exp(metric->bssn_phi*6.0)>params->psi6threshold) {
      diagnostics->failures_inhoriz++;
      diagnostics->pointcount_inhoriz++;
    }
  }
  diagnostics->pointcount++;
}
