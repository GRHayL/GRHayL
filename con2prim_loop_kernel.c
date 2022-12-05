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
#include "con2prim_header.h"

void con2prim_loop_kernel(const GRMHD_parameters *restrict params, const eos_parameters *restrict eos,
                          metric_quantities *restrict metric, conservative_quantities *restrict cons,
                          primitive_quantities *restrict prims, con2prim_diagnostics *restrict diagnostics, stress_energy *restrict Tmunu) {

  //FIXME: might slow down the code. Was formerly a CCTK_WARN
  if(isnan(cons->rho*cons->S_x*cons->S_y*cons->S_z*cons->tau*prims->Bx*prims->By*prims->Bz)) {
    printf("NaN found at start of C2P kernel: st_i = %e %e %e, rho_* = %e, ~tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e\n",
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
      apply_inequality_fixes(params, eos, metric, prims, cons, diagnostics);

/*************************************************************************/
    // declare some variables for the C2P routine.
    conservative_quantities cons_undens;
    primitive_quantities prims_guess;

    // Set the conserved variables required by the con2prim routine
    undensitize_conservatives(metric, cons, &cons_undens);
  
    /************* Conservative-to-primitive recovery ************/

    for(int which_guess=1;which_guess<3;which_guess++) {

      // Set primitive guesses
      guess_primitives(eos, which_guess, metric, prims, cons, &prims_guess);
      int check = C2P_Select_Hybrid_Method(params, eos, params->main_routine, metric, &cons_undens, &prims_guess, diagnostics);

      if( (check != 0) && (params->backup_routine[0] != None) ) {
        // Backup 1 triggered
        diagnostics->backup[0] = 1;
        // Recompute guesses
        guess_primitives(eos, which_guess, metric, prims, cons, &prims_guess);
        // Backup routine #1
        check = C2P_Select_Hybrid_Method(params, eos, params->backup_routine[0], metric, &cons_undens, &prims_guess, diagnostics);

        if( (check != 0) && (params->backup_routine[1] != None) ) {
          // Backup 2 triggered
          diagnostics->backup[1] = 1;
          // Recompute guesses
          guess_primitives(eos, which_guess, metric, prims, cons, &prims_guess);
          // Backup routine #2
          check = C2P_Select_Hybrid_Method(params, eos, params->backup_routine[1], metric, &cons_undens, &prims_guess, diagnostics);

          if( (check != 0) && (params->backup_routine[2] != None) ) {
            // Backup 3 triggered
            diagnostics->backup[2] = 1;
            // Recompute guesses
            guess_primitives(eos, which_guess, metric, prims, cons, &prims_guess);
            // Backup routine #3
            check = C2P_Select_Hybrid_Method(params, eos, params->backup_routine[2], metric, &cons_undens, &prims_guess, diagnostics);
          }
        }
      }
      /*************************************************************/
  
      if(check!=0)
        check = font_fix(eos, metric, cons, prims, &prims_guess, diagnostics);
  
      if(check==0) {
  //       Check for NAN!
        if( isnan(prims_guess.rho*prims_guess.press*prims_guess.eps*prims_guess.vx*prims_guess.vy*prims_guess.vz) ) {
          printf("***********************************************************\n");
          printf("NAN found in function %s (file: %s)\n",__func__,__FILE__);
          printf("Input conserved variables:\n");
          printf("rho_*, ~tau, ~S_{i}: %e %e %e %e %e\n", cons->rho, cons->tau, cons->S_x, cons->S_y, cons->S_z);
          printf("Undensitized conserved variables:\n");
          printf("D, tau, S_{i}: %e %e %e %e %e\n", cons_undens.rho, cons_undens.tau, cons_undens.S_x, cons_undens.S_y, cons_undens.S_z);
          printf("Output primitive variables:\n");
          printf("rho, P: %e %e\n", prims_guess.rho, prims_guess.press);
          printf("v: %e %e %e\n", prims_guess.vx, prims_guess.vy, prims_guess.vz);
          printf("***********************************************************");
        }
  
        *prims = prims_guess;
        which_guess=3; //TODO: can we make the multiple guess loop cleaner?
      } else {
        printf("Con2Prim and Font fix failed!");
        printf("diagnostics->failure_checker = %d st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                diagnostics->failure_checker, cons_orig.S_x, cons_orig.S_y, cons_orig.S_z, cons_orig.rho, prims->Bx, prims->By, prims->Bz,
                metric->adm_gxx, metric->adm_gxy, metric->adm_gxz, metric->adm_gyy, metric->adm_gyz, metric->adm_gzz, metric->psi6);
      }
    } //If we didn't find a root, then try again with a different guess.

  } else {
    diagnostics->failure_checker+=1;
    reset_prims_to_atmosphere(eos, prims, diagnostics);
    diagnostics->rho_star_fix_applied++;
  } // if rho_star>0

  if( check != 0 ) {
    //--------------------------------------------------
    //----------- Primitive recovery failed ------------
    //--------------------------------------------------
    // Sigh, reset to atmosphere
    reset_prims_to_atmosphere( eos, prims, diagnostics);
    diagnostics->failure_checker+=100000;
    diagnostics->atm_resets++;
    // Then flag this point as a "success"
    check = 0;
//TODO: change to prinf
    printf("Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e\n",
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
//printf("Cons: tau rho S %.16e %.16e %.16e %.16e %.16e\n", cons->tau, cons->rho, cons->S_x, cons->S_y, cons->S_z);
//printf("Orig: tau rho S %.16e %.16e %.16e %.16e %.16e\n", cons_orig.tau, cons_orig.rho, cons_orig.S_x, cons_orig.S_y, cons_orig.S_z);
  diagnostics->error_int_numer += fabs(cons->tau - cons_orig.tau) + fabs(cons->rho - cons_orig.rho) + fabs(cons->S_x - cons_orig.S_x)
                    + fabs(cons->S_y - cons_orig.S_y) + fabs(cons->S_z - cons_orig.S_z);
  diagnostics->error_int_denom += cons_orig.tau + cons_orig.rho + fabs(cons_orig.S_x) + fabs(cons_orig.S_y) + fabs(cons_orig.S_z);

  if(check!=0) {
    diagnostics->failures++;
    if(metric->psi6 > params->psi6threshold) {
      diagnostics->failures_inhoriz++;
      diagnostics->pointcount_inhoriz++;
    }
  }
  diagnostics->pointcount++;
}
