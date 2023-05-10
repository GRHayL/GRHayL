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


#include "IGH.h"
#include "Symmetry.h"

void GRHayL_IGH_conserv_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGH_conserv_to_prims;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr=0;
    ierr+=CartSymGN(cctkGH,"GRHayL_IGH::grmhd_conservatives");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr+=CartSymGN(cctkGH,"lapse::lapse_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_AH");
    ierr+=CartSymGN(cctkGH,"shift::shift_vars");
    if(ierr!=0) CCTK_VERROR("GRHayL_IGH ERROR (grep for it, foo!)  :(");
  }

  //Start the timer, so we can benchmark the primitives solver during evolution.
  //  Slower solver -> harder to find roots -> things may be going crazy!
  //FIXME: Replace this timing benchmark with something more meaningful, like the avg # of Newton-Raphson iterations per gridpoint!
  /*
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);
  */

  const double poison = 0.0/0.0;

  // Diagnostic variables.
  int failures=0;
  int font_fixes=0;
  int vel_limited_ptcount=0;
  int atm_resets=0;
  int rho_star_fix_applied=0;
  int pointcount=0;
  int failures_inhoriz=0;
  int pointcount_inhoriz=0;
  int backup0=0;
  int backup1=0;
  int backup2=0;
  int nan_found=0;
  double error_int_numer=0;
  double error_int_denom=0;
  int n_iter=0;
  double dummy1, dummy2, dummy3;
  double dummy4, dummy5, dummy6;

#pragma omp parallel for reduction(+:failures,font_fixes,vel_limited_ptcount,atm_resets,rho_star_fix_applied,pointcount,failures_inhoriz,pointcount_inhoriz,backup0,backup1,backup2,nan_found,error_int_numer,error_int_denom,n_iter) schedule(static)
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    //    const int index0 = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
    //    const int index1 = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
    //    const int index2 = CCTK_GFINDEX4D(cctkGH,i,j,k,2);

    con2prim_diagnostics diagnostics;
    initialize_diagnostics(&diagnostics);
//    diagnostics.c2p_fail_flag = con2prim_failed_flag[index]; from Leo's IGH

    // Read in BSSN metric quantities from gridfunctions and
    // set auxiliary and ADM metric quantities
    metric_quantities metric;
    GRHayL_enforce_detgtij_and_initialize_metric(
          alp[index],
          gxx[index], gxy[index], gxz[index],
          gyy[index], gyz[index], gzz[index],
          betax[index], betay[index], betaz[index],
          &metric);

    // Read in primitive variables from gridfunctions
    primitive_quantities prims;
    initialize_primitives(rho_b[index],
          pressure[index], eps[index],
          vx[index], vy[index], vz[index],
          0.0, 0.0, 0.0,
          poison, poison, poison, &prims);
          //wont have storage for these vars for hybrid
          //entropy[index], Y_e[index], temperature[index], &prims);

    // Read in conservative variables from gridfunctions
    conservative_quantities cons, cons_orig;
    initialize_conservatives(
          rho_star[index], tau[index],
          Stildex[index], Stildey[index], Stildez[index],
          poison, poison, &cons);
          //Y_e[index], entropy[index],

    // Here we save the original values of conservative variables in cons_orig for debugging purposes.
    cons_orig = cons;

    //FIXME: might slow down the code. Was formerly a CCTK_WARN
    if(isnan(cons.rho*cons.S_x*cons.S_y*cons.S_z*cons.tau*prims.Bx*prims.By*prims.Bz)) {
      CCTK_VINFO("NaN found at start of C2P kernel: index %d %d %d, st_i = %e %e %e, rho_* = %e, ~tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e\n",
                 i,j,k, cons.S_x,cons.S_y,cons.S_z,cons.rho,cons.tau,prims.Bx,prims.By,prims.Bz,
                 metric.adm_gxx,metric.adm_gxy,metric.adm_gxz,metric.adm_gyy,metric.adm_gyy,metric.adm_gzz,metric.psi6);
      diagnostics.nan_found++;
    }

    /************* Main conservative-to-primitive logic ************/
    int check=0;
    if(cons.rho>0.0) {
      // Apply the tau floor
      if( grhayl_eos->eos_type == grhayl_eos_hybrid )
        apply_inequality_fixes(grhayl_params, grhayl_eos, &metric, &prims, &cons, &diagnostics);

      // declare some variables for the C2P routine.
      conservative_quantities cons_undens;

      // Set the conserved variables required by the con2prim routine
      undensitize_conservatives(&metric, &cons, &cons_undens);

      /************* Conservative-to-primitive recovery ************/
      int check = grhayl_con2prim_multi_method(grhayl_params, grhayl_eos, &metric, &cons_undens, &prims, &diagnostics);

      if(check==0) {
        //Check for NAN!
        if( isnan(prims.rho*prims.press*prims.eps*prims.vx*prims.vy*prims.vz) ) {
          CCTK_VINFO("***********************************************************\n");
          CCTK_VINFO("NAN found in function %s (file: %s)\n",__func__,__FILE__);
          CCTK_VINFO("Input conserved variables:\n");
          CCTK_VINFO("rho_*, ~tau, ~S_{i}: %e %e %e %e %e\n", cons.rho, cons.tau, cons.S_x, cons.S_y, cons.S_z);
          CCTK_VINFO("Undensitized conserved variables:\n");
          CCTK_VINFO("D, tau, S_{i}: %e %e %e %e %e\n", cons_undens.rho, cons_undens.tau, cons_undens.S_x, cons_undens.S_y, cons_undens.S_z);
          CCTK_VINFO("Output primitive variables:\n");
          CCTK_VINFO("rho, P: %e %e\n", prims.rho, prims.press);
          CCTK_VINFO("v: %e %e %e\n", prims.vx, prims.vy, prims.vz);
          CCTK_VINFO("***********************************************************");
        }
      } else {
        CCTK_VINFO("Con2Prim and Font fix failed!");
        CCTK_VINFO("diagnostics->failure_checker = %d st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                diagnostics.failure_checker, cons_orig.S_x, cons_orig.S_y, cons_orig.S_z, cons_orig.rho, prims.Bx, prims.By, prims.Bz,
                metric.adm_gxx, metric.adm_gxy, metric.adm_gxz, metric.adm_gyy, metric.adm_gyz, metric.adm_gzz, metric.psi6);
      }
    } else {
      diagnostics.failure_checker+=1;
      reset_prims_to_atmosphere(grhayl_eos, &prims);
      rho_star_fix_applied++;
    } // if rho_star>0
    /***************************************************************/

    if( check != 0 ) {
      //--------------------------------------------------
      //----------- Primitive recovery failed ------------
      //--------------------------------------------------
      // Sigh, reset to atmosphere
      reset_prims_to_atmosphere(grhayl_eos, &prims);
      diagnostics.failure_checker+=100000;
      atm_resets++;
      // Then flag this point as a "success"
      check = 0;
      CCTK_VINFO("Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e\n",
                 cons_orig.tau, cons_orig.rho, cons_orig.S_x, cons_orig.S_y, cons_orig.S_z, cons_orig.rho/metric.psi6, grhayl_eos->rho_atm,
                 prims.Bx, prims.By, prims.Bz, metric.adm_gxx, metric.adm_gxy, metric.adm_gxz, metric.adm_gyy, metric.adm_gyy, metric.adm_gzz, metric.lapse);
    }

    //--------------------------------------------------
    //---------- Primitive recovery succeeded ----------
    //--------------------------------------------------
    // Enforce limits on primitive variables and recompute conservatives.
    stress_energy Tmunu;
    enforce_primitive_limits_and_compute_u0(grhayl_params, grhayl_eos, &metric, &prims, &diagnostics.failure_checker);
    compute_conservs_and_Tmunu(&metric, &prims, &cons, &Tmunu);

    //Now we compute the difference between original & new conservatives, for diagnostic purposes:
    error_int_numer += fabs(cons.tau - cons_orig.tau) + fabs(cons.rho - cons_orig.rho) + fabs(cons.S_x - cons_orig.S_x)
                       + fabs(cons.S_y - cons_orig.S_y) + fabs(cons.S_z - cons_orig.S_z);
    error_int_denom += cons_orig.tau + cons_orig.rho + fabs(cons_orig.S_x) + fabs(cons_orig.S_y) + fabs(cons_orig.S_z);

    if(check!=0) {
      diagnostics.failures++;
      if(metric.psi6 > grhayl_params->psi6threshold) {
        failures_inhoriz++;
        pointcount_inhoriz++;
      }
    }

    failure_checker[index] = diagnostics.failure_checker;

    return_primitives(&prims,
          &rho_b[index], &pressure[index], &eps[index],
          &vx[index], &vy[index], &vz[index],
          //&Bvec[index0], &Bvec[index1], &Bvec[index2],
          &dummy4, &dummy5, &dummy6,
          &dummy1, &dummy2, &dummy3);
          //wont have storage for these vars for hybrid
          //&entropy[index], &Y_e[index], &temperature[index]);

    return_conservatives(&cons,
          &rho_star[index], &tau[index],
          &Stildex[index], &Stildey[index], &Stildez[index],
          &dummy1, &dummy2);
          //&Y_e[index], &entropy[index]);

    if(update_Tmunu) {
      return_stress_energy(&Tmunu, &eTtt[index], &eTtx[index],
            &eTty[index], &eTtz[index], &eTxx[index],
            &eTxy[index], &eTxz[index], &eTyy[index],
            &eTyz[index], &eTzz[index]);
    }

    pointcount++;
    failures += diagnostics.failures;
    font_fixes += diagnostics.font_fix;
    vel_limited_ptcount += diagnostics.speed_limited;
    backup0 += diagnostics.backup[0];
    backup1 += diagnostics.backup[1];
    backup2 += diagnostics.backup[2];
    nan_found += diagnostics.nan_found;
    n_iter += diagnostics.n_iter;
}

  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VINFO("C2P: Lev: %d NumPts= %d | Fixes: Font= %d VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e | %.2f iters/gridpt",
               (int)GetRefinementLevel(cctkGH), pointcount,
               font_fixes, vel_limited_ptcount, rho_star_fix_applied,
               failures, failures_inhoriz, pointcount_inhoriz,
               error_int_numer/error_int_denom, error_int_denom,
               (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ));
  }
}
