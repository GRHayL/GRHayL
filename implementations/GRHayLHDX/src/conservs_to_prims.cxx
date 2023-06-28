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

#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_conservs_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_conservs_to_prims;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL poison = 0.0/0.0;

  // Diagnostic variables.
  int failures=0;
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
  CCTK_REAL error_int_numer=0;
  CCTK_REAL error_int_denom=0;
  int n_iter=0;

  constexpr std::array<int, Loop::dim> indextype = {1, 1, 1};
  const Loop::GF3D2layout layout(cctkGH, indextype);

//I don't think I can do this sort of diagnostic summing in carpetx
//#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,atm_resets,rho_star_fix_applied,pointcount,failures_inhoriz,pointcount_inhoriz,backup0,backup1,backup2,nan_found,error_int_numer,error_int_denom,n_iter) schedule(static)
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::GF3D2index index(layout, p.I);

    con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);

    // Read in ADM metric quantities from gridfunctions and
    // set auxiliary and ADM metric quantities
    metric_quantities ADM_metric;
    ghl_initialize_metric(
          ccc_lapse(index),
          ccc_betax(index), ccc_betay(index), ccc_betaz(index),
          ccc_gxx(index), ccc_gxy(index), ccc_gxz(index),
          ccc_gyy(index), ccc_gyz(index), ccc_gzz(index),
          &ADM_metric);

    ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    // Read in primitive variables from gridfunctions
    primitive_quantities prims;
    ghl_initialize_primitives(
          rho_b(index), pressure(index), eps(index),
          vx(index), vy(index), vz(index),
          0.0, 0.0, 0.0,
          poison, poison, poison, &prims);

    // Read in conservative variables from gridfunctions
    conservative_quantities cons, cons_orig;
    ghl_initialize_conservatives(
          rho_star(index), tau(index),
          Stildex(index), Stildey(index), Stildez(index),
          poison, poison, &cons);

    // Here we save the original values of conservative variables in cons_orig for debugging purposes.
    cons_orig = cons;

    //FIXME: might slow down the code. Was formerly a CCTK_WARN
    if(isnan(cons.rho*cons.tau*cons.SD[0]*cons.SD[1]*cons.SD[2])) {
      CCTK_VERROR("NaN found at start of C2P kernel: index %d %d %d, rho_* = %e, ~tau = %e, ~S_i = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e\n",
                 p.i, p.j, p.k, cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2],
                 ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2], ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2], ADM_metric.sqrt_detgamma);
    }

    /************* Main conservative-to-primitive logic ************/
    int check=0;
    if(cons.rho>0.0) {
      // Apply the tau floor
      if( ghl_eos->eos_type == ghl_eos_hybrid )
        ghl_apply_conservative_limits(
              ghl_params, ghl_eos, &ADM_metric,
              &prims, &cons, &diagnostics);

      // declare some variables for the C2P routine.
      conservative_quantities cons_undens;

      // Set the conserved variables required by the con2prim routine
      ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

      /************* Conservative-to-primitive recovery ************/
      int check = ghl_con2prim_multi_method(
            ghl_params, ghl_eos, &ADM_metric, &metric_aux,
            &cons_undens, &prims, &diagnostics);

      if(check==0) {
        //Check for NAN!
        if( isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]) ) {
          CCTK_VERROR("***********************************************************\n"
                      "NAN found in function %s (file: %s)\n"
                      "Input conserved variables:\n"
                      "rho_*, ~tau, ~S_{i}: %e %e %e %e %e\n"
                      "Undensitized conserved variables:\n"
                      "D, tau, S_{i}: %e %e %e %e %e\n"
                      "Output primitive variables:\n"
                      "rho, P: %e %e\n"
                      "v: %e %e %e\n"
                      "***********************************************************",
                      __func__,__FILE__,
                      cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2],
                      cons_undens.rho, cons_undens.tau, cons_undens.SD[0], cons_undens.SD[1], cons_undens.SD[2],
                      prims.rho, prims.press,
                      prims.vU[0], prims.vU[1], prims.vU[2]);
        }
      } else {
        CCTK_VINFO("Con2Prim failed!");
        CCTK_VINFO("diagnostics->failure_checker = %d rho_* = %e, ~tau = %e, ~S_i = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                diagnostics.failure_checker, cons_orig.rho, cons_orig.tau, cons_orig.SD[0], cons_orig.SD[1], cons_orig.SD[2],
                ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2], ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2], ADM_metric.sqrt_detgamma);
      }
    } else {
      diagnostics.failure_checker+=1;
      ghl_set_prims_to_constant_atm(ghl_eos, &prims);
      //rho_star_fix_applied++;
    } // if rho_star>0
    /***************************************************************/

    if( check != 0 ) {
      //--------------------------------------------------
      //----------- Primitive recovery failed ------------
      //--------------------------------------------------
      // Sigh, reset to atmosphere
      ghl_set_prims_to_constant_atm(ghl_eos, &prims);
      diagnostics.failure_checker+=100000;
      //atm_resets++;
      // Then flag this point as a "success"
      check = 0;
      CCTK_VINFO("Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, gij=%e %e %e %e %e %e, alpha=%e\n",
                 cons_orig.rho, cons_orig.tau, cons_orig.SD[0], cons_orig.SD[1], cons_orig.SD[2], cons_orig.rho/ADM_metric.sqrt_detgamma, ghl_eos->rho_atm,
                 ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2], ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2], ADM_metric.lapse);
    }

    //--------------------------------------------------
    //---------- Primitive recovery succeeded ----------
    //--------------------------------------------------
    // Enforce limits on primitive variables and recompute conservatives.
    ghl_enforce_primitive_limits_and_compute_u0(
          ghl_params, ghl_eos, &ADM_metric,
          &prims, &diagnostics.failure_checker);
    ghl_compute_conservs(
          &ADM_metric, &metric_aux, &prims, &cons);

    CCTK_REAL dummy1, dummy2, dummy3;
    CCTK_REAL dummy4, dummy5, dummy6;
    ghl_return_primitives(
          &prims,
          &rho_b(index), &pressure(index), &eps(index),
          &vx(index), &vy(index), &vz(index),
          &dummy4, &dummy5, &dummy6,
          &dummy1, &dummy2, &dummy3);
    u0(index) = prims.u0;

    ghl_return_conservatives(
          &cons,
          &rho_star(index), &tau(index),
          &Stildex(index), &Stildey(index), &Stildez(index),
          &dummy1, &dummy2);

    // We tally up our diagnostic information
    //if(check!=0) {
    //  diagnostics.failures++;
    //  if(ADM_metric.sqrt_detgamma > ghl_params->psi6threshold) {
    //    failures_inhoriz++;
    //    pointcount_inhoriz++;
    //  }
    //}

    failure_checker(index) = diagnostics.failure_checker;

    //We compute the difference between original & new conservatives, for diagnostic purposes:
    //error_int_numer += fabs(cons.tau - cons_orig.tau) + fabs(cons.rho - cons_orig.rho) + fabs(cons.SD[0] - cons_orig.SD[0])
    //                   + fabs(cons.SD[1] - cons_orig.SD[1]) + fabs(cons.SD[2] - cons_orig.SD[2]);
    //error_int_denom += cons_orig.tau + cons_orig.rho + fabs(cons_orig.SD[0]) + fabs(cons_orig.SD[1]) + fabs(cons_orig.SD[2]);

    //pointcount++;
    //failures += diagnostics.failures;
    //vel_limited_ptcount += diagnostics.speed_limited;
    //backup0 += diagnostics.backup[0];
    //backup1 += diagnostics.backup[1];
    //backup2 += diagnostics.backup[2];
    //nan_found += diagnostics.nan_found;
    //n_iter += diagnostics.n_iter;
  }); // ccc loop everywhere

  //if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
  //  CCTK_VINFO("C2P: Lev: %d NumPts= %d | Fixes: Font= %d VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e | %.2f iters/gridpt",
  //             (int)GetRefinementLevel(cctkGH), pointcount,
  //             backup0, vel_limited_ptcount, rho_star_fix_applied,
  //             failures, failures_inhoriz, pointcount_inhoriz,
  //             error_int_numer/error_int_denom, error_int_denom,
  //             (CCTK_REAL)n_iter/( (CCTK_REAL)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ));
  //}
}
