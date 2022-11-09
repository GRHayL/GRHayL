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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "con2prim_header.h"

//#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER
//#include "standalone_conserv_to_prims_main_function.h"
//#else
//#include "cctk.h"
//#include "cctk_Arguments.h"
//#include "cctk_Parameters.h"
//#include "Symmetry.h"
//
////#include "IllinoisGRMHD_headers.h"
//#include "con2prim_headers.h"
//#include "new_header.h"
//#include "inlined_functions.h"
//#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"
//
//extern "C" void IllinoisGRMHD_conserv_to_prims(CCTK_ARGUMENTS) {
//  DECLARE_CCTK_ARGUMENTS;
//  DECLARE_CCTK_PARAMETERS;
//
//  // We use proper C++ here, for file I/O later.
//  using namespace std;
//#endif

extern "C" void con2prim_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // We use proper C++ here, for file I/O later.
  using namespace std;

// Should probably be in schedule function, not here.
#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER
  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr=0;
    ierr+=CartSymGN(cctkGH,"IllinoisGRMHD::grmhd_conservatives");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr+=CartSymGN(cctkGH,"lapse::lapse_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_AH");
    ierr+=CartSymGN(cctkGH,"shift::shift_vars");
    if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"IllinoisGRMHD ERROR (grep for it, foo!)  :(");
  }
#endif

  double poison = 1e200;
  main = Noble2D;
  backup_routines = {None,None,None};

  GRMHD_parameters params;
  initialize_parameters(&params, main, backup[3], false, false, true, Psi6threshold, update_Tmunu);

  int eos_type = 0; // Hybrid=0, Tabulated=1;
  eos_parameters eos;
  initialize_general_eos(&eos, 0, tau_atm, GAMMA_SPEED_LIMIT,
             poison, poison,poison, //epsilon
             poison, poison, poison, //pressure
             poison, poison, poison, //entropy
             rho_b_atm, rho_atm, rho_max);

  initialize_hybrid_eos( &eos, 
             neos, rho_ppoly_tab,
             Gamma_ppoly_tab[], K_ppoly_tab[],
             eps_integ_const[], gamma_th);

  con2prim_diagnostics diagnostics;
  initialize_diagnostics(&diagnostics);

  int imin=0,imax=cctk_lsh[0];
  int jmin=0,jmax=cctk_lsh[1];
  int kmin=0,kmax=cctk_lsh[2];
  int npoints = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  // Whenever we get a conservative-to-primitive major failure, i.e. all
  // the routines and backups failed to recover the primitives from the
  // input conservatives, we will introduce a new fix, in which we will
  // reset the conservative variables at the given point by a weighted
  // average of the conservative variables at the neighboring points.
  // After that, the con2prim attempt will be retried. This mask allows
  // us to flag points in which the averaging procedure must be performed.

  // Initialization of the masks above. Flag meaning:
  //
  // 0 -> Primitive recovery succeeded
  // 1 -> Primitive recovery failed
#pragma omp parallel for
  for(int i=0;i<npoints;i++)
    con2prim_failed_flag[i] = 0;

  // We now add an integer to count the number of
  // points in which the averaging fix is required.
  // We initialize it to a nonzero value so that
  // the while condition below is triggered at least
  // once.
//  int num_of_conservative_averagings_needed = 1;
//  int cons_avgs = 0;
//  int loop_count = 0;
//  int nan_found = 0;

//  while( num_of_conservative_averagings_needed > 0 ) {

    // Now set the number of conserved averages to zero
//    num_of_conservative_averagings_needed = 0;

#pragma omp parallel for reduction(+:diagnostics.failures,diagnostics.vel_limited_ptcount,diagnostics.font_fixes,diagnostics.pointcount,diagnostics.failures_inhoriz,diagnostics.pointcount_inhoriz,diagnostics.error_int_numer,diagnostics.error_int_denom,diagnostics.rho_star_fix_applied,diagnostics.atm_resets,diagnostics.backup[0],diagnostics.backup[1],diagnostics.backup[2],num_of_conservative_averagings_needed,diagnostics.nan_found) schedule(static)
for (int index=0;index<npoints;index++) { //Since everything is point-wise, just make the loop 1D for OMP
//  for(int k=kmin;k<kmax;k++) {
//    for(int j=jmin;j<jmax;j++) {
//      for(int i=imin;i<imax;i++) {
//        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        diagnostics.c2p_fail_flag = con2prim_failed_flag[index];

        // Read in BSSN metric quantities from gridfunctions and
        // set auxiliary and ADM metric quantities
        struct metric_quantities metric;
        initialize_metric(&metric, phi[index], psi[index], lapse[index],
                          gxx[index], gxy[index], gxz[index],
                          gyy[index], gyz[index], gzz[index],
                          betax[index], betay[index], betaz[index]);
           
        // Read in primitive variables from gridfunctions
        struct primitive_quantities prims;
        initialize_primitives(&eos, &metric, rho[index]/metric.psi6,
                              press[index], epsilon[index],
                              vx[index], vy[index], vz[index],
                              Bx[index], By[index], Bz[index],
                              entropy[index], Y_e[index], temp[index],
                              &prims);

        // Read in conservative variables from gridfunctions
        struct conservative_quantities cons;
        initialize_conservatives(&params, &eos, rho[index], tau[index],
                                 S_x[index], S_y[index], S_z[index],
                                 Y_e[index], entropy[index] &cons);

        con2prim_loop_kernel(&param, &eos, &metric, &cons, &prims, &diagnostics);

        failure_checker[index] = diagnostics.failure_checker;
        con2prim_failed_flag[index] = diagnostics.c2p_fail_flag;
        igm_c2p_mask[index] = diagnostics.which_routine;

        return_primitives(&eos, &prims, &rho[index], &press[index], &epsilon[index],
                      &vx[index], &vy[index], &vz[index], &Bx[index], &By[index], &Bz[index],
                      &entropy[index], &Y_e[index], &temp[index]);

        return_conservatives(&params, &eos, &cons, &rho[index], &tau[index],
                            &S_x[index], &S_y[index], S_z[index],
                            &Y_e[index],  &entropy[index]);

//      } // for(int i=imin;i<imax;i++)
//    } // for(int j=jmin;j<jmax;j++)
//  } // for(int k=kmin;k<kmax;k++)
  }

//    if( num_of_conservative_averagings_needed > cons_avgs ) cons_avgs = num_of_conservative_averagings_needed;

    // Increment the loop counter. This counter is used to avoid
    // attempting to recover primitives at points in which we
    // have already done so successfully.
//    loop_count++;

//  } // while( num_of_conservative_averagings_needed > 0 )

  // fclose(c2pmaskfile);

  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VInfo(CCTK_THORNSTRING,"C2P: Lev: %d NumPts= %d | Fixes: BU: %d %d %d Font= %d VL= %d rho*= %d AVG= %d ATM= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e",
               (int)GetRefinementLevel(cctkGH),pointcount,
               backup1,backup2,backup3,
               font_fixes,vel_limited_ptcount,rho_star_fix_applied,cons_avgs,atm_resets,
               failures,
               failures_inhoriz,pointcount_inhoriz,
               error_int_numer/error_int_denom,error_int_denom);
  }
  if( nan_found ) {
    if( GetRefinementLevel(cctkGH) > 6 ) {
      CCTK_ERROR("Found NAN during con2prim driver. See error messages above. ABORTING!");
    }
    else {
      CCTK_VWARN(CCTK_WARN_ALERT,"Found NAN during con2prim driver, but not at finest level. Proceeding with caution...");
    }
  }

//  // Very useful con2prim debugger. If the primitives (con2prim) solver fails, this will output all data needed to
//  //     debug where and why the solver failed. Strongly suggested for experimenting with new fixes.
//  if( ( (conserv_to_prims_debug==1) && (error_int_numer/error_int_denom > 0.05) ) ||
//      ( atm_resets != 0 ) ) {
//
//    ofstream myfile;
//    char filename[100];
//    srand(time(NULL));
//    sprintf(filename,"primitives_debug-%e.dat",error_int_numer/error_int_denom);
//    //Alternative, for debugging purposes as well:
//    //srand(time(NULL));
//    //sprintf(filename,"primitives_debug-%d.dat",rand());
//    myfile.open (filename, ios::out | ios::binary);
//    //myfile.open ("data.bin", ios::out | ios::binary);
//
//    // This checker value will be printed last, and will
//    // allow us to make sure we have read the dump file
//    // correctly when debugging it.
//    int checker=1063;
//
//    // Grid information
//    int fullsize=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
//    myfile.write((char*)cctk_lsh,                             3*sizeof(int));
//    myfile.write((char*)x,                           (fullsize)*sizeof(double));
//    myfile.write((char*)y,                           (fullsize)*sizeof(double));
//    myfile.write((char*)z,                           (fullsize)*sizeof(double));
//
//    // IllinoisGRMHD parameters
//    myfile.write((char*)&GAMMA_SPEED_LIMIT,                   1*sizeof(double));
//    myfile.write((char*)&rho_b_max,                           1*sizeof(double));
//    myfile.write((char*)&rho_b_atm,                           1*sizeof(double));
//    myfile.write((char*)&tau_atm,                             1*sizeof(double));
//    myfile.write((char*)&Psi6threshold,                       1*sizeof(double));
//    myfile.write((char*)&update_Tmunu,                        1*sizeof(bool));
//    myfile.write((char*)&neos,                                1*sizeof(int));
//    myfile.write((char*)&Gamma_th,                            1*sizeof(double));
//    myfile.write((char*)&K_ppoly_tab0,                        1*sizeof(double));
//    myfile.write((char*)Gamma_ppoly_tab_in,                neos*sizeof(double));
//    myfile.write((char*)rho_ppoly_tab_in,              (neos-1)*sizeof(double));
//    myfile.write((char*)&igm_Ye_atm,                          1*sizeof(double));
//    myfile.write((char*)&igm_T_atm,                           1*sizeof(double));
//    myfile.write((char*)&igm_eos_table_ceiling_safety_factor, 1*sizeof(double));
//    myfile.write((char*)&igm_eos_table_floor_safety_factor,   1*sizeof(double));
//    myfile.write((char*)&igm_eos_root_finding_precision,      1*sizeof(double));
//    myfile.write((char*)&igm_eos_root_finding_precision,      1*sizeof(double));
//    myfile.write((char*)&igm_evolve_temperature,              1*sizeof(bool));
//    myfile.write((char*)&igm_evolve_entropy,                  1*sizeof(bool));
//
//    // Failure checker gridfunction
//    myfile.write((char*)failure_checker,             (fullsize)*sizeof(double));
//
//    // Energy momentum tensor
//    myfile.write((char*)eTtt,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTtx,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTty,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTtz,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTxx,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTxy,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTxz,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTyy,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTyz,                        (fullsize)*sizeof(double));
//    myfile.write((char*)eTzz,                        (fullsize)*sizeof(double));
//
//    // Metric quantities
//    myfile.write((char*)alp,                         (fullsize)*sizeof(double));
//    myfile.write((char*)gxx,                         (fullsize)*sizeof(double));
//    myfile.write((char*)gxy,                         (fullsize)*sizeof(double));
//    myfile.write((char*)gxz,                         (fullsize)*sizeof(double));
//    myfile.write((char*)gyy,                         (fullsize)*sizeof(double));
//    myfile.write((char*)gyz,                         (fullsize)*sizeof(double));
//    myfile.write((char*)gzz,                         (fullsize)*sizeof(double));
//    myfile.write((char*)psi_bssn,                    (fullsize)*sizeof(double));
//    myfile.write((char*)phi_bssn,                    (fullsize)*sizeof(double));
//    myfile.write((char*)gtxx,                        (fullsize)*sizeof(double));
//    myfile.write((char*)gtxy,                        (fullsize)*sizeof(double));
//    myfile.write((char*)gtxz,                        (fullsize)*sizeof(double));
//    myfile.write((char*)gtyy,                        (fullsize)*sizeof(double));
//    myfile.write((char*)gtyz,                        (fullsize)*sizeof(double));
//    myfile.write((char*)gtzz,                        (fullsize)*sizeof(double));
//    myfile.write((char*)gtupxx,                      (fullsize)*sizeof(double));
//    myfile.write((char*)gtupxy,                      (fullsize)*sizeof(double));
//    myfile.write((char*)gtupxz,                      (fullsize)*sizeof(double));
//    myfile.write((char*)gtupyy,                      (fullsize)*sizeof(double));
//    myfile.write((char*)gtupyz,                      (fullsize)*sizeof(double));
//    myfile.write((char*)gtupzz,                      (fullsize)*sizeof(double));
//    myfile.write((char*)lapm1,                       (fullsize)*sizeof(double));
//    myfile.write((char*)betax,                       (fullsize)*sizeof(double));
//    myfile.write((char*)betay,                       (fullsize)*sizeof(double));
//    myfile.write((char*)betaz,                       (fullsize)*sizeof(double));
//
//    // Original conservative variables (we used the flux gridfunctions to store these)
//    myfile.write((char*)rho_star_flux,               (fullsize)*sizeof(double));
//    myfile.write((char*)tau_flux,                    (fullsize)*sizeof(double));
//    myfile.write((char*)st_x_flux,                   (fullsize)*sizeof(double));
//    myfile.write((char*)st_y_flux,                   (fullsize)*sizeof(double));
//    myfile.write((char*)st_z_flux,                   (fullsize)*sizeof(double));
//    myfile.write((char*)Ye_star_flux,                (fullsize)*sizeof(double));
//    myfile.write((char*)S_star_flux,                 (fullsize)*sizeof(double));
//
//    // Magnetic fields
//    myfile.write((char*)Bx,                          (fullsize)*sizeof(double));
//    myfile.write((char*)By,                          (fullsize)*sizeof(double));
//    myfile.write((char*)Bz,                          (fullsize)*sizeof(double));
//
//    // Primitive variables
//    myfile.write((char*)rho_b,                       (fullsize)*sizeof(double));
//    myfile.write((char*)P,                           (fullsize)*sizeof(double));
//    myfile.write((char*)vx,                          (fullsize)*sizeof(double));
//    myfile.write((char*)vy,                          (fullsize)*sizeof(double));
//    myfile.write((char*)vz,                          (fullsize)*sizeof(double));
//    myfile.write((char*)igm_Ye,                      (fullsize)*sizeof(double));
//    myfile.write((char*)igm_temperature,             (fullsize)*sizeof(double));
//    myfile.write((char*)igm_eps,                     (fullsize)*sizeof(double));
//    myfile.write((char*)igm_entropy,                 (fullsize)*sizeof(double));
//
//    // Checker value
//    myfile.write((char*)&checker,                             1*sizeof(int));
//
//    // All done! Close the file.
//    myfile.close();
//    CCTK_VInfo(CCTK_THORNSTRING,"Finished writing %s",filename);
//  }

#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER
  return 0; // int main() requires an integer be returned
#endif

}

#include "harm_u2p_util.h"
#include "con2prim_wrapper.h"
