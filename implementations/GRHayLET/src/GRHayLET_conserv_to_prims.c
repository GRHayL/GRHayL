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


#include "cctk.h"
//#include <iostream>
//#include <iomanip>
#include <fstream>
//#include <sys/time.h>
//#include <cmath>
//#include <ctime>
//#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "con2prim.h"
#include "IllinoisGRMHD_headers.h"
#include "harm_primitives_headers.h"
#include "inlined_functions.C"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"

void IllinoisGRMHD_conserv_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  double poison = 1e200;
  int main = Noble2D;
  int backup_routines[3] = {None,None,None};

  GRHayL_parameters params;
  initialize_parameters_ETK(main, backup_routines, false, false, true, Psi6threshold, update_Tmunu, Cupp_Fix, (void*)&params);

  eos_parameters eos;
  int eos_type = 0; // Hybrid=0, Tabulated=1;
  initialize_general_eos_ETK(eos_type, GAMMA_SPEED_LIMIT,
             rho_b_atm, rho_b_atm, rho_b_max,
             (void*)&eos);

  initialize_hybrid_eos_ETK(
             neos, rho_tab,
             gamma_tab, k_tab[0],
             gamma_th, (void*)&eos);

  // These BSSN-based variables are not evolved, and so are not defined anywhere that the grid has moved.
  // Here we convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                phi_bssn,psi_bssn,lapm1);

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

  //Start the timer, so we can benchmark the primitives solver during evolution.
  //  Slower solver -> harder to find roots -> things may be going crazy!
  //FIXME: Replace this timing benchmark with something more meaningful, like the avg # of Newton-Raphson iterations per gridpoint!
  /*
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);
  */

//  int pressure_cap_hit=0;

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

#pragma omp parallel for reduction(+:failures,font_fixes,vel_limited_ptcount,atm_resets,rho_star_fix_applied,pointcount,failures_inhoriz,pointcount_inhoriz,backup0,backup1,backup2,nan_found,error_int_numer,error_int_denom,n_iter) schedule(static)
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    con2prim_diagnostics diagnostics_private;
    initialize_diagnostics_ETK((void*)&diagnostics_private);
//    diagnostics_private.c2p_fail_flag = con2prim_failed_flag[index]; This is also from Leo's IGM

    // Read in BSSN metric quantities from gridfunctions and
    // set auxiliary and ADM metric quantities
    metric_quantities metric;
    initialize_metric_ETK(alp[index],
                      gxx[index], gxy[index], gxz[index],
                      gyy[index], gyz[index], gzz[index],
                      betax[index], betay[index], betaz[index],
                      (void*)&metric);
       
    // Read in primitive variables from gridfunctions
    primitive_quantities prims, prims_orig;
    initialize_primitives_ETK(rho_b[index],
                          P[index], eps[index],
                          vx[index], vy[index], vz[index],
                          Bx[index], By[index], Bz[index],
                          poison, poison, poison,
                          //entropy[index], Y_e[index], temp[index],
                          (void*)&prims);

    // Read in conservative variables from gridfunctions
    conservative_quantities cons, cons_orig;
    initialize_conservatives_ETK(
                             rho_star[index], tau[index],
                             mhd_st_x[index], mhd_st_y[index], mhd_st_z[index],
                             poison, poison,
                             //Y_e[index], entropy[index],
                             (void*)&cons);

    prims_orig = prims;
    cons_orig = cons;

    stress_energy Tmunu;
    Con2Prim_Kernel((void*)&params, (void*)&eos, (void*)&metric, (void*)&cons, (void*)&prims, (void*)&diagnostics_private, (void*)&Tmunu);

    failure_checker[index] = diagnostics_private.failure_checker;

    double dummy;
    return_primitives_ETK((void*)&prims, &rho_b[index], &P[index], &eps[index],
                  &vx[index], &vy[index], &vz[index], &Bx[index], &By[index], &Bz[index],
                  &dummy, &dummy, &dummy);
                  //&entropy[index], &Y_e[index], &temp[index]);

    return_conservatives_ETK((void*)&cons, &rho_star[index], &tau[index],
                        &mhd_st_x[index], &mhd_st_y[index], &mhd_st_z[index],
                        &dummy, &dummy);
                        //&Y_e[index], &entropy[index]);

    if(update_Tmunu) {
      eTtt[index] = Tmunu.Ttt;
      eTtx[index] = Tmunu.Ttx;
      eTty[index] = Tmunu.Tty;
      eTtz[index] = Tmunu.Ttz;
      eTxx[index] = Tmunu.Txx;
      eTxy[index] = Tmunu.Txy;
      eTxz[index] = Tmunu.Txz;
      eTyy[index] = Tmunu.Tyy;
      eTyz[index] = Tmunu.Tyz;
      eTzz[index] = Tmunu.Tzz;
    }

    failures += diagnostics_private.failures;
    font_fixes += diagnostics_private.font_fixes;
    vel_limited_ptcount += diagnostics_private.vel_limited_ptcount;
    atm_resets += diagnostics_private.atm_resets;
    rho_star_fix_applied += diagnostics_private.rho_star_fix_applied;
    pointcount += diagnostics_private.pointcount;
    failures_inhoriz += diagnostics_private.failures_inhoriz;
    pointcount_inhoriz += diagnostics_private.pointcount_inhoriz;
    backup0 += diagnostics_private.backup[0];
    backup1 += diagnostics_private.backup[1];
    backup2 += diagnostics_private.backup[2];
    nan_found += diagnostics_private.nan_found;
    error_int_numer += diagnostics_private.error_int_numer;
    error_int_denom += diagnostics_private.error_int_denom;
    n_iter += diagnostics_private.n_iter;
}

  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VINFO("C2P: Lev: %d NumPts= %d | Fixes: Font= %d VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e | %.2f iters/gridpt",
               (int)GetRefinementLevel(cctkGH), pointcount,
               font_fixes, vel_limited_ptcount, rho_star_fix_applied,
               failures, failures_inhoriz, pointcount_inhoriz,
               error_int_numer/error_int_denom, error_int_denom,
               (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ));
  }
//
//  if(pressure_cap_hit!=0) {
//    //CCTK_VInfo(CCTK_THORNSTRING,"PRESSURE CAP HIT %d TIMES!  Outputting debug file!",pressure_cap_hit);
//  }
//
//  // Very useful con2prim debugger. If the primitives (con2prim) solver fails, this will output all data needed to
//  //     debug where and why the solver failed. Strongly suggested for experimenting with new fixes.
//  if(conserv_to_prims_debug==1) {
//
//    ofstream myfile;
//    char filename[100];
//    srand(time(NULL));
//    sprintf(filename,"primitives_debug-%e.dat",rho_star[CCTK_GFINDEX3D(cctkGH,3,15,6)]);
//    //Alternative, for debugging purposes as well: sprintf(filename,"primitives_debug-%d.dat",rand());
//    myfile.open (filename, ios::out | ios::binary);
//    //myfile.open ("data.bin", ios::out | ios::binary);
//    myfile.write((char*)cctk_lsh, 3*sizeof(int));
//
//    myfile.write((char*)&rho_b_atm, 1*sizeof(CCTK_REAL));
//    myfile.write((char*)&tau_atm, 1*sizeof(CCTK_REAL));
//
//    myfile.write((char*)&Psi6threshold, 1*sizeof(CCTK_REAL));
//
//    CCTK_REAL gamma_th=2.0;
//    myfile.write((char*)&gamma_th, 1*sizeof(CCTK_REAL));
//    myfile.write((char*)&neos,     1*sizeof(int));
//
//    myfile.write((char*)gamma_tab, (neos+1)*sizeof(CCTK_REAL));
//    myfile.write((char*)k_tab,     (neos+1)*sizeof(CCTK_REAL));
//
//    myfile.write((char*)eps_tab,   neos*sizeof(CCTK_REAL));
//    myfile.write((char*)rho_tab,   neos*sizeof(CCTK_REAL));
//    myfile.write((char*)P_tab,     neos*sizeof(CCTK_REAL));
//
//    int fullsize=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
//    myfile.write((char*)x,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)y,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)z,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)phi_bssn, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtxx, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtxy, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtxz, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtyy, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtyz, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtzz, (fullsize)*sizeof(CCTK_REAL));
//
//    myfile.write((char*)gtupxx, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtupxy, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtupxz, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtupyy, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtupyz, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)gtupzz, (fullsize)*sizeof(CCTK_REAL));
//
//    myfile.write((char*)betax, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)betay, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)betaz, (fullsize)*sizeof(CCTK_REAL));
//
//    myfile.write((char*)lapm1, (fullsize)*sizeof(CCTK_REAL));
//
//    // HERE WE USE _flux variables as temp storage for original values of conservative variables.. This is used for debugging purposes only.
//    myfile.write((char*)tau_flux,      (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)st_x_flux, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)st_y_flux, (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)st_z_flux, (fullsize)*sizeof(CCTK_REAL));
//
//    myfile.write((char*)rho_star_flux, (fullsize)*sizeof(CCTK_REAL));
//
//    myfile.write((char*)Bx,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)By,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)Bz,   (fullsize)*sizeof(CCTK_REAL));
//
//    myfile.write((char*)vx,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)vy,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)vz,   (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)P,    (fullsize)*sizeof(CCTK_REAL));
//    myfile.write((char*)rho_b,(fullsize)*sizeof(CCTK_REAL));
//
//    int checker=1063; myfile.write((char*)&checker,sizeof(int));
//
//    myfile.close();
//    CCTK_VInfo(CCTK_THORNSTRING,"Finished writing...");
//  }
}

#include "harm_primitives_lowlevel.C"
