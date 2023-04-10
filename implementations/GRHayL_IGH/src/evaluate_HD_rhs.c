/*********************************************
 * Evaluate RHS of GRHD & induction equations
 * (vector potential prescription), using the
 * generalized Lorenz gauge condition for the
 * EM gauge.
 *
 * Based originally on the Illinois GRHD code,
 * written by Matt Duez, Yuk Tung Liu, and Branson
 * Stephens (original version), and then developed
 * primarily by Zachariah Etienne, Yuk Tung Liu,
 * and Vasileios Paschalidis.
 *
 * Rewritten for public release in 2013
 *      by Zachariah B. Etienne
 *
 * References:
 * Original unigrid GRHD evolution prescription:
 *    http://arxiv.org/abs/astro-ph/0503420
 * Vector potential formulation in full GR:
 *    http://arxiv.org/abs/1007.2848
 * Improved EM gauge conditions for AMR grids:
 *    http://arxiv.org/abs/1110.4633
 * Generalized Lorenz gauge prescription:
 *    http://arxiv.org/abs/1207.3354
 *
 * Note that the Generalized Lorenz gauge strength
 *  parameter has units of 1/M, just like the \eta
 *  parameter in the gamma-driving shift condition,
 *  so setting it too large will result in violation
 *  of the CFL condition.
 *
 * This version of PPM implements the standard
 * Colella & Woodward PPM, though modified as in GRHydro
 * to have 3 ghostzones instead of 4.
 *********************************************/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "IGH.h"

void GRHayL_IGH_evaluate_HD_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGH_evaluate_HD_rhs;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_Equals(verbose, "essential+iteration output")) {
    const int levelnumber = GetRefinementLevel(cctkGH);
    CCTK_VINFO("***** Iter. # %d, Lev: %d, Integrating to time: %e *****",cctk_iteration,levelnumber,cctk_delta_time/cctk_levfac[0]+cctk_time);
  }

  if( sizeof(CCTK_REAL) < 8 ) CCTK_VERROR("Error: GRHayL_IGH assumes that CCTK_REAL is a double precision number. Setting otherwise will likely cause havoc with the conserv_to_prims solver.");

  if(cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3) { CCTK_VERROR("ERROR. Need at least 3 ghostzones for GRHayL_IGH evolutions."); }

  CCTK_REAL dX[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };


  // in_prims,out_prims_r, and out_prims_l are arrays of pointers to the actual gridfunctions.
  // Can't use restrict because pointers for e.g. vxr are in in_prims and out_prims_r.
  // TODO: consider a better way of storing/passing these pointers since we no longer have
  //       to track the ghostzone information
  const double *in_prims[MAXNUMVARS];
  double *out_prims_r[MAXNUMVARS];
  double *out_prims_l[MAXNUMVARS];

  /* SET POINTERS TO GRHD GRIDFUNCTIONS */
  // The order here MATTERS, and must be consistent with the global variable declarations in
  //   evaluate_HD_rhs_headers.h (look for RHOB=0, etc.)
  //   For example, in_prims[0] _must_ be rho_b.
  int ww=0;
  in_prims[ww]=rho_b;      out_prims_r[ww]=rhor;        out_prims_l[ww]=rhol;        ww++;
  in_prims[ww]=pressure;   out_prims_r[ww]=pressr;      out_prims_l[ww]=pressl;      ww++;
  in_prims[ww]=vx;         out_prims_r[ww]=vxr;         out_prims_l[ww]=vxl;         ww++;
  in_prims[ww]=vy;         out_prims_r[ww]=vyr;         out_prims_l[ww]=vyl;         ww++;
  in_prims[ww]=vz;         out_prims_r[ww]=vzr;         out_prims_l[ww]=vzl;         ww++;

  double *cmin[3] = {cmin_x, cmin_y, cmin_z};
  double *cmax[3] = {cmax_x, cmax_y, cmax_z};

  // Convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  GRHayL_IGH_convert_ADM_to_BSSN(cctkGH,
                             gxx, gxy, gxz, gyy, gyz, gzz,
                             phi_bssn, psi_bssn,
                             gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,
                             gtupxx, gtupxy, gtupxz, gtupyy, gtupyz, gtupzz);

//  /* SET POINTERS TO METRIC GRIDFUNCTIONS */
  const double *metric[10]; // "metric" here is array of pointers to the actual gridfunctions.
  metric[LAPSE] = alp;
  metric[BETAX] = betax;
  metric[BETAY] = betay;
  metric[BETAZ] = betaz;
  metric[GXX]   = gxx;
  metric[GXY]   = gxy;
  metric[GXZ]   = gxz;
  metric[GYY]   = gyy;
  metric[GYZ]   = gyz;
  metric[GZZ]   = gzz;

  const double *curv[6];
  curv[KXX] = kxx;
  curv[KXY] = kxy;
  curv[KXZ] = kxz;
  curv[KYY] = kyy;
  curv[KYZ] = kyz;
  curv[KZZ] = kzz;

  // 1) First initialize RHS variables to zero
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        tau_rhs[index]      = 0.0;
        rho_star_rhs[index] = 0.0;
        Stildex_rhs[index]  = 0.0;
        Stildey_rhs[index]  = 0.0;
        Stildez_rhs[index]  = 0.0;
  }

  // Here, we:
  // 1) Compute tau_rhs extrinsic curvature terms, and
  // 2) Compute TUPmunu.
  // This function is housed in the file: "compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu.C"
  GRHayL_IGH_calculate_tau_source_rhs(cctkGH, grhayl_eos, metric, curv, in_prims, tau_rhs);

  int flux_dir;
  const int var_indices[3] = {VX, VY, VZ};
  flux_dir=0;
  /*
   *  Computation of \partial_x on RHS of \partial_t {rho_star,tau,mhd_st_{x,y,z}},
   *  via PPM reconstruction onto (i-1/2,j,k), so that
   *  \partial_x F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   */
  GRHayL_IGH_reconstruction_loop(cctkGH, flux_dir, 3, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);

  GRHayL_IGH_compute_characteristic_speeds(cctkGH, flux_dir, grhayl_eos, metric,
                                out_prims_r, out_prims_l, cmin[flux_dir], cmax[flux_dir]);

  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  GRHayL_IGH_calculate_HD_dirn_rhs(cctkGH, flux_dir, dX, grhayl_eos, metric, in_prims,
                         out_prims_r, out_prims_l, cmin[flux_dir], cmax[flux_dir],
                         rho_star_flux, tau_flux, Stildex_flux, Stildey_flux, Stildez_flux,
                         rho_star_rhs, tau_rhs, Stildex_rhs, Stildey_rhs, Stildez_rhs);

  /*
   *  Computation of \partial_y on RHS of \partial_t {rho_star,tau,mhd_st_{x,y,z}},
   *  via PPM reconstruction onto (i,j-1/2,k), so that
   *  \partial_y F = [ F(i,j+1/2,k) - F(i,j-1/2,k) ] / dy
   */
  flux_dir=1;
  GRHayL_IGH_reconstruction_loop(cctkGH, flux_dir, 3, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);

  GRHayL_IGH_compute_characteristic_speeds(cctkGH, flux_dir, grhayl_eos, metric,
                                out_prims_r, out_prims_l, cmin[flux_dir], cmax[flux_dir]);

  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  GRHayL_IGH_calculate_HD_dirn_rhs(cctkGH, flux_dir, dX, grhayl_eos, metric, in_prims,
                         out_prims_r, out_prims_l, cmin[flux_dir], cmax[flux_dir],
                         rho_star_flux, tau_flux, Stildex_flux, Stildey_flux, Stildez_flux,
                         rho_star_rhs, tau_rhs, Stildex_rhs, Stildey_rhs, Stildez_rhs);

  flux_dir=2;

  /*
   *  Computation of \partial_z on RHS of \partial_t {rho_star,tau,mhd_st_{x,y,z}},
   *  via PPM reconstruction onto (i,j,k-1/2), so that
   *  \partial_z F = [ F(i,j,k+1/2) - F(i,j,k-1/2) ] / dz
   */
  GRHayL_IGH_reconstruction_loop(cctkGH, flux_dir, 3, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);

  GRHayL_IGH_compute_characteristic_speeds(cctkGH, flux_dir, grhayl_eos, metric,
                                out_prims_r, out_prims_l, cmin[flux_dir], cmax[flux_dir]);

  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  GRHayL_IGH_calculate_HD_dirn_rhs(cctkGH, flux_dir, dX, grhayl_eos, metric, in_prims,
                         out_prims_r, out_prims_l, cmin[flux_dir], cmax[flux_dir],
                         rho_star_flux, tau_flux, Stildex_flux, Stildey_flux, Stildez_flux,
                         rho_star_rhs, tau_rhs, Stildex_rhs, Stildey_rhs, Stildez_rhs);
}
