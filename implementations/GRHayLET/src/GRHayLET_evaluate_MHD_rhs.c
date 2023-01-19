/*********************************************
 * Evaluate RHS of GRMHD & induction equations
 * (vector potential prescription), using the
 * generalized Lorenz gauge condition for the
 * EM gauge.
 *
 * Based originally on the Illinois GRMHD code,
 * written by Matt Duez, Yuk Tung Liu, and Branson
 * Stephens (original version), and then developed
 * primarily by Zachariah Etienne, Yuk Tung Liu,
 * and Vasileios Paschalidis.
 *
 * Rewritten for public release in 2013
 *      by Zachariah B. Etienne
 *
 * References:
 * Original unigrid GRMHD evolution prescription:
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
#include "GRHayLET.h"

void GRHayLET_evaluate_MHD_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_Equals(verbose, "essential+iteration output")) {
    const int levelnumber = GetRefinementLevel(cctkGH);
    CCTK_VINFO("***** Iter. # %d, Lev: %d, Integrating to time: %e *****",cctk_iteration,levelnumber,cctk_delta_time/cctk_levfac[0]+cctk_time);
  }

  if( sizeof(CCTK_REAL) < 8 ) CCTK_VERROR("Error: IllinoisGRMHD assumes that CCTK_REAL is a double precision number. Setting otherwise will likely cause havoc with the conserv_to_prims solver.");

  if(cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3) { CCTK_VERROR("ERROR. Need at least 3 ghostzones for IllinoisGRMHD evolutions."); }

  CCTK_REAL dX[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };


  // in_prims,out_prims_r, and out_prims_l are arrays of pointers to the actual gridfunctions.
  // Can't use restrict because pointers for e.g. vxr are in in_prims and out_prims_r.
  // TODO: consider a better way of storing/passing these pointers since we no longer have
  //       to track the ghostzone information
  const double *in_prims[MAXNUMVARS];
  double *out_prims_r[MAXNUMVARS];
  double *out_prims_l[MAXNUMVARS];
  int num_prims_to_reconstruct;

  /* SET POINTERS TO GRMHD GRIDFUNCTIONS */
  // The order here MATTERS, and must be consistent with the global variable declarations in
  //   evaluate_MHD_rhs_headers.h (look for RHOB=0, etc.)
  //   For example, in_prims[0] _must_ be rho_b.
  int ww=0;
  in_prims[ww]=rho;        out_prims_r[ww]=rhor;        out_prims_l[ww]=rhol;        ww++;
  in_prims[ww]=press;      out_prims_r[ww]=pressr;      out_prims_l[ww]=pressl;      ww++;
  in_prims[ww]=&vel[0];    out_prims_r[ww]=vxr;         out_prims_l[ww]=vxl;         ww++;
  in_prims[ww]=&vel[1];    out_prims_r[ww]=vyr;         out_prims_l[ww]=vyl;         ww++;
  in_prims[ww]=&vel[2];    out_prims_r[ww]=vzr;         out_prims_l[ww]=vzl;         ww++;
  in_prims[ww]=&Bvec[0];   out_prims_r[ww]=Bxr;         out_prims_l[ww]=Bxl;         ww++;
  in_prims[ww]=&Bvec[1];   out_prims_r[ww]=Byr;         out_prims_l[ww]=Byl;         ww++;
  in_prims[ww]=&Bvec[2];   out_prims_r[ww]=Bzr;         out_prims_l[ww]=Bzl;         ww++;
  in_prims[ww]=Bx_stagger; out_prims_r[ww]=Bx_staggerr; out_prims_l[ww]=Bx_staggerl; ww++;
  in_prims[ww]=By_stagger; out_prims_r[ww]=By_staggerr; out_prims_l[ww]=By_staggerl; ww++;
  in_prims[ww]=Bz_stagger; out_prims_r[ww]=Bz_staggerr; out_prims_l[ww]=Bz_staggerl; ww++;
  in_prims[ww]=vxr;        out_prims_r[ww]=vxrr;        out_prims_l[ww]=vxrl;        ww++;
  in_prims[ww]=vyr;        out_prims_r[ww]=vyrr;        out_prims_l[ww]=vyrl;        ww++;
  in_prims[ww]=vzr;        out_prims_r[ww]=vzrr;        out_prims_l[ww]=vzrl;        ww++;
  in_prims[ww]=vxl;        out_prims_r[ww]=vxlr;        out_prims_l[ww]=vxll;        ww++;
  in_prims[ww]=vyl;        out_prims_r[ww]=vylr;        out_prims_l[ww]=vyll;        ww++;
  in_prims[ww]=vzl;        out_prims_r[ww]=vzlr;        out_prims_l[ww]=vzll;        ww++;

//TODO: this function needs to be pulled into GRHayLET first
  // Convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
//  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
//                                                                gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
//                                                                gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
//                                                                phi_bssn,psi_bssn,lapm1);

//  /* SET POINTERS TO METRIC GRIDFUNCTIONS */
//  CCTK_REAL *metric[NUMVARS_FOR_METRIC_FACEVALS]; // "metric" here is array of pointers to the actual gridfunctions.
//  ww=0;
//  metric[ww]=phi_bssn;ww++;
//  metric[ww]=psi_bssn;ww++;
//  metric[ww]=gtxx;    ww++;
//  metric[ww]=gtxy;    ww++;
//  metric[ww]=gtxz;    ww++;
//  metric[ww]=gtyy;    ww++;
//  metric[ww]=gtyz;    ww++;
//  metric[ww]=gtzz;    ww++;
//  metric[ww]=lapm1;   ww++;
//  metric[ww]=betax;   ww++;
//  metric[ww]=betay;   ww++;
//  metric[ww]=betaz;   ww++;
//  metric[ww]=gtupxx;  ww++;
//  metric[ww]=gtupyy;  ww++;
//  metric[ww]=gtupzz;  ww++;
//
//  /* SET POINTERS TO STRESS-ENERGY TENSOR GRIDFUNCTIONS */
//  CCTK_REAL *TUPmunu[10];// "TUPmunu" here is array of pointers to the actual gridfunctions.
//  ww=0;
//  TUPmunu[ww]=TUPtt; ww++;
//  TUPmunu[ww]=TUPtx; ww++;
//  TUPmunu[ww]=TUPty; ww++;
//  TUPmunu[ww]=TUPtz; ww++;
//  TUPmunu[ww]=TUPxx; ww++;
//  TUPmunu[ww]=TUPxy; ww++;
//  TUPmunu[ww]=TUPxz; ww++;
//  TUPmunu[ww]=TUPyy; ww++;
//  TUPmunu[ww]=TUPyz; ww++;
//  TUPmunu[ww]=TUPzz; ww++;

  // 1) First initialize {rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs} to zero
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        Ax_rhs[index]=0.0;
        Ay_rhs[index]=0.0;
        Az_rhs[index]=0.0;
        phitilde_rhs[index]=0.0;

//        tau_rhs[index]=0.0;
//        rho_star_rhs[index]=0.0;
//        st_x_rhs[index]=0.0;
//        st_y_rhs[index]=0.0;
//        st_z_rhs[index]=0.0;
  }

//TODO: pull compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu into GRHayLET
  // Here, we:
  // 1) Compute tau_rhs extrinsic curvature terms, and
  // 2) Compute TUPmunu.
  // This function is housed in the file: "compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu.C"
//  compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu(cctkGH,cctk_lsh,cctk_nghostzones,dX,metric,in_prims,TUPmunu,grhayl_eos,
//                                                        gtupxy,gtupxz,gtupyz,
//                                                        kxx,kxy,kxz,kyy,kyz,kzz,
//                                                        tau_rhs);

  int flux_dir;
  flux_dir=1;
//  /* There are two stories going on here:
//   * 1) Computation of \partial_x on RHS of \partial_t {rho_star,tau,mhd_st_{x,y,z}},
//   *    via PPM reconstruction onto (i-1/2,j,k), so that
//   *    \partial_x F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
//   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
//   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
//   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} v^j B^k,
//   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
//   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
//   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
//   *     staggered points. For example:
//   * 2Aa) vx and vy are at (i,j,k), and we reconstruct them to (i-1/2,j,k) below. After
//   *      this, we'll reconstruct again in the y-dir'n to get {vx,vy} at (i-1/2,j-1/2,k)
//   * 2Ab) By_stagger is at (i,j+1/2,k), and we reconstruct below to (i-1/2,j+1/2,k). */
/**************************************new code*******************************************/
  {const int num_vars = 6;
  num_prims_to_reconstruct=num_vars+2;
  const int var_indices[6] = {VX, VY, VZ, BY_CENTER, BZ_CENTER, BY_STAGGER};
  reconstruction_loop(cctkGH, flux_dir, num_vars, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);
  }
/*****************************************************************************************/

  //Right and left face values of BI_CENTER are used in mhdflux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^x face values to be consistent with BX_STAGGER.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k), indexim1=CCTK_GFINDEX3D(cctkGH,i-1+(i==0),j,k); /* indexim1=0 when i=0 */
        out_prims_r[BX_CENTER][index]=out_prims_l[BX_CENTER][index]=in_prims[BX_STAGGER][indexim1]; }

//TODO: add in Terrence's code
  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
//  add_fluxes_and_source_terms_to_hydro_rhss(flux_dir,cctkGH,cctk_lsh,cctk_nghostzones,dX,   metric,in_prims,TUPmunu,
//                                            num_prims_to_reconstruct,out_prims_r,out_prims_l,eos,
//                                            cmax_x,cmin_x,
//                                            rho_star_flux,tau_flux,st_x_flux,st_y_flux,st_z_flux,
//                                            rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs);

  // Note that we have already reconstructed vx and vy along the x-direction,
  //   at (i-1/2,j,k). That result is stored in v{x,y}{r,l}.  Bx_stagger data
  //   are defined at (i+1/2,j,k).
  // Next goal: reconstruct Bx, vx and vy at (i+1/2,j+1/2,k).
  flux_dir=2;

  /* There are two stories going on here:
   * 1) Computation of \partial_y on RHS of \partial_t {rho_star,tau,mhd_st_{x,y,z}},
   *    via PPM reconstruction onto (i,j-1/2,k), so that
   *    \partial_y F = [ F(i,j+1/2,k) - F(i,j-1/2,k) ] / dy
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} v^j B^k,
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) VXR = [right-face of vx reconstructed along x-direction above] is at (i-1/2,j,k),
   *      and we reconstruct it to (i-1/2,j-1/2,k) below. Similarly for {VXL,VYR,VYL}
   * 2Ab) Bx_stagger is at (i+1/2,j,k), and we reconstruct to (i+1/2,j-1/2,k) below
   * 2Ac) By_stagger is at (i-1/2,j+1/2,k) already for Az_rhs, from the previous step.
   * 2B) Ax_rhs is defined at (i,j+1/2,k+1/2), and it depends on {By,Bz,vy,vz}.
   *     Again the trick is to reconstruct these onto these staggered points.
   * 2Ba) Bz_stagger is at (i,j,k+1/2), and we reconstruct to (i,j-1/2,k+1/2) below */
  //// NOTE! The order of variable reconstruction is important here,
  ////   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
/**************************************new code*******************************************/
  {const int num_vars = 4;
  const int var_indices[4] = {VXR, VYR, VXL, VYL};
  reconstruction_loop_no_rho_P(cctkGH, flux_dir, num_vars, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);
  }
  {const int num_vars = 7;
  num_prims_to_reconstruct = num_vars+2;
  const int var_indices[7] = {VX, VY, VZ, BX_CENTER, BZ_CENTER, BX_STAGGER, BZ_STAGGER};
  reconstruction_loop(cctkGH, flux_dir, num_vars, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);
  }
/*****************************************************************************************/

  //Right and left face values of BI_CENTER are used in mhdflux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^y face values to be consistent with BY_STAGGER.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k), indexjm1=CCTK_GFINDEX3D(cctkGH,i,j-1+(j==0),k); /* indexjm1=0 when j=0 */
        out_prims_r[BY_CENTER][index]=out_prims_l[BY_CENTER][index]=in_prims[BY_STAGGER][indexjm1]; }

//TODO: add in Terrence's code
  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
//  add_fluxes_and_source_terms_to_hydro_rhss(flux_dir,cctkGH,cctk_lsh,cctk_nghostzones,dX,   metric,in_prims,TUPmunu,
//                                            num_prims_to_reconstruct,out_prims_r,out_prims_l,eos,
//                                            cmax_y,cmin_y,
//                                            rho_star_flux,tau_flux,st_x_flux,st_y_flux,st_z_flux,
//                                            rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs);

  /*****************************************
   * COMPUTING RHS OF A_z, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_z - [gauge terms] = \psi^{6} (v^x B^y - v^y B^x).
   * A_z is defined at (i+1/2,j+1/2,k).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j-1/2,k)| {vxrr,vxrl,vxlr,vxll,vyrr,vyrl,vylr,vyll}
   * (i+1/2,j-1/2,k)| {Bx_stagger_r,Bx_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j+1/2,k)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Interpolates to i+1/2
  flux_dir=3;
  // cmin/max could be done internally using the same indices as v and B if all c var pointers were collected into
  // single array
  A_no_gauge_rhs(cctkGH, flux_dir, out_prims_r, out_prims_l, phi_bssn, cmin_x, cmax_x, cmin_y, cmax_y, Az_rhs);

  /* There are two stories going on here:
   * 1) Single reconstruction to (i,j,k-1/2) for {rho,P,vx,vy,vz,Bx,By,Bz} to compute
   *    z-dir'n advection terms in \partial_t {rho_star,tau,mhd_st_{x,y,z}} at (i,j,k)
   * 2) Multiple reconstructions for *staggered* gridfunctions A_i:
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} v^j B^k,
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Ax_rhs is defined at (i,j+1/2,k+1/2), depends on v{y,z} and B{y,z}
   * 2Aa) v{y,z}{r,l} are at (i,j-1/2,k), so we reconstruct here to (i,j-1/2,k-1/2)
   * 2Ab) Bz_stagger{r,l} are at (i,j-1/2,k+1/2) already.
   * 2Ac) By_stagger is at (i,j+1/2,k), and below we reconstruct its value at (i,j+1/2,k-1/2)
   * 2B) Ay_rhs is defined at (i+1/2,j,k+1/2), depends on v{z,x} and B{z,x}.
   * 2Ba) v{x,z} are reconstructed to (i,j,k-1/2). Later we'll reconstruct again to (i-1/2,j,k-1/2).
   * 2Bb) Bz_stagger is at (i,j,k+1/2). Later we will reconstruct to (i-1/2,j,k+1/2).
   * 2Bc) Bx_stagger is at (i+1/2,j,k), and below we reconstruct its value at (i+1/2,j,k-1/2)
   */
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vyr,vyl,vzr,vzl}!
/**************************************new code*******************************************/
  {const int num_vars = 4;
  const int var_indices[4] = {VYR, VZR, VYL, VZL};
  reconstruction_loop_no_rho_P(cctkGH, flux_dir, num_vars, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);
  }
  {const int num_vars = 7;
  num_prims_to_reconstruct = num_vars+2;
  const int var_indices[7] = {VX, VY, VZ, BX_CENTER, BY_CENTER, BX_STAGGER, BY_STAGGER};
  reconstruction_loop(cctkGH, flux_dir, num_vars, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);
  }
/*****************************************************************************************/

  //Right and left face values of BI_CENTER are used in mhdflux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^z face values to be consistent with BZ_STAGGER.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k), indexkm1=CCTK_GFINDEX3D(cctkGH,i,j,k-1+(k==0)); /* indexkm1=0 when k=0 */
        out_prims_r[BZ_CENTER][index]=out_prims_l[BZ_CENTER][index]=in_prims[BZ_STAGGER][indexkm1]; }

//TODO: add in Terrence's code
  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
//  add_fluxes_and_source_terms_to_hydro_rhss(flux_dir,cctkGH,cctk_lsh,cctk_nghostzones,dX,   metric,in_prims,TUPmunu,
//                                            num_prims_to_reconstruct,out_prims_r,out_prims_l,eos,
//                                            cmax_z,cmin_z,
//                                            rho_star_flux,tau_flux,st_x_flux,st_y_flux,st_z_flux,
//                                            rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs);

  /*****************************************
   * COMPUTING RHS OF A_x, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_x - [gauge terms] = \psi^{6} (v^y B^z - v^z B^y).
   * A_x is defined at (i,j+1/2,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i,j-1/2,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i,j+1/2,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j-1/2,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  flux_dir=1;
  // cmin/max could be done internally using the same indices as v and B if all c var pointers were collected into
  // single array
  A_no_gauge_rhs(cctkGH, flux_dir, out_prims_r, out_prims_l, phi_bssn, cmin_y, cmax_y, cmin_z, cmax_z, Ax_rhs);

  // We reprise flux_dir=1 reconstruction to finish up computations of Ai_rhs's!
/**************************************new code*******************************************/
  {const int num_vars = 5;
  const int var_indices[5] = {VXR, VZR, VXL, VZL, BZ_STAGGER};
  reconstruction_loop_no_rho_P(cctkGH, flux_dir, num_vars, var_indices, grhayl_eos, in_prims, out_prims_r, out_prims_l);
  }
/*****************************************************************************************/


  /*****************************************
   * COMPUTING RHS OF A_y, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_y - [gauge terms] = \psi^{6} (v^z B^x - v^x B^z).
   * A_y is defined at (i+1/2,j,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j,k-1/2)| {vxrr,vxrl,vxlr,vxll,vzrr,vzrl,vzlr,vzll}
   * (i+1/2,j,k-1/2)| {Bx_stagger_r,Bx_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  flux_dir=2;
  // cmin/max could be done internally using the same indices as v and B if all c var pointers were collected into
  // single array
  A_no_gauge_rhs(cctkGH, flux_dir, out_prims_r, out_prims_l, phi_bssn, cmin_z, cmax_z, cmin_x, cmax_x, Ay_rhs);

  // Next compute phitilde_rhs, and add gauge terms to A_i_rhs terms!
  //   Note that in the following function, we don't bother with reconstruction, instead interpolating.
  // We need A^i, but only have A_i. So we add gtupij to the list of input variables.
  // We are FINISHED with v{x,y,z}{r,l} and P{r,l} so we use these 8 gridfunctions' worth of space as temp storage.
  phitilde_and_A_gauge_rhs(cctkGH, dX, gtupxx, gtupxy, gtupxz, gtupyy, gtupyz, gtupzz,
              psi_bssn, lapm1, betax, betay, betaz, Ax, Ay, Az, phitilde,
              damp_lorenz, vxr, vyr, vzr, vxl, vyl, vzl, pressr, pressl,
              phitilde_rhs, Ax_rhs, Ay_rhs, Az_rhs);

  return;
  /*
  // FUN DEBUGGING TOOL (trust me!):
  #pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
  int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
  //st_x_rhs[index]=0.0;
  //st_y_rhs[index]=0.0;
  //st_z_rhs[index]=0.0;
  //rho_star_rhs[index]=0.0;
  //tau_rhs[index]=0.0;

  phitilde_rhs[index] = 0.0;
  Ax_rhs[index] = 0.0;
  Ay_rhs[index] = 0.0;
  Az_rhs[index] = 0.0;
  }
  */
}
