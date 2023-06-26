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

#include "GRHayLMHD.h"

void GRHayLMHD_evaluate_MHD_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_evaluate_MHD_rhs;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dX[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };

  // in_prims,out_prims_r, and out_prims_l are arrays of pointers to the actual gridfunctions.
  // Most pointers are passed explicitly. However, we need to programmatically choose gridfunctions
  // for the A_i reconstructions.
  const double *in_prims[MAXNUMVARS];
  double *out_prims_r[MAXNUMVARS];
  double *out_prims_l[MAXNUMVARS];
  in_prims[BX_STAGGER]=Bx_stagger; out_prims_r[BX_STAGGER]=Bx_staggerr; out_prims_l[BX_STAGGER]=Bx_staggerl;
  in_prims[BY_STAGGER]=By_stagger; out_prims_r[BY_STAGGER]=By_staggerr; out_prims_l[BY_STAGGER]=By_staggerl;
  in_prims[BZ_STAGGER]=Bz_stagger; out_prims_r[BZ_STAGGER]=Bz_staggerr; out_prims_l[BZ_STAGGER]=Bz_staggerl;
  in_prims[VXR       ]=vxr;        out_prims_r[VXR       ]=vxrr;        out_prims_l[VXR       ]=vxrl;
  in_prims[VYR       ]=vyr;        out_prims_r[VYR       ]=vyrr;        out_prims_l[VYR       ]=vyrl;
  in_prims[VZR       ]=vzr;        out_prims_r[VZR       ]=vzrr;        out_prims_l[VZR       ]=vzrl;
  in_prims[VXL       ]=vxl;        out_prims_r[VXL       ]=vxlr;        out_prims_l[VXL       ]=vxll;
  in_prims[VYL       ]=vyl;        out_prims_r[VYL       ]=vylr;        out_prims_l[VYL       ]=vyll;
  in_prims[VZL       ]=vzl;        out_prims_r[VZL       ]=vzlr;        out_prims_l[VZL       ]=vzll;

  const double *vel[3] = {vx, vy, vz};
  const double *B_center[3] = {Bx_center, By_center, Bz_center};
  double *vel_r[3] = {vxr, vyr, vzr};
  double *vel_l[3] = {vxl, vyl, vzl};
  const double *B_stagger[3] = {Bx_stagger, By_stagger, Bz_stagger};
  double *cmin[3] = {cmin_x, cmin_y, cmin_z};
  double *cmax[3] = {cmax_x, cmax_y, cmax_z};

  /*
     We first compute the x-direction RHS for the conservative hydro variables using PPM
     reconstruction from (i, j, k) to (i-1/2, j, k). This reconstructs rho_b, P, v^i, and
     B^i. B^x is not reconstructed, as this is given by the Bx_stagger grid function. This
     is all done to be able to compute
         \partial_x F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
     We also save the v^i reconstructions for the A_i RHS.
  */
  int flux_dir = 0;
  GRHayLMHD_calculate_MHD_dirn_rhs(cctkGH, flux_dir, dX,
                         alp,  betax,  betay,  betaz,  gxx,  gxy,  gxz,  gyy,  gyz,  gzz,
                         rho_b, pressure, vx, vy, vz, B_center, B_stagger[flux_dir],
                         vel_r, vel_l, cmin[flux_dir], cmax[flux_dir],
                         rho_star_flux, tau_flux, Stildex_flux, Stildey_flux, Stildez_flux,
                         rho_star_rhs, tau_rhs, Stildex_rhs, Stildey_rhs, Stildez_rhs);

  /*
     Here we perform reconstructions in preparation for computing the A_i RHS. First, we aim to compute
         \partial_t A_z - [gauge terms] = \psi^{6} (v^x B^y - v^y B^x)
     where A_z is defined at (i+1/2, j+1/2, k). For this, we need vx, vy, Bx, and By. Since
     we have B_stagger, the main issue is the double reconstructed velocities. However, we
     saved the previous reconstructions for this. So, the following reconstructs
         Bx_stagger @ (i+1/2, j, k) -> (i+1/2, j-1/2, k)
         By_stagger @ (i, j+1/2, k) -> (i-1/2, j+1/2, k)
         v{x,y}_{r,l} @ (i-1/2, j, k) -> (i-1/2, j-1/2, k)

      We also perform the reconstruction
         Bz_stagger @ (i, j, k+1/2) -> (i, j-1/2, k+1/2)
      in preparation for computing the RHS of A_x.
  */
  {
    const int var_indices[1] = {BY_STAGGER};
    GRHayLMHD_reconstruction_loop(cctkGH, flux_dir, 1, var_indices, rho_b, pressure, vel[flux_dir], in_prims, out_prims_r, out_prims_l);
  }

  flux_dir=1;

  {
    const int var_indices[6] = {VXR, VYR, VXL, VYL, BX_STAGGER, BZ_STAGGER};
    GRHayLMHD_reconstruction_loop(cctkGH, flux_dir, 6, var_indices, rho_b, pressure, vel[flux_dir], in_prims, out_prims_r, out_prims_l);
  }

  /*
     We compute the y-direction RHS for the conservative hydro variables using PPM
     reconstruction from (i, j, k) to (i, j-1/2, k). This reconstructs rho_b, P, v^i, and
     B^i. B^x is not reconstructed, as this is given by the Bx_stagger grid function. This
     is all done to be able to compute
         \partial_x F = [ F(i,j+1/2,k) - F(i,j-1/2,k) ] / dx
     We also save the v^i reconstructions for the A_i RHS.
  */
  GRHayLMHD_calculate_MHD_dirn_rhs(cctkGH, flux_dir, dX,
                         alp,  betax,  betay,  betaz,  gxx,  gxy,  gxz,  gyy,  gyz,  gzz,
                         rho_b, pressure, vx, vy, vz, B_center, B_stagger[flux_dir],
                         vel_r, vel_l, cmin[flux_dir], cmax[flux_dir],
                         rho_star_flux, tau_flux, Stildex_flux, Stildey_flux, Stildez_flux,
                         rho_star_rhs, tau_rhs, Stildex_rhs, Stildey_rhs, Stildez_rhs);

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
  flux_dir=2;
  GRHayLMHD_A_flux_rhs(cctkGH, flux_dir, out_prims_r, out_prims_l, phi_bssn, cmin, cmax, Az_rhs);

  /*
     Here we perform reconstructions in preparation for computing the A_i RHS. We aim to compute
         \partial_t A_x - [gauge terms] = \psi^{6} (v^y B^z - v^z B^y)
     where A_x is defined at (i, j+1/2, k+1/2). For this, we need vy, vz, By, and Bz. We previously reconstructed
         Bz_stagger @ (i, j, k+1/2) -> (i, j-1/2, k+1/2)

     The following reconstructs
         By_stagger @ (i, j+1/2, k) -> (i, j+1/2, k-1/2)
         v{y,z}_{r,l} @ (i-1/2, j, k) -> (i-1/2, j-1/2, k)

      We also perform the reconstruction
         Bx_stagger @ (i+1/2, j, k) -> (i+1/2, j, k-1/2)
      in preparation for computing the RHS of A_y.
  */
  {
    const int var_indices[6] = {VYR, VZR, VYL, VZL, BX_STAGGER, BY_STAGGER};
    GRHayLMHD_reconstruction_loop(cctkGH, flux_dir, 6, var_indices, rho_b, pressure, vel[flux_dir], in_prims, out_prims_r, out_prims_l);
  }

  /*
     We compute the z-direction RHS for the conservative hydro variables using PPM
     reconstruction from (i, j, k) to (i, j-1/2, k). This reconstructs rho_b, P, v^i, and
     B^i. B^x is not reconstructed, as this is given by the Bx_stagger grid function. This
     is all done to be able to compute
       \partial_z F = [ F(i,j,k+1/2) - F(i,j,k-1/2) ] / dx
     We also save the v^i reconstructions for the A_i RHS.
  */
  GRHayLMHD_calculate_MHD_dirn_rhs(cctkGH, flux_dir, dX,
                         alp,  betax,  betay,  betaz,  gxx,  gxy,  gxz,  gyy,  gyz,  gzz,
                         rho_b, pressure, vx, vy, vz, B_center, B_stagger[flux_dir],
                         vel_r, vel_l, cmin[flux_dir], cmax[flux_dir],
                         rho_star_flux, tau_flux, Stildex_flux, Stildey_flux, Stildez_flux,
                         rho_star_rhs, tau_rhs, Stildex_rhs, Stildey_rhs, Stildez_rhs);

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
  flux_dir=0;
  // cmin/max could be done internally using the same indices as v and B if all c var pointers were collected into
  // single array
  GRHayLMHD_A_flux_rhs(cctkGH, flux_dir, out_prims_r, out_prims_l, phi_bssn, cmin, cmax, Ax_rhs);

  /*
     Here we perform reconstructions in preparation for computing the A_i RHS. We aim to compute
         \partial_t A_y - [gauge terms] = \psi^{6} (v^z B^x - v^x B^z)
     where A_y is defined at (i+1/2, j, k+1/2). For this, we need vx, vz, Bx, and Bz. We previously reconstructed
         Bx_stagger @ (i+1/2, j, k) -> (i+1/2, j, k-1/2)
     The following reconstructs
         Bz_stagger @ (i, j, k+1/2) -> (i-1/2, j, k+1/2)
         v{x,z}_{r,l} @ (i, j, k-1/2) -> (i-1/2, j, k-1/2)
  */
  {
    const int var_indices[5] = {VXR, VZR, VXL, VZL, BZ_STAGGER};
    GRHayLMHD_reconstruction_loop(cctkGH, flux_dir, 5, var_indices, rho_b, pressure, vel[flux_dir], in_prims, out_prims_r, out_prims_l);
  }


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
  flux_dir=1;
  GRHayLMHD_A_flux_rhs(cctkGH, flux_dir, out_prims_r, out_prims_l, phi_bssn, cmin, cmax, Ay_rhs);
}
