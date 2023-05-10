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

#include "IGH.h"

void GRHayL_IGH_evaluate_tau_curvature_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGH_evaluate_tau_curvature_rhs;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_Equals(verbose, "essential+iteration output")) {
    const int levelnumber = GetRefinementLevel(cctkGH);
    CCTK_VINFO("***** Iter. # %d, Lev: %d, Integrating to time: %e *****",cctk_iteration,levelnumber,cctk_delta_time/cctk_levfac[0]+cctk_time);
  }

  // First initialize RHS variables to zero
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

  const double poison = 0.0/0.0;

  const int imin = cctk_nghostzones[0];
  const int jmin = cctk_nghostzones[1];
  const int kmin = cctk_nghostzones[2];
  const int imax = cctk_lsh[0] - cctk_nghostzones[0];
  const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
  const int kmax = cctk_lsh[2] - cctk_nghostzones[2];

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++)
    for(int j=jmin; j<jmax; j++)
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

        metric_quantities metric;
        initialize_metric(alp[index],
                          gxx[index], gxy[index], gxz[index],
                          gyy[index], gyz[index], gzz[index],
                          betax[index], betay[index], betaz[index],
                          &metric);

        extrinsic_curvature curv;
        initialize_extrinsic_curvature(
                          kxx[index], kxy[index], kxz[index],
                          kyy[index], kyz[index], kzz[index],
                          &curv);

        primitive_quantities prims;
        initialize_primitives(rho_b[index], pressure[index], poison,
                              vx[index], vy[index], vz[index],
                              0.0, 0.0, 0.0,
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims);

        int speed_limited = 0;
        limit_v_and_compute_u0(grhayl_eos, &metric, &prims, &speed_limited);

        conservative_quantities cons_source;
        cons_source.tau = 0;
        calculate_tau_tilde_source_term_extrinsic_curv(&prims, grhayl_eos, &metric, &curv, &cons_source);
        tau_rhs[index] += cons_source.tau;
  }
}
