#include "GRHayLMHD.h"

void GRHayLMHD_evaluate_tau_curvature_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLMHD_evaluate_tau_curvature_rhs;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_Equals(verbose, "essential+iteration output")) {
    const int levelnumber = GetRefinementLevel(cctkGH);
    CCTK_VINFO("***** Iter. # %d, Lev: %d, Integrating to time: %e *****",cctk_iteration,levelnumber,cctk_delta_time/cctk_levfac[0]+cctk_time);
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

        // Initialize RHS variables to zero
        tau_rhs[index]      = 0.0;
        rho_star_rhs[index] = 0.0;
        Stildex_rhs[index]  = 0.0;
        Stildey_rhs[index]  = 0.0;
        Stildez_rhs[index]  = 0.0;
        phitilde_rhs[index] = 0.0;
        Ax_rhs[index]       = 0.0;
        Ay_rhs[index]       = 0.0;
        Az_rhs[index]       = 0.0;

        metric_quantities ADM_metric;
        grhayl_initialize_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        extrinsic_curvature curv;
        grhayl_initialize_extrinsic_curvature(
              kxx[index], kxy[index], kxz[index],
              kyy[index], kyz[index], kzz[index],
              &curv);

        primitive_quantities prims;
        grhayl_initialize_primitives(
              rho_b[index], pressure[index], poison,
              vx[index], vy[index], vz[index],
              Bx_center[index], By_center[index], Bz_center[index],
              poison, poison, poison, // entropy, Y_e, temp
              &prims);

        int speed_limited = 0;
        grhayl_limit_v_and_compute_u0(
              grhayl_eos, &ADM_metric, &prims, &speed_limited);

        conservative_quantities cons_source;
        cons_source.tau = 0;
        calculate_tau_tilde_source_term_extrinsic_curv(
                 &prims, grhayl_eos, &ADM_metric, &curv, &cons_source);
        tau_rhs[index] += cons_source.tau;
  }
}
