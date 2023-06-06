#include "GRHayLMHD.h"

// We could reduce computational cost by introducing grid functions for the
// metric face interpolations (used in both loops
void GRHayLMHD_calculate_MHD_dirn_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const double *restrict dX,
      const eos_parameters *restrict eos,
      const double **in_metric,
      const double **in_prims,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      double *restrict cmin,
      double *restrict cmax,
      double *restrict rho_star_flux,
      double *restrict tau_flux,
      double *restrict Stildex_flux,
      double *restrict Stildey_flux,
      double *restrict Stildez_flux,
      double *restrict rho_star_rhs,
      double *restrict tau_rhs,
      double *restrict Stildex_rhs,
      double *restrict Stildey_rhs,
      double *restrict Stildez_rhs) {

  const double dxi[3] = { 1.0/dX[0],1.0/dX[1],1.0/dX[2] };
  const double poison = 0.0/0.0;

  // Function pointer to allow for loop over fluxes and sources
  void (*calculate_characteristic_speed)(const primitive_quantities *restrict prims_r,
                                         const primitive_quantities *restrict prims_l,
                                         struct eos_parameters const *restrict eos,
                                         const metric_quantities *restrict ADM_metric_face,
                                         double *cmin, double *cmax);

  void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict prims_r,
                                const primitive_quantities *restrict prims_l,
                                const eos_parameters *restrict eos,
                                const metric_quantities *restrict ADM_metric_face,
                                const double cmin,
                                const double cmax,
                                conservative_quantities *restrict cons_fluxes);

  void (*calculate_source_terms)(const primitive_quantities *restrict prims,
                                 const eos_parameters *restrict eos,
                                 const metric_quantities *restrict ADM_metric,
                                 const metric_quantities *restrict metric_derivs,
                                 conservative_quantities *restrict cons_sources);

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  // Set function pointer to specific function for a given direction
  switch(flux_dir) {
    case 0:
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
      calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn0;
      calculate_source_terms = &ghl_calculate_source_terms_dirn0;
      break;
    case 1:
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
      calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn1;
      calculate_source_terms = &ghl_calculate_source_terms_dirn1;
      break;
    case 2:
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
      calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn2;
      calculate_source_terms = &ghl_calculate_source_terms_dirn2;
      break;
    default:
      CCTK_VERROR("Warning: invalid flux_dir value (not 0, 1, or 2) has been passed to calculate_MHD_rhs.");
  }

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];

  // This loop includes 1 ghostzone because the RHS calculation for e.g. the x direction
  // requires (i,j,k) and (i+1,j,k); if cmin/max weren't also needed for A_i, we could
  // technically have the loop only go 1 extra point in the flux_dir direction
#pragma omp parallel for
  for(int k=kmin; k<kmax+1; k++) {
    for(int j=jmin; j<jmax+1; j++) {
      for(int i=imin; i<imax+1; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

        metric_quantities ADM_metric_face;
        GRHayLMHD_interpolate_to_face_and_initialize_metric(
              cctkGH, i, j, k,
              flux_dir, in_metric[LAPSE],
              in_metric[BETAX], in_metric[BETAY], in_metric[BETAZ],
              in_metric[GXX], in_metric[GXY], in_metric[GXZ],
              in_metric[GYY], in_metric[GYZ], in_metric[GZZ],
              &ADM_metric_face);

        primitive_quantities prims_r, prims_l;
        ghl_initialize_primitives(
              in_prims_r[RHOB][index], in_prims_r[PRESSURE][index], poison,
              in_prims_r[VX][index], in_prims_r[VY][index], in_prims_r[VZ][index],
              in_prims_r[BX_CENTER][index], in_prims_r[BY_CENTER][index], in_prims_r[BZ_CENTER][index],
              poison, poison, poison, // entropy, Y_e, temp
              &prims_r);

        ghl_initialize_primitives(
              in_prims_l[RHOB][index], in_prims_l[PRESSURE][index], poison,
              in_prims_l[VX][index], in_prims_l[VY][index], in_prims_l[VZ][index],
              in_prims_l[BX_CENTER][index], in_prims_l[BY_CENTER][index], in_prims_l[BZ_CENTER][index],
              poison, poison, poison, // entropy, Y_e, temp
              &prims_l);

        int speed_limited = 0;
        ghl_limit_v_and_compute_u0(
              eos, &ADM_metric_face, &prims_r, &speed_limited);
        ghl_limit_v_and_compute_u0(
              eos, &ADM_metric_face, &prims_l, &speed_limited);

        conservative_quantities cons_fluxes;
        calculate_characteristic_speed(&prims_r, &prims_l, eos, &ADM_metric_face, &cmin[index], &cmax[index]);
        calculate_HLLE_fluxes(&prims_r, &prims_l, eos, &ADM_metric_face, cmin[index], cmax[index], &cons_fluxes);

        rho_star_flux[index] = cons_fluxes.rho;
        tau_flux[index]      = cons_fluxes.tau;
        Stildex_flux[index]  = cons_fluxes.SD[0];
        Stildey_flux[index]  = cons_fluxes.SD[1];
        Stildez_flux[index]  = cons_fluxes.SD[2];
      }
    }
  }

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);
        const int indp1 = CCTK_GFINDEX3D(cctkGH, i+xdir, j+ydir, k+zdir);

        metric_quantities ADM_metric_face, ADM_metric_facep1;
        GRHayLMHD_interpolate_to_face_and_initialize_metric(
              cctkGH, i, j, k,
              flux_dir, in_metric[LAPSE],
              in_metric[BETAX], in_metric[BETAY], in_metric[BETAZ],
              in_metric[GXX], in_metric[GXY], in_metric[GXZ],
              in_metric[GYY], in_metric[GYZ], in_metric[GZZ],
              &ADM_metric_face);

        GRHayLMHD_interpolate_to_face_and_initialize_metric(
              cctkGH, i+xdir, j+ydir, k+zdir,
              flux_dir, in_metric[LAPSE],
              in_metric[BETAX], in_metric[BETAY], in_metric[BETAZ],
              in_metric[GXX], in_metric[GXY], in_metric[GXZ],
              in_metric[GYY], in_metric[GYZ], in_metric[GZZ],
              &ADM_metric_facep1);

        rho_star_rhs[index] += dxi[flux_dir]*(rho_star_flux[index] - rho_star_flux[indp1]);
        tau_rhs[index]      += dxi[flux_dir]*(tau_flux[index]      - tau_flux[indp1]);
        Stildex_rhs[index]  += dxi[flux_dir]*(Stildex_flux[index]  - Stildex_flux[indp1]);
        Stildey_rhs[index]  += dxi[flux_dir]*(Stildey_flux[index]  - Stildey_flux[indp1]);
        Stildez_rhs[index]  += dxi[flux_dir]*(Stildez_flux[index]  - Stildez_flux[indp1]);

        metric_quantities ADM_metric;
        ghl_initialize_metric(
              in_metric[LAPSE][index],
              in_metric[BETAX][index], in_metric[BETAY][index], in_metric[BETAZ][index],
              in_metric[GXX][index], in_metric[GXY][index], in_metric[GXZ][index],
              in_metric[GYY][index], in_metric[GYZ][index], in_metric[GZZ][index],
              &ADM_metric);

        primitive_quantities prims;
        ghl_initialize_primitives(
              in_prims[RHOB][index], in_prims[PRESSURE][index], poison,
              in_prims[VX][index], in_prims[VY][index], in_prims[VZ][index],
              in_prims[BX_CENTER][index], in_prims[BY_CENTER][index], in_prims[BZ_CENTER][index],
              poison, poison, poison, // entropy, Y_e, temp
              &prims);

        int speed_limited = 0;
        ghl_limit_v_and_compute_u0(
              eos, &ADM_metric, &prims, &speed_limited);

        metric_quantities ADM_metric_derivs;

        ADM_metric_derivs.lapse      = dxi[flux_dir]*(ADM_metric_facep1.lapse - ADM_metric_face.lapse);
        ADM_metric_derivs.betaU[0]   = dxi[flux_dir]*(ADM_metric_facep1.betaU[0] - ADM_metric_face.betaU[0]);
        ADM_metric_derivs.betaU[1]   = dxi[flux_dir]*(ADM_metric_facep1.betaU[1] - ADM_metric_face.betaU[1]);
        ADM_metric_derivs.betaU[2]   = dxi[flux_dir]*(ADM_metric_facep1.betaU[2] - ADM_metric_face.betaU[2]);

        ADM_metric_derivs.gammaDD[0][0] = dxi[flux_dir]*(ADM_metric_facep1.gammaDD[0][0] - ADM_metric_face.gammaDD[0][0]);
        ADM_metric_derivs.gammaDD[0][1] = dxi[flux_dir]*(ADM_metric_facep1.gammaDD[0][1] - ADM_metric_face.gammaDD[0][1]);
        ADM_metric_derivs.gammaDD[0][2] = dxi[flux_dir]*(ADM_metric_facep1.gammaDD[0][2] - ADM_metric_face.gammaDD[0][2]);
        ADM_metric_derivs.gammaDD[1][1] = dxi[flux_dir]*(ADM_metric_facep1.gammaDD[1][1] - ADM_metric_face.gammaDD[1][1]);
        ADM_metric_derivs.gammaDD[1][2] = dxi[flux_dir]*(ADM_metric_facep1.gammaDD[1][2] - ADM_metric_face.gammaDD[1][2]);
        ADM_metric_derivs.gammaDD[2][2] = dxi[flux_dir]*(ADM_metric_facep1.gammaDD[2][2] - ADM_metric_face.gammaDD[2][2]);

        conservative_quantities cons_source;
        cons_source.tau = 0.0;
        cons_source.SD[0] = 0.0;
        cons_source.SD[1] = 0.0;
        cons_source.SD[2] = 0.0;

        calculate_source_terms(&prims, grhayl_eos, &ADM_metric, &ADM_metric_derivs, &cons_source);
        tau_rhs[index]     += cons_source.tau;
        Stildex_rhs[index] += cons_source.SD[0];
        Stildey_rhs[index] += cons_source.SD[1];
        Stildez_rhs[index] += cons_source.SD[2];
      }
    }
  }
}
