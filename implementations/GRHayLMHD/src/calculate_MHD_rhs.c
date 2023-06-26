#include "GRHayLMHD.h"

CCTK_REAL eos_Gamma_eff(const CCTK_REAL rho_in, const CCTK_REAL press_in) {
  CCTK_REAL K, Gamma;
  ghl_hybrid_get_K_and_Gamma(ghl_eos, rho_in, &K, &Gamma);
  const CCTK_REAL P_cold = K*pow(rho_in, Gamma);
  return ghl_eos->Gamma_th + (Gamma - ghl_eos->Gamma_th)*P_cold/press_in;
}

// We could reduce computational cost by introducing grid functions for the
// metric face interpolations (used in both loops
void GRHayLMHD_calculate_MHD_dirn_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const CCTK_REAL *restrict dX,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      const CCTK_REAL *restrict rho_b,
      const CCTK_REAL *restrict pressure,
      const CCTK_REAL *vx,
      const CCTK_REAL *vy,
      const CCTK_REAL *vz,
      const CCTK_REAL **B_center,
      const CCTK_REAL *restrict B_stagger,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax,
      CCTK_REAL *restrict rho_star_flux,
      CCTK_REAL *restrict tau_flux,
      CCTK_REAL *restrict Stildex_flux,
      CCTK_REAL *restrict Stildey_flux,
      CCTK_REAL *restrict Stildez_flux,
      CCTK_REAL *restrict rho_star_rhs,
      CCTK_REAL *restrict tau_rhs,
      CCTK_REAL *restrict Stildex_rhs,
      CCTK_REAL *restrict Stildey_rhs,
      CCTK_REAL *restrict Stildez_rhs) {

  const CCTK_REAL dxi = 1.0/dX[flux_dir];
  const CCTK_REAL poison = 0.0/0.0;

  // Function pointer to allow for loop over fluxes and sources
  void (*calculate_characteristic_speed)(const primitive_quantities *restrict prims_r,
                                         const primitive_quantities *restrict prims_l,
                                         struct eos_parameters const *restrict eos,
                                         const metric_quantities *restrict ADM_metric_face,
                                         CCTK_REAL *cmin, CCTK_REAL *cmax);

  void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict prims_r,
                                const primitive_quantities *restrict prims_l,
                                const eos_parameters *restrict eos,
                                const metric_quantities *restrict ADM_metric_face,
                                const CCTK_REAL cmin,
                                const CCTK_REAL cmax,
                                conservative_quantities *restrict cons_fluxes);

  void (*calculate_source_terms)(const primitive_quantities *restrict prims,
                                 const eos_parameters *restrict eos,
                                 const metric_quantities *restrict ADM_metric,
                                 const metric_quantities *restrict metric_derivs,
                                 conservative_quantities *restrict cons_sources);

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const CCTK_REAL *v_flux;
  int B_recon[3];
  // Set function pointer to specific function for a given direction
  switch(flux_dir) {
    case 0:
      v_flux = vx;
      B_recon[0] = 0;
      B_recon[1] = 1;
      B_recon[2] = 2;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
      calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn0;
      calculate_source_terms = &ghl_calculate_source_terms_dirn0;
      break;
    case 1:
      v_flux = vy;
      B_recon[0] = 1;
      B_recon[1] = 2;
      B_recon[2] = 0;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
      calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn1;
      calculate_source_terms = &ghl_calculate_source_terms_dirn1;
      break;
    case 2:
      v_flux = vz;
      B_recon[0] = 2;
      B_recon[1] = 0;
      B_recon[2] = 1;
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
        const int indm1 = CCTK_GFINDEX3D(cctkGH, i-xdir, j-ydir, k-zdir); /* indexim1=0 when i=0 */
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL rho_stencil[6], press_stencil[6], v_flux_dir[6];
        CCTK_REAL rhor, rhol, pressr, pressl, B_r[3], B_l[3];
        CCTK_REAL var_data[5][6], vars_r[5], vars_l[5];

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind] = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          rho_stencil[ind] = rho_b[stencil];
          press_stencil[ind] = pressure[stencil];
          var_data[0][ind] = vx[stencil];
          var_data[1][ind] = vy[stencil];
          var_data[2][ind] = vz[stencil];
          var_data[3][ind] = B_center[B_recon[1]][stencil];
          var_data[4][ind] = B_center[B_recon[2]][stencil];
        }

        // Compute Gamma
        const CCTK_REAL Gamma = eos_Gamma_eff(rho_b[index], pressure[index]);

        ghl_ppm(
              rho_stencil, press_stencil, var_data,
              5, v_flux_dir, Gamma,
              &rhor, &rhol, &pressr, &pressl, vars_r, vars_l);

        vel_r[0][index] = vars_r[0];
        vel_r[1][index] = vars_r[1];
        vel_r[2][index] = vars_r[2];

        vel_l[0][index] = vars_l[0];
        vel_l[1][index] = vars_l[1];
        vel_l[2][index] = vars_l[2];

        B_r[B_recon[0]] = B_stagger[indm1];
        B_r[B_recon[1]] = vars_r[3];
        B_r[B_recon[2]] = vars_r[4];

        B_l[B_recon[0]] = B_stagger[indm1];
        B_l[B_recon[1]] = vars_l[3];
        B_l[B_recon[2]] = vars_l[4];

        metric_quantities ADM_metric_face;
        GRHayLMHD_interpolate_metric_to_face(
              cctkGH, i, j, k,
              flux_dir, lapse,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_face);

        primitive_quantities prims_r, prims_l;
        ghl_initialize_primitives(
              rhor, pressr, poison,
              vel_r[0][index], vel_r[1][index], vel_r[2][index],
              B_r[0], B_r[1], B_r[2],
              poison, poison, poison, // entropy, Y_e, temp
              &prims_r);

        ghl_initialize_primitives(
              rhol, pressl, poison,
              vel_l[0][index], vel_l[1][index], vel_l[2][index],
              B_l[0], B_l[1], B_l[2],
              poison, poison, poison, // entropy, Y_e, temp
              &prims_l);

        int speed_limited = 0;
        ghl_limit_v_and_compute_u0(
              ghl_eos, &ADM_metric_face, &prims_r, &speed_limited);
        ghl_limit_v_and_compute_u0(
              ghl_eos, &ADM_metric_face, &prims_l, &speed_limited);

        conservative_quantities cons_fluxes;
        calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin[index], &cmax[index]);
        calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin[index], cmax[index], &cons_fluxes);

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

        rho_star_rhs[index] += dxi*(rho_star_flux[index] - rho_star_flux[indp1]);
        tau_rhs[index]      += dxi*(tau_flux[index]      - tau_flux[indp1]);
        Stildex_rhs[index]  += dxi*(Stildex_flux[index]  - Stildex_flux[indp1]);
        Stildey_rhs[index]  += dxi*(Stildey_flux[index]  - Stildey_flux[indp1]);
        Stildez_rhs[index]  += dxi*(Stildez_flux[index]  - Stildez_flux[indp1]);

        metric_quantities ADM_metric;
        ghl_initialize_metric(
              lapse[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        primitive_quantities prims;
        ghl_initialize_primitives(
              rho_b[index], pressure[index], poison,
              vx[index], vy[index], vz[index],
              B_center[0][index], B_center[1][index], B_center[2][index],
              poison, poison, poison, // entropy, Y_e, temp
              &prims);

        int speed_limited = 0;
        ghl_limit_v_and_compute_u0(
              ghl_eos, &ADM_metric, &prims, &speed_limited);

        metric_quantities ADM_metric_derivs;
        GRHayLMHD_compute_metric_derivs(
              cctkGH, i, j, k,
              flux_dir, dxi, lapse,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs);

        conservative_quantities cons_source;
        cons_source.tau = 0.0;
        cons_source.SD[0] = 0.0;
        cons_source.SD[1] = 0.0;
        cons_source.SD[2] = 0.0;

        calculate_source_terms(&prims, ghl_eos, &ADM_metric, &ADM_metric_derivs, &cons_source);
        tau_rhs[index]     += cons_source.tau;
        Stildex_rhs[index] += cons_source.SD[0];
        Stildey_rhs[index] += cons_source.SD[1];
        Stildez_rhs[index] += cons_source.SD[2];
      }
    }
  }
}
