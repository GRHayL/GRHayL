#include "GRHayLHD.h"

double GRHayLHD_eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in) {
  double K, Gamma;
  ghl_hybrid_get_K_and_Gamma(eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return eos->Gamma_th + (Gamma - eos->Gamma_th)*P_cold/press_in;
}

void GRHayLHD_evaluate_flux_source_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_evaluate_flux_source_rhs;
  DECLARE_CCTK_PARAMETERS;

  /*
   *  Computation of \partial_i on RHS of \partial_t {rho_star,tau,Stilde{x,y,z}},
   *  via PPM reconstruction onto e.g. (i+1/2,j,k), so that
   *  \partial_x F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   */
  const double poison = 0.0/0.0;

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];

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

  for(int flux_dir=0; flux_dir<3; flux_dir++) {
    const int xdir = (flux_dir == 0);
    const int ydir = (flux_dir == 1);
    const int zdir = (flux_dir == 2);
    const double *v_flux_dir;

    // Set function pointer to specific function for a given direction
    switch(flux_dir) {
      case 0:
        v_flux_dir = vx;
        calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
        calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn0;
        calculate_source_terms = &ghl_calculate_source_terms_dirn0;
        break;
      case 1:
        v_flux_dir = vy;
        calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
        calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn1;
        calculate_source_terms = &ghl_calculate_source_terms_dirn1;
        break;
      case 2:
        v_flux_dir = vz;
        calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
        calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn2;
        calculate_source_terms = &ghl_calculate_source_terms_dirn2;
        break;
      default:
        CCTK_VERROR("Warning: invalid flux_dir value (not 0, 1, or 2) has been passed to calculate_MHD_rhs.");
    }

    // This loop includes 1 ghostzone because the RHS calculation for e.g. the x direction
    // requires (i,j,k) and (i+1,j,k)
#pragma omp parallel for
    for(int k=kmin-1; k<kmax; k++) {
      for(int j=jmin-1; j<jmax; j++) {
        for(int i=imin-1; i<imax; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

          double rho_stencil[6], press_stencil[6], v_flux[6];
          double rhor, rhol, pressr, pressl;
          double vel_stencil[3][6], vel_r[3], vel_l[3];

          for(int ind=0; ind<6; ind++) {
            // Stencil from -2 to +3 reconstructs to e.g. i+1/2
            const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-2), j+ydir*(ind-2), k+zdir*(ind-2));
            v_flux[ind] = v_flux_dir[stencil]; // Could be smaller; doesn't use full stencil
            rho_stencil[ind] = rho_b[stencil];
            press_stencil[ind] = pressure[stencil];
            vel_stencil[0][ind] = vx[stencil];
            vel_stencil[1][ind] = vy[stencil];
            vel_stencil[2][ind] = vz[stencil];
          }

          // Compute Gamma
          const double Gamma = GRHayLHD_eos_Gamma_eff(ghl_eos, rho_b[index], pressure[index]);

          ghl_ppm(
                rho_stencil, press_stencil, vel_stencil,
                3, v_flux, Gamma,
                &rhor, &rhol, &pressr, &pressl, vel_r, vel_l);

          metric_quantities ADM_metric_face;
          GRHayLHD_interpolate_metric_to_face(
                cctkGH, i, j, k,
                flux_dir, alp,
                betax, betay, betaz,
                gxx, gxy, gxz,
                gyy, gyz, gzz,
                &ADM_metric_face);

          primitive_quantities prims_r, prims_l;
          ghl_initialize_primitives(
                rhor, pressr, poison,
                vel_r[0], vel_r[1], vel_r[2],
                0.0, 0.0, 0.0,
                poison, poison, poison, // entropy, Y_e, temp
                &prims_r);

          ghl_initialize_primitives(
                rhol, pressl, poison,
                vel_l[0], vel_l[1], vel_l[2],
                0.0, 0.0, 0.0,
                poison, poison, poison, // entropy, Y_e, temp
                &prims_l);

          int speed_limited = 0;
          ghl_limit_v_and_compute_u0(
                ghl_eos, &ADM_metric_face, &prims_r, &speed_limited);
          ghl_limit_v_and_compute_u0(
                ghl_eos, &ADM_metric_face, &prims_l, &speed_limited);

          double cmin, cmax;
          conservative_quantities cons_fluxes;
          calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin, &cmax);
          calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin, cmax, &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux[index]      = cons_fluxes.tau;
          Stildex_flux[index]  = cons_fluxes.SD[0];
          Stildey_flux[index]  = cons_fluxes.SD[1];
          Stildez_flux[index]  = cons_fluxes.SD[2];
        }
      }
    }

    CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

#pragma omp parallel for
    for(int k=kmin; k<kmax; k++) {
      for(int j=jmin; j<jmax; j++) {
        for(int i=imin; i<imax; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i-xdir, j-ydir, k-zdir);

          rho_star_rhs[index] += dxi*(rho_star_flux[indm1] - rho_star_flux[index]);
          tau_rhs[index]      += dxi*(tau_flux[indm1]      - tau_flux[index]);
          Stildex_rhs[index]  += dxi*(Stildex_flux[indm1]  - Stildex_flux[index]);
          Stildey_rhs[index]  += dxi*(Stildey_flux[indm1]  - Stildey_flux[index]);
          Stildez_rhs[index]  += dxi*(Stildez_flux[indm1]  - Stildez_flux[index]);

          metric_quantities ADM_metric;
          ghl_initialize_metric(alp[index],
                betax[index], betay[index], betaz[index],
                gxx[index], gxy[index], gxz[index],
                gyy[index], gyz[index], gzz[index],
                &ADM_metric);

          primitive_quantities prims;
          ghl_initialize_primitives(
                rho_b[index], pressure[index], eps[index],
                vx[index], vy[index], vz[index],
                0.0, 0.0, 0.0,
                poison, poison, poison, // entropy, Y_e, temp
                &prims);

          int speed_limited = 0;
          ghl_limit_v_and_compute_u0(
                ghl_eos, &ADM_metric, &prims, &speed_limited);

          metric_quantities ADM_metric_derivs;
          GRHayLHD_compute_metric_derivs(
                cctkGH, i, j, k,
                flux_dir, dxi, alp,
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
}
