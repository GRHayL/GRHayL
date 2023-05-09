#include "IGH.h"

double GRHayL_IGH_eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in) {
  double K, Gamma;
  eos->hybrid_get_K_and_Gamma(eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return eos->Gamma_th + (Gamma - eos->Gamma_th)*P_cold/press_in;
}

void GRHayL_IGH_evaluate_flux_source_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGH_evaluate_flux_source_rhs;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dxi[3] = { 1.0/CCTK_DELTA_SPACE(0), 1.0/CCTK_DELTA_SPACE(1), 1.0/CCTK_DELTA_SPACE(2) };

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


  for(int flux_dir=0; flux_dir<3; flux_dir++) {
    const int xdir = (flux_dir == 0);
    const int ydir = (flux_dir == 1);
    const int zdir = (flux_dir == 2);
    const double *v_flux_dir;

    void (*calculate_characteristic_speed)(const primitive_quantities *restrict prims_r,
                                           const primitive_quantities *restrict prims_l,
                                           struct eos_parameters const *restrict eos,
                                           const metric_quantities *restrict metric_face,
                                           double *cmin, double *cmax);

    void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict, const primitive_quantities *restrict,
                                  const eos_parameters *restrict, const metric_quantities *restrict, const double, const double,
                                  conservative_quantities *restrict);

    void (*calculate_source_terms)(const primitive_quantities *restrict, const eos_parameters *restrict eos,
                                   const metric_quantities *restrict, const metric_derivatives *restrict, conservative_quantities *restrict);

    // Set function pointer to specific function for a given direction
    switch(flux_dir) {
      case 0:
        v_flux_dir = vx;
        calculate_characteristic_speed = &calculate_characteristic_speed_dirn0;
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn0;
        calculate_source_terms = &calculate_source_terms_dirn0;
        break;
      case 1:
        v_flux_dir = vy;
        calculate_characteristic_speed = &calculate_characteristic_speed_dirn1;
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn1;
        calculate_source_terms = &calculate_source_terms_dirn1;
        break;
      case 2:
        v_flux_dir = vz;
        calculate_characteristic_speed = &calculate_characteristic_speed_dirn2;
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn2;
        calculate_source_terms = &calculate_source_terms_dirn2;
        break;
      default:
        CCTK_VERROR("Warning: invalid flux_dir value (not 0, 1, or 2) has been passed to calculate_MHD_rhs.");
    }

    // This loop includes 1 ghostzone because the RHS calculation for e.g. the x direction
    // requires (i,j,k) and (i+1,j,k)
#pragma omp parallel for
    for(int k=kmin-1; k<kmax; k++)
      for(int j=jmin-1; j<jmax; j++)
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
          const double Gamma = GRHayL_IGH_eos_Gamma_eff(grhayl_eos, rho_b[index], pressure[index]);

          simple_ppm(
            rho_stencil, press_stencil, vel_stencil, 3,
            v_flux, Gamma,
            &rhor, &rhol, &pressr, &pressl, vel_r, vel_l);

          metric_quantities metric_face;
          GRHayL_IGH_interpolate_to_face_and_initialize_metric(
                            cctkGH, i, j, k,
                            flux_dir, alp,
                            betax, betay, betaz,
                            gxx, gxy, gxz,
                            gyy, gyz, gzz,
                            &metric_face);

          primitive_quantities prims_r, prims_l;
          initialize_primitives(rhor, pressr, poison,
                                vel_r[0], vel_r[1], vel_r[2],
                                0.0, 0.0, 0.0,
                                poison, poison, poison, // entropy, Y_e, temp
                                &prims_r);

          initialize_primitives(rhol, pressl, poison,
                                vel_l[0], vel_l[1], vel_l[2],
                                0.0, 0.0, 0.0,
                                poison, poison, poison, // entropy, Y_e, temp
                                &prims_l);

          int speed_limited = 0;
          limit_v_and_compute_u0(grhayl_eos, &metric_face, &prims_r, &speed_limited);
          limit_v_and_compute_u0(grhayl_eos, &metric_face, &prims_l, &speed_limited);

          double cmin, cmax;
          conservative_quantities cons_fluxes;
          calculate_characteristic_speed(&prims_r, &prims_l, grhayl_eos, &metric_face, &cmin, &cmax);
          calculate_HLLE_fluxes(&prims_r, &prims_l, grhayl_eos, &metric_face, cmin, cmax, &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux[index]      = cons_fluxes.tau;
          Stildex_flux[index]  = cons_fluxes.S_x;
          Stildey_flux[index]  = cons_fluxes.S_y;
          Stildez_flux[index]  = cons_fluxes.S_z;
    }

#pragma omp parallel for
    for(int k=kmin; k<kmax; k++)
      for(int j=jmin; j<jmax; j++)
        for(int i=imin; i<imax; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i-xdir, j-ydir, k-zdir);

          metric_quantities metric_face, metric_facem1;
          GRHayL_IGH_interpolate_to_face_and_initialize_metric(
                            cctkGH, i, j, k,
                            flux_dir, alp,
                            betax, betay, betaz,
                            gxx, gxy, gxz,
                            gyy, gyz, gzz,
                            &metric_face);

          GRHayL_IGH_interpolate_to_face_and_initialize_metric(
                            cctkGH, i-xdir, j-ydir, k-zdir,
                            flux_dir, alp,
                            betax, betay, betaz,
                            gxx, gxy, gxz,
                            gyy, gyz, gzz,
                            &metric_facem1);

          rho_star_rhs[index] += dxi[flux_dir]*(rho_star_flux[indm1] - rho_star_flux[index]);
          tau_rhs[index]      += dxi[flux_dir]*(tau_flux[indm1]      - tau_flux[index]);
          Stildex_rhs[index]  += dxi[flux_dir]*(Stildex_flux[indm1]  - Stildex_flux[index]);
          Stildey_rhs[index]  += dxi[flux_dir]*(Stildey_flux[indm1]  - Stildey_flux[index]);
          Stildez_rhs[index]  += dxi[flux_dir]*(Stildez_flux[indm1]  - Stildez_flux[index]);

          metric_quantities metric;
          initialize_metric(alp[index],
                            gxx[index], gxy[index], gxz[index],
                            gyy[index], gyz[index], gzz[index],
                            betax[index], betay[index], betaz[index],
                            &metric);

          primitive_quantities prims;
          initialize_primitives(rho_b[index], pressure[index], poison,
                                vx[index], vy[index], vz[index],
                                0.0, 0.0, 0.0,
                                poison, poison, poison, // entropy, Y_e, temp
                                &prims);

          int speed_limited = 0;
          limit_v_and_compute_u0(grhayl_eos, &metric, &prims, &speed_limited);

          metric_derivatives metric_derivs;
      
          metric_derivs.lapse[flux_dir]   = dxi[flux_dir]*(metric_face.lapse - metric_facem1.lapse);
          metric_derivs.betax[flux_dir]   = dxi[flux_dir]*(metric_face.betax - metric_facem1.betax);
          metric_derivs.betay[flux_dir]   = dxi[flux_dir]*(metric_face.betay - metric_facem1.betay);
          metric_derivs.betaz[flux_dir]   = dxi[flux_dir]*(metric_face.betaz - metric_facem1.betaz);
          metric_derivs.adm_gxx[flux_dir] = dxi[flux_dir]*(metric_face.adm_gxx - metric_facem1.adm_gxx);
          metric_derivs.adm_gxy[flux_dir] = dxi[flux_dir]*(metric_face.adm_gxy - metric_facem1.adm_gxy);
          metric_derivs.adm_gxz[flux_dir] = dxi[flux_dir]*(metric_face.adm_gxz - metric_facem1.adm_gxz);
          metric_derivs.adm_gyy[flux_dir] = dxi[flux_dir]*(metric_face.adm_gyy - metric_facem1.adm_gyy);
          metric_derivs.adm_gyz[flux_dir] = dxi[flux_dir]*(metric_face.adm_gyz - metric_facem1.adm_gyz);
          metric_derivs.adm_gzz[flux_dir] = dxi[flux_dir]*(metric_face.adm_gzz - metric_facem1.adm_gzz);

          conservative_quantities cons_source;
          cons_source.tau = 0.0;
          cons_source.S_x = 0.0;
          cons_source.S_y = 0.0;
          cons_source.S_z = 0.0;

          calculate_source_terms(&prims, grhayl_eos, &metric, &metric_derivs, &cons_source);
          tau_rhs[index]     += cons_source.tau;
          Stildex_rhs[index] += cons_source.S_x;
          Stildey_rhs[index] += cons_source.S_y;
          Stildez_rhs[index] += cons_source.S_z;
    }
  }
}
