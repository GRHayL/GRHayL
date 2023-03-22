#include "cctk.h"
#include "IGH.h"

#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

//static inline void calculate_face_value(
//      const cGH *cctkGH,
//      const int flux_dir,
//      const double *restrict cell_var,
//      double *restrict face_var) {
//
//  const int xdir = (flux_dir == 0);
//  const int ydir = (flux_dir == 1);
//  const int zdir = (flux_dir == 2);
//
//#pragma omp parallel for
//  for(int k=cctkGH->cctk_nghostzones[2]-1; k<cctkGH->cctk_lsh[0]-(cctkGH->cctk_nghostzones[2]-2); k++)
//    for(int j=cctkGH->cctk_nghostzones[1]-1; j<cctkGH->cctk_lsh[1]-(cctkGH->cctk_nghostzones[1]-2); j++)
//      for(int i=cctkGH->cctk_nghostzones[0]-1; i<cctkGH->cctk_lsh[2]-(cctkGH->cctk_nghostzones[0]-2); i++) {
//        const int indm2  = CCTK_GFINDEX3D(cctkGH, i-2*xdir, j-2*ydir, k-2*zdir);
//        const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
//        const int index  = CCTK_GFINDEX3D(cctkGH, i,        j ,       k);
//        const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);
//
//        face_var[index] = COMPUTE_FCVAL(cell_var[indm2],
//                                        cell_var[indm1],
//                                        cell_var[index],
//                                        cell_var[indp1]);
//  }
//}

void GRHayL_IGH_calculate_tau_source_rhs(
      const cGH *restrict cctkGH,
      const eos_parameters *restrict eos,
      const double **in_metric,
      const double **in_curv,
      const double **in_prims,
      double *restrict tau_rhs) {

  const double poison = 0.0/0.0;

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++)
    for(int j=jmin; j<jmax; j++)
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

        metric_quantities metric;
        initialize_metric(in_metric[LAPSE][index],
                          in_metric[GXX][index], in_metric[GXY][index], in_metric[GXZ][index],
                          in_metric[GYY][index], in_metric[GYZ][index], in_metric[GZZ][index],
                          in_metric[BETAX][index], in_metric[BETAY][index], in_metric[BETAZ][index],
                          &metric);

        extrinsic_curvature curv;
        initialize_extrinsic_curvature(
                          in_curv[KXX][index], in_curv[KXY][index], in_curv[KXZ][index],
                          in_curv[KYY][index], in_curv[KYZ][index], in_curv[KZZ][index],
                          &curv);

        primitive_quantities prims;
        initialize_primitives(in_prims[RHOB][index], in_prims[PRESSURE][index], poison,
                              in_prims[VX][index], in_prims[VY][index], in_prims[VZ][index],
                              0.0, 0.0, 0.0,
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims);

        int speed_limited = 0;
        limit_v_and_compute_u0(eos, &metric, &prims, &speed_limited);

        conservative_quantities cons_source;
        cons_source.tau = 0;
        calculate_tau_tilde_source_term_extrinsic_curv(&prims, grhayl_eos, &metric, &curv, &cons_source);
        tau_rhs[index] += cons_source.tau;
  }
}

void GRHayL_IGH_calculate_HD_dirn_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const double *restrict dX,
      const eos_parameters *restrict eos,
      const double **in_metric,
      const double **in_prims,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      const double *restrict cmin,
      const double *restrict cmax,
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
  void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict, const primitive_quantities *restrict,
                                const eos_parameters *restrict, const metric_quantities *restrict, const double, const double,
                                conservative_quantities *restrict);
  void (*calculate_source_terms)(const primitive_quantities *restrict, const eos_parameters *restrict eos,
                                 const metric_quantities *restrict, const metric_derivatives *restrict, conservative_quantities *restrict);

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  // Set function pointer to specific function for a given direction
  switch(flux_dir) {
    case 0:
      calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn0;
      calculate_source_terms = &calculate_source_terms_dirn0;
      break;
    case 1:
      calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn1;
      calculate_source_terms = &calculate_source_terms_dirn1;
      break;
    case 2:
      calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn2;
      calculate_source_terms = &calculate_source_terms_dirn2;
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

//This section needs a lot of extra memory (all metric quantities have face GFs)
//However, computing them once and saving them would prevent the repeated computation
//of the +1 point for the derivative
//  // Calculate the values of the metric quantities at the faces by interpolating
//  // from the cell-centered quantities
//  calculate_face_value(cctkGH, flux_dir, lapse, face_lapse);
//  calculate_face_value(cctkGH, flux_dir, betax, face_betax);
//  calculate_face_value(cctkGH, flux_dir, betay, face_betay);
//  calculate_face_value(cctkGH, flux_dir, betaz, face_betaz);
//  calculate_face_value(cctkGH, flux_dir, gxx, face_gxx);
//  calculate_face_value(cctkGH, flux_dir, gxy, face_gxy);
//  calculate_face_value(cctkGH, flux_dir, gxz, face_gxz);
//  calculate_face_value(cctkGH, flux_dir, gyy, face_gyy);
//  calculate_face_value(cctkGH, flux_dir, gyz, face_gyz);
//  calculate_face_value(cctkGH, flux_dir, gzz, face_gzz);

  // This loop includes 1 ghostzone because the RHS calculation for e.g. the x direction
  // requires (i,j,k) and (i+1,j,k)
#pragma omp parallel for
  for(int k=kmin; k<kmax+1; k++)
    for(int j=jmin; j<jmax+1; j++)
      for(int i=imin; i<imax+1; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

        metric_quantities metric_face;
        GRHayL_IGH_interpolate_to_face_and_initialize_metric(
                          cctkGH, i, j, k,
                          flux_dir, in_metric[LAPSE],
                          in_metric[BETAX], in_metric[BETAY], in_metric[BETAZ],
                          in_metric[GXX], in_metric[GXY], in_metric[GXZ],
                          in_metric[GYY], in_metric[GYZ], in_metric[GZZ],
                          &metric_face);

        primitive_quantities prims_r, prims_l;
        initialize_primitives(in_prims_r[RHOB][index], in_prims_r[PRESSURE][index], poison,
                              in_prims_r[VX][index], in_prims_r[VY][index], in_prims_r[VZ][index],
                              0.0, 0.0, 0.0,
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_r);

        initialize_primitives(in_prims_l[RHOB][index], in_prims_l[PRESSURE][index], poison,
                              in_prims_l[VX][index], in_prims_l[VY][index], in_prims_l[VZ][index],
                              0.0, 0.0, 0.0,
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_l);

        int speed_limited = 0;
        limit_v_and_compute_u0(eos, &metric_face, &prims_r, &speed_limited);
        limit_v_and_compute_u0(eos, &metric_face, &prims_l, &speed_limited);

        conservative_quantities cons_fluxes;
        calculate_HLLE_fluxes(&prims_r, &prims_l, eos, &metric_face, cmin[index], cmax[index], &cons_fluxes);

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
        const int indp1 = CCTK_GFINDEX3D(cctkGH, i+xdir, j+ydir, k+zdir);

        metric_quantities metric_face, metric_facep1;
        GRHayL_IGH_interpolate_to_face_and_initialize_metric(
                          cctkGH, i, j, k,
                          flux_dir, in_metric[LAPSE],
                          in_metric[BETAX], in_metric[BETAY], in_metric[BETAZ],
                          in_metric[GXX], in_metric[GXY], in_metric[GXZ],
                          in_metric[GYY], in_metric[GYZ], in_metric[GZZ],
                          &metric_face);

        GRHayL_IGH_interpolate_to_face_and_initialize_metric(
                          cctkGH, i+xdir, j+ydir, k+zdir,
                          flux_dir, in_metric[LAPSE],
                          in_metric[BETAX], in_metric[BETAY], in_metric[BETAZ],
                          in_metric[GXX], in_metric[GXY], in_metric[GXZ],
                          in_metric[GYY], in_metric[GYZ], in_metric[GZZ],
                          &metric_facep1);

        rho_star_rhs[index] += dxi[flux_dir]*(rho_star_flux[index] - rho_star_flux[indp1]);
        tau_rhs[index]      += dxi[flux_dir]*(tau_flux[index]      - tau_flux[indp1]);
        Stildex_rhs[index]  += dxi[flux_dir]*(Stildex_flux[index]  - Stildex_flux[indp1]);
        Stildey_rhs[index]  += dxi[flux_dir]*(Stildey_flux[index]  - Stildey_flux[indp1]);
        Stildez_rhs[index]  += dxi[flux_dir]*(Stildez_flux[index]  - Stildez_flux[indp1]);

        metric_quantities metric;
        initialize_metric(in_metric[LAPSE][index],
                          in_metric[GXX][index], in_metric[GXY][index], in_metric[GXZ][index],
                          in_metric[GYY][index], in_metric[GYZ][index], in_metric[GZZ][index],
                          in_metric[BETAX][index], in_metric[BETAY][index], in_metric[BETAZ][index],
                          &metric);

        primitive_quantities prims;
        initialize_primitives(in_prims[RHOB][index], in_prims[PRESSURE][index], poison,
                              in_prims[VX][index], in_prims[VY][index], in_prims[VZ][index],
                              0.0, 0.0, 0.0,
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims);

        int speed_limited = 0;
        limit_v_and_compute_u0(eos, &metric, &prims, &speed_limited);

        metric_derivatives metric_derivs;
    
        metric_derivs.lapse[flux_dir]   = dxi[flux_dir]*(metric_facep1.lapse - metric_face.lapse);
        metric_derivs.betax[flux_dir]   = dxi[flux_dir]*(metric_facep1.betax - metric_face.betax);
        metric_derivs.betay[flux_dir]   = dxi[flux_dir]*(metric_facep1.betay - metric_face.betay);
        metric_derivs.betaz[flux_dir]   = dxi[flux_dir]*(metric_facep1.betaz - metric_face.betaz);
        metric_derivs.adm_gxx[flux_dir] = dxi[flux_dir]*(metric_facep1.adm_gxx - metric_face.adm_gxx);
        metric_derivs.adm_gxy[flux_dir] = dxi[flux_dir]*(metric_facep1.adm_gxy - metric_face.adm_gxy);
        metric_derivs.adm_gxz[flux_dir] = dxi[flux_dir]*(metric_facep1.adm_gxz - metric_face.adm_gxz);
        metric_derivs.adm_gyy[flux_dir] = dxi[flux_dir]*(metric_facep1.adm_gyy - metric_face.adm_gyy);
        metric_derivs.adm_gyz[flux_dir] = dxi[flux_dir]*(metric_facep1.adm_gyz - metric_face.adm_gyz);
        metric_derivs.adm_gzz[flux_dir] = dxi[flux_dir]*(metric_facep1.adm_gzz - metric_face.adm_gzz);

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
