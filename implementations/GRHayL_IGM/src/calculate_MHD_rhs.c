#include "cctk.h"
#include "IGM.h"

//#define AM2 -0.0625
//#define AM1  0.5625
//#define A0   0.5625
//#define A1  -0.0625
//#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))
//
//static inline void calculate_face_value(
//      const cGH *cctkGH,
//      const int flux_dirn,
//      const double *restrict cell_var,
//      double *restrict face_var) {
//
//  const int xdir = (flux_dirn == 0);
//  const int ydir = (flux_dirn == 1);
//  const int zdir = (flux_dirn == 2);
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

void calculate_MHD_rhs(const cGH *cctkGH, const int flux_dirn, double *restrict dX,
                       const eos_parameters *restrict eos,
                       const double **metric,
                       const double **in_prims,
                       /*const*/ double **in_prims_r,
                       /*const*/ double **in_prims_l,
                       double *restrict rho_star_flux, double *restrict tau_flux,
                       double *restrict Stildex_flux, double *restrict Stildey_flux, double *restrict Stildez_flux,
                       double *restrict rho_star_rhs, double *restrict tau_rhs,
                       double *restrict Stildex_rhs, double *restrict Stildey_rhs, double *restrict Stildez_rhs) {

  const double dxi[3] = { 1.0/dX[0],1.0/dX[1],1.0/dX[2] };
  const double poison = 0.0/0.0;

  // Function pointer to allow for loop over fluxes
  void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict, const primitive_quantities *restrict,
                              const eos_parameters *restrict, const metric_quantities *restrict, conservative_quantities *restrict);

  const int xdir = (flux_dirn == 0);
  const int ydir = (flux_dirn == 1);
  const int zdir = (flux_dirn == 2);

  // Set function pointer to specific function for a given direction
  switch(flux_dirn) {
    case 0:
      calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn0;
      break;
    case 1:
      calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn1;
      break;
    case 2:
      calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn2;
      break;
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
//  calculate_face_value(cctkGH, flux_dirn, lapse, face_lapse);
//  calculate_face_value(cctkGH, flux_dirn, betax, face_betax);
//  calculate_face_value(cctkGH, flux_dirn, betay, face_betay);
//  calculate_face_value(cctkGH, flux_dirn, betaz, face_betaz);
//  calculate_face_value(cctkGH, flux_dirn, gxx, face_gxx);
//  calculate_face_value(cctkGH, flux_dirn, gxy, face_gxy);
//  calculate_face_value(cctkGH, flux_dirn, gxz, face_gxz);
//  calculate_face_value(cctkGH, flux_dirn, gyy, face_gyy);
//  calculate_face_value(cctkGH, flux_dirn, gyz, face_gyz);
//  calculate_face_value(cctkGH, flux_dirn, gzz, face_gzz);
//
//  // Upper bound includes 1 ghostzone for RHS calculation in following
//  // loop.
//#pragma omp parallel for
//  for(int k=kmin; k<kmax+1; k++)
//    for(int j=jmin; j<jmax+1; j++)
//      for(int i=imin; i<imax+1; i++) {
//        const int index  = CCTK_GFINDEX3D(cctkGH, i, j ,k);
//
//        metric_quantities metric_face;
//        initialize_metric(face_lapse[index],
//                          face_gxx[index], face_gxy[index], face_gxz[index],
//                          face_gyy[index], face_gyz[index], face_gzz[index],
//                          face_betax[index], face_betay[index], face_betaz[index],
//                          &metric_face);
//
//        primitive_quantities prims_r, prims_l;
//        initialize_primitives(rho_r[index], press_r[index], poison,
//                              vx_r[index], vy_r[index], vz_r[index],
//                              Bx_r[index], By_r[index], Bz_r[index],
//                              poison, poison, poison, // entropy, Y_e, temp
//                              &prims_r);
//
//        initialize_primitives(rho_l[index], press_l[index], poison,
//                              vx_l[index], vy_l[index], vz_l[index],
//                              Bx_l[index], By_l[index], Bz_l[index],
//                              poison, poison, poison, // entropy, Y_e, temp
//                              &prims_l);
//
//        // Generate randomized u^0
//        prims_r.u0 = rho_r[index]*Bx_r[index]/vy_r[index];
//        prims_l.u0 = rho_l[index]*Bx_l[index]/vy_l[index];
//
//        conservative_quantities cons_fluxes;
//        calculate_HLLE_fluxes(&prims_r, 
//                              &prims_l,
//                              &eos,
//                              &metric_face, 
//                              &cons_fluxes);
//
//        rho_star_flux[index]  = cons_fluxes.rho;
//        tau_flux[index]       = cons_fluxes.tau;
//        Stildex_flux[index]       = cons_fluxes.S_x;
//        Stildey_flux[index]       = cons_fluxes.S_y;
//        Stildez_flux[index]       = cons_fluxes.S_z;
//  }
//
//#pragma omp parallel for
//  for(int k=kmin; k<kmax; k++)
//    for(int j=jmin; j<jmax; j++)
//      for(int i=imin; i<imax; i++) {
//        const int index  = CCTK_GFINDEX3D(cctkGH, i, j ,k);
//        const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir, j+ydir, k+zdir);
//
//        rho_star_rhs[index] += dxi[flux_dirn]*(rho_star_flux[index] - rho_star_flux[indp1]);
//        tau_rhs[index] += dxi[flux_dirn]*(tau_flux[index] - tau_flux[indp1]);
//        Stildex_rhs[index] += dxi[flux_dirn]*(Stildex_flux[index] - Stildex_flux[indp1]);
//        Stildey_rhs[index] += dxi[flux_dirn]*(Stildey_flux[index] - Stildey_flux[indp1]);
//        Stildez_rhs[index] += dxi[flux_dirn]*(Stildez_flux[index] - Stildez_flux[indp1]);
//
//        D_lapse[flux_dirn][index] = dxi[flux_dirn]*(face_lapse[indp1] - face_lapse[index]);
//        D_betax[flux_dirn][index] = dxi[flux_dirn]*(face_betax[indp1] - face_betax[index]);
//        D_betay[flux_dirn][index] = dxi[flux_dirn]*(face_betay[indp1] - face_betay[index]);
//        D_betaz[flux_dirn][index] = dxi[flux_dirn]*(face_betaz[indp1] - face_betaz[index]);
//        D_gxx[flux_dirn][index] = dxi[flux_dirn]*(face_gxx[indp1] - face_gxx[index]);
//        D_gxy[flux_dirn][index] = dxi[flux_dirn]*(face_gxy[indp1] - face_gxy[index]);
//        D_gxz[flux_dirn][index] = dxi[flux_dirn]*(face_gxz[indp1] - face_gxz[index]);
//        D_gyy[flux_dirn][index] = dxi[flux_dirn]*(face_gyy[indp1] - face_gyy[index]);
//        D_gyz[flux_dirn][index] = dxi[flux_dirn]*(face_gyz[indp1] - face_gyz[index]);
//        D_gzz[flux_dirn][index] = dxi[flux_dirn]*(face_gzz[indp1] - face_gzz[index]);
//  }

#pragma omp parallel for
  for(int k=kmin; k<kmax+1; k++)
    for(int j=jmin; j<jmax+1; j++)
      for(int i=imin; i<imax+1; i++) {
        const int index  = CCTK_GFINDEX3D(cctkGH, i, j ,k);

        metric_quantities metric_face;
        interpolate_to_face_and_initialize_metric(
                          cctkGH, i, j, k,
                          flux_dirn, metric[LAPSE],
                          metric[BETAX], metric[BETAY], metric[BETAZ],
                          metric[GXX], metric[GXY], metric[GXZ],
                          metric[GYY], metric[GYZ], metric[GZZ],
                          &metric_face);

        primitive_quantities prims_r, prims_l;
        initialize_primitives(in_prims_r[RHOB][index], in_prims_r[PRESSURE][index], poison,
                              in_prims_r[VX][index], in_prims_r[VY][index], in_prims_r[VZ][index],
                              in_prims_r[BX_CENTER][index], in_prims_r[BY_CENTER][index], in_prims_r[BZ_CENTER][index],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_r);

        initialize_primitives(in_prims_l[RHOB][index], in_prims_l[PRESSURE][index], poison,
                              in_prims_l[VX][index], in_prims_l[VY][index], in_prims_l[VZ][index],
                              in_prims_l[BX_CENTER][index], in_prims_l[BY_CENTER][index], in_prims_l[BZ_CENTER][index],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims_l);

        int speed_limit = 0;
        limit_v_and_compute_u0(eos, &metric_face, &prims_r, &speed_limit);
        limit_v_and_compute_u0(eos, &metric_face, &prims_l, &speed_limit);

        conservative_quantities cons_fluxes;
        calculate_HLLE_fluxes(&prims_r, &prims_l, eos, &metric_face, &cons_fluxes);

        rho_star_flux[index]  = cons_fluxes.rho;
        tau_flux[index]       = cons_fluxes.tau;
        Stildex_flux[index]       = cons_fluxes.S_x;
        Stildey_flux[index]       = cons_fluxes.S_y;
        Stildez_flux[index]       = cons_fluxes.S_z;
  }

  for(int k=kmin; k<kmax; k++)
    for(int j=jmin; j<jmax; j++)
      for(int i=imin; i<imax; i++) {
        const int index  = CCTK_GFINDEX3D(cctkGH, i, j ,k);
        const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir, j+ydir, k+zdir);

        metric_quantities metric_face, metric_facep1;
        interpolate_to_face_and_initialize_metric(
                          cctkGH, i, j, k,
                          flux_dirn, metric[LAPSE],
                          metric[BETAX], metric[BETAY], metric[BETAZ],
                          metric[GXX], metric[GXY], metric[GXZ],
                          metric[GYY], metric[GYZ], metric[GZZ],
                          &metric_face);

        interpolate_to_face_and_initialize_metric(
                          cctkGH, i+xdir, j+ydir, k+zdir,
                          flux_dirn, metric[LAPSE],
                          metric[BETAX], metric[BETAY], metric[BETAZ],
                          metric[GXX], metric[GXY], metric[GXZ],
                          metric[GYY], metric[GYZ], metric[GZZ],
                          &metric_facep1);

        rho_star_rhs[index] += dxi[flux_dirn]*(rho_star_flux[index] - rho_star_flux[indp1]);
        tau_rhs[index] += dxi[flux_dirn]*(tau_flux[index] - tau_flux[indp1]);
        Stildex_rhs[index] += dxi[flux_dirn]*(Stildex_flux[index] - Stildex_flux[indp1]);
        Stildey_rhs[index] += dxi[flux_dirn]*(Stildey_flux[index] - Stildey_flux[indp1]);
        Stildez_rhs[index] += dxi[flux_dirn]*(Stildez_flux[index] - Stildez_flux[indp1]);

//        D_lapse[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.lapse - metric_face.lapse);
//        D_betax[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.betax - metric_face.betax);
//        D_betay[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.betay - metric_face.betay);
//        D_betaz[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.betaz - metric_face.betaz);
//        D_gxx[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.gxx - metric_face.gxx);
//        D_gxy[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.gxy - metric_face.gxy);
//        D_gxz[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.gxz - metric_face.gxz);
//        D_gyy[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.gyy - metric_face.gyy);
//        D_gyz[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.gyz - metric_face.gyz);
//        D_gzz[flux_dirn][index] = dxi[flux_dirn]*(metric_facep1.gzz - metric_face.gzz);
  }

//#pragma omp parallel for
//  for(int k=ghostzone; k<dirlength-ghostzone; k++)
//    for(int j=ghostzone; j<dirlength-ghostzone; j++)
//      for(int i=ghostzone; i<dirlength-ghostzone; i++) {
//        const int index  = CCTK_GFINDEX3D(cctkGH, i, j ,k);
//
//        metric_quantities metric;
//        initialize_metric(lapse[index],
//                          gxx[index], gxy[index], gxz[index],
//                          gyy[index], gyz[index], gzz[index],
//                          betax[index], betay[index], betaz[index],
//                          &metric);
//
//        extrinsic_curvature curv;
//        initialize_extrinsic_curvature(
//                          kxx[index], kxy[index], kxz[index],
//                          kyy[index], kyz[index], kzz[index],
//                          &curv);
//
//        metric_derivatives metric_derivs;
//        for(int dirn=0; dirn<3; dirn++) {
//          metric_derivs.lapse[dirn] = D_lapse[dirn][index];
//          metric_derivs.betax[dirn] = D_betax[dirn][index];
//          metric_derivs.betay[dirn] = D_betay[dirn][index];
//          metric_derivs.betaz[dirn] = D_betaz[dirn][index];
//          metric_derivs.adm_gxx[dirn] = D_gxx[dirn][index];
//          metric_derivs.adm_gxy[dirn] = D_gxy[dirn][index];
//          metric_derivs.adm_gxz[dirn] = D_gxz[dirn][index];
//          metric_derivs.adm_gyy[dirn] = D_gyy[dirn][index];
//          metric_derivs.adm_gyz[dirn] = D_gyz[dirn][index];
//          metric_derivs.adm_gzz[dirn] = D_gzz[dirn][index];
//        }
//
//        primitive_quantities prims;
//        initialize_primitives(rho[index], press[index], poison,
//                              vx[index], vy[index], vz[index],
//                              Bx[index], By[index], Bz[index],
//                              poison, poison, poison, // entropy, Y_e, temp
//                              &prims);
//        prims.u0  = rho[index]*Bx[index] / vy[index];
//
//
//        conservative_quantities cons_sources;
//        calculate_all_source_terms(&prims,
//                                   &eos,
//                                   &metric,
//                                   &curv,
//                                   &metric_derivs,
//                                   &cons_sources);
//                                    
//        tau_rhs[index] += cons_sources.tau;
//        Stildex_rhs[index] += cons_sources.S_x;
//        Stildey_rhs[index] += cons_sources.S_y;
//        Stildez_rhs[index] += cons_sources.S_z;
//  }
}
