#include "GRHayLMHD.h"

void GRHayLMHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
      const CCTK_REAL dxi,
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
      metric_quantities *restrict metric_derivs) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int indm2  = CCTK_GFINDEX3D(cctkGH, i-2*xdir, j-2*ydir, k-2*zdir);
  const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
  const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);
  const int indp2  = CCTK_GFINDEX3D(cctkGH, i+2*xdir, j+2*ydir, k+2*zdir);

  const CCTK_REAL d_lapse = dxi*COMPUTE_DERIV(lapse[indm2], lapse[indm1], lapse[indp1], lapse[indp2]);
  const CCTK_REAL d_betax = dxi*COMPUTE_DERIV(betax[indm2], betax[indm1], betax[indp1], betax[indp2]);
  const CCTK_REAL d_betay = dxi*COMPUTE_DERIV(betay[indm2], betay[indm1], betay[indp1], betay[indp2]);
  const CCTK_REAL d_betaz = dxi*COMPUTE_DERIV(betaz[indm2], betaz[indm1], betaz[indp1], betaz[indp2]);

  const CCTK_REAL d_gxx = dxi*COMPUTE_DERIV(gxx[indm2], gxx[indm1], gxx[indp1], gxx[indp2]);
  const CCTK_REAL d_gxy = dxi*COMPUTE_DERIV(gxy[indm2], gxy[indm1], gxy[indp1], gxy[indp2]);
  const CCTK_REAL d_gxz = dxi*COMPUTE_DERIV(gxz[indm2], gxz[indm1], gxz[indp1], gxz[indp2]);
  const CCTK_REAL d_gyy = dxi*COMPUTE_DERIV(gyy[indm2], gyy[indm1], gyy[indp1], gyy[indp2]);
  const CCTK_REAL d_gyz = dxi*COMPUTE_DERIV(gyz[indm2], gyz[indm1], gyz[indp1], gyz[indp2]);
  const CCTK_REAL d_gzz = dxi*COMPUTE_DERIV(gzz[indm2], gzz[indm1], gzz[indp1], gzz[indp2]);

  ghl_initialize_metric(
        d_lapse,
        d_betax, d_betay, d_betaz,
        d_gxx, d_gxy, d_gxz,
        d_gyy, d_gyz, d_gzz,
        metric_derivs);
}
