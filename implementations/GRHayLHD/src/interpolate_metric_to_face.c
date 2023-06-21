#include "GRHayLHD.h"

void GRHayLHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
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
      metric_quantities *restrict metric) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
  const int index  = CCTK_GFINDEX3D(cctkGH, i,        j ,       k);
  const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);
  const int indp2  = CCTK_GFINDEX3D(cctkGH, i+2*xdir, j+2*ydir, k+2*zdir);

  const CCTK_REAL face_lapse = COMPUTE_FCVAL(lapse[indm1], lapse[index], lapse[indp1], lapse[indp2]);
  const CCTK_REAL face_betax = COMPUTE_FCVAL(betax[indm1], betax[index], betax[indp1], betax[indp2]);
  const CCTK_REAL face_betay = COMPUTE_FCVAL(betay[indm1], betay[index], betay[indp1], betay[indp2]);
  const CCTK_REAL face_betaz = COMPUTE_FCVAL(betaz[indm1], betaz[index], betaz[indp1], betaz[indp2]);

  const CCTK_REAL face_gxx = COMPUTE_FCVAL(gxx[indm1], gxx[index], gxx[indp1], gxx[indp2]);
  const CCTK_REAL face_gxy = COMPUTE_FCVAL(gxy[indm1], gxy[index], gxy[indp1], gxy[indp2]);
  const CCTK_REAL face_gxz = COMPUTE_FCVAL(gxz[indm1], gxz[index], gxz[indp1], gxz[indp2]);
  const CCTK_REAL face_gyy = COMPUTE_FCVAL(gyy[indm1], gyy[index], gyy[indp1], gyy[indp2]);
  const CCTK_REAL face_gyz = COMPUTE_FCVAL(gyz[indm1], gyz[index], gyz[indp1], gyz[indp2]);
  const CCTK_REAL face_gzz = COMPUTE_FCVAL(gzz[indm1], gzz[index], gzz[indp1], gzz[indp2]);

  ghl_initialize_metric(
        face_lapse,
        face_betax, face_betay, face_betaz,
        face_gxx, face_gxy, face_gxz,
        face_gyy, face_gyz, face_gzz,
        metric);
}
