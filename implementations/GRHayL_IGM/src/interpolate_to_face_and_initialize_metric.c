#include "cctk.h"
#include "IGM.h"

#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

void interpolate_to_face_and_initialize_metric(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const double *restrict lapse,
      const double *restrict betax,
      const double *restrict betay,
      const double *restrict betaz,
      const double *restrict gxx,
      const double *restrict gxy,
      const double *restrict gxz,
      const double *restrict gyy,
      const double *restrict gyz,
      const double *restrict gzz,
      metric_quantities *restrict metric) {

  const int xdir = (flux_dirn == 0);
  const int ydir = (flux_dirn == 1);
  const int zdir = (flux_dirn == 2);

  const int indm2  = CCTK_GFINDEX3D(cctkGH, i-2*xdir, j-2*ydir, k-2*zdir);
  const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
  const int index  = CCTK_GFINDEX3D(cctkGH, i,        j ,       k);
  const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);

  const double face_lapse = COMPUTE_FCVAL(lapse[indm2], lapse[indm1], lapse[index], lapse[indp1]);
  const double face_betax = COMPUTE_FCVAL(betax[indm2], betax[indm1], betax[index], betax[indp1]);
  const double face_betay = COMPUTE_FCVAL(betay[indm2], betay[indm1], betay[index], betay[indp1]);
  const double face_betaz = COMPUTE_FCVAL(betaz[indm2], betaz[indm1], betaz[index], betaz[indp1]);

  const double face_gxx = COMPUTE_FCVAL(gxx[indm2], gxx[indm1], gxx[index], gxx[indp1]);
  const double face_gxy = COMPUTE_FCVAL(gxy[indm2], gxy[indm1], gxy[index], gxy[indp1]);
  const double face_gxz = COMPUTE_FCVAL(gxz[indm2], gxz[indm1], gxz[index], gxz[indp1]);
  const double face_gyy = COMPUTE_FCVAL(gyy[indm2], gyy[indm1], gyy[index], gyy[indp1]);
  const double face_gyz = COMPUTE_FCVAL(gyz[indm2], gyz[indm1], gyz[index], gyz[indp1]);
  const double face_gzz = COMPUTE_FCVAL(gzz[indm2], gzz[indm1], gzz[index], gzz[indp1]);

  initialize_metric(face_lapse,
                    face_gxx, face_gxy, face_gxz,
                    face_gyy, face_gyz, face_gzz,
                    face_betax, face_betay, face_betaz,
                    metric);
}
