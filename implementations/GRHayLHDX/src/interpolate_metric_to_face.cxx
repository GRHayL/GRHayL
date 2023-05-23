#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_interpolate_metric_to_face(
      const Loop::GF3D2index indm1,
      const Loop::GF3D2index index,
      const Loop::GF3D2index indp1,
      const Loop::GF3D2index indp2,
      const Loop::GF3D2<const CCTK_REAL> lapse,
      const Loop::GF3D2<const CCTK_REAL> betax,
      const Loop::GF3D2<const CCTK_REAL> betay,
      const Loop::GF3D2<const CCTK_REAL> betaz,
      const Loop::GF3D2<const CCTK_REAL> gxx,
      const Loop::GF3D2<const CCTK_REAL> gxy,
      const Loop::GF3D2<const CCTK_REAL> gxz,
      const Loop::GF3D2<const CCTK_REAL> gyy,
      const Loop::GF3D2<const CCTK_REAL> gyz,
      const Loop::GF3D2<const CCTK_REAL> gzz,
      metric_quantities *restrict metric) {

  const CCTK_REAL face_lapse = COMPUTE_FCVAL(lapse(indm1), lapse(index), lapse(indp1), lapse(indp2));
  const CCTK_REAL face_betax = COMPUTE_FCVAL(betax(indm1), betax(index), betax(indp1), betax(indp2));
  const CCTK_REAL face_betay = COMPUTE_FCVAL(betay(indm1), betay(index), betay(indp1), betay(indp2));
  const CCTK_REAL face_betaz = COMPUTE_FCVAL(betaz(indm1), betaz(index), betaz(indp1), betaz(indp2));

  const CCTK_REAL face_gxx = COMPUTE_FCVAL(gxx(indm1), gxx(index), gxx(indp1), gxx(indp2));
  const CCTK_REAL face_gxy = COMPUTE_FCVAL(gxy(indm1), gxy(index), gxy(indp1), gxy(indp2));
  const CCTK_REAL face_gxz = COMPUTE_FCVAL(gxz(indm1), gxz(index), gxz(indp1), gxz(indp2));
  const CCTK_REAL face_gyy = COMPUTE_FCVAL(gyy(indm1), gyy(index), gyy(indp1), gyy(indp2));
  const CCTK_REAL face_gyz = COMPUTE_FCVAL(gyz(indm1), gyz(index), gyz(indp1), gyz(indp2));
  const CCTK_REAL face_gzz = COMPUTE_FCVAL(gzz(indm1), gzz(index), gzz(indp1), gzz(indp2));

  ghl_initialize_metric(
        face_lapse,
        face_betax, face_betay, face_betaz,
        face_gxx, face_gxy, face_gxz,
        face_gyy, face_gyz, face_gzz,
        metric);
}
