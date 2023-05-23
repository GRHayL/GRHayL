#include "GRHayLHDX.h"

extern "C" void GRHayLHDX_compute_metric_derivs(
      const CCTK_REAL dxi,
      const Loop::GF3D2index indm2,
      const Loop::GF3D2index indm1,
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
      metric_quantities *restrict metric_derivs) {

  const CCTK_REAL d_lapse = dxi*COMPUTE_DERIV(lapse(indm2), lapse(indm1), lapse(indp1), lapse(indp2));
  const CCTK_REAL d_betax = dxi*COMPUTE_DERIV(betax(indm2), betax(indm1), betax(indp1), betax(indp2));
  const CCTK_REAL d_betay = dxi*COMPUTE_DERIV(betay(indm2), betay(indm1), betay(indp1), betay(indp2));
  const CCTK_REAL d_betaz = dxi*COMPUTE_DERIV(betaz(indm2), betaz(indm1), betaz(indp1), betaz(indp2));

  const CCTK_REAL d_gxx = dxi*COMPUTE_DERIV(gxx(indm2), gxx(indm1), gxx(indp1), gxx(indp2));
  const CCTK_REAL d_gxy = dxi*COMPUTE_DERIV(gxy(indm2), gxy(indm1), gxy(indp1), gxy(indp2));
  const CCTK_REAL d_gxz = dxi*COMPUTE_DERIV(gxz(indm2), gxz(indm1), gxz(indp1), gxz(indp2));
  const CCTK_REAL d_gyy = dxi*COMPUTE_DERIV(gyy(indm2), gyy(indm1), gyy(indp1), gyy(indp2));
  const CCTK_REAL d_gyz = dxi*COMPUTE_DERIV(gyz(indm2), gyz(indm1), gyz(indp1), gyz(indp2));
  const CCTK_REAL d_gzz = dxi*COMPUTE_DERIV(gzz(indm2), gzz(indm1), gzz(indp1), gzz(indp2));

  ghl_initialize_metric(
        d_lapse,
        d_betax, d_betay, d_betaz,
        d_gxx, d_gxy, d_gxz,
        d_gyy, d_gyz, d_gzz,
        metric_derivs);
}
