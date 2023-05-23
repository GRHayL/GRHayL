#ifndef GRHAYLHDX_H_
#define GRHAYLHDX_H_

#include "loop_device.hxx"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

// Interpolates to the +1/2 face
#define AM1 -0.0625
#define A0   0.5625
#define A1   0.5625
#define A2  -0.0625
#define COMPUTE_FCVAL(Varm1,Var,Varp1,Varp2) (AM1*(Varm1) + A0*(Var) + A1*(Varp1) + A2*(Varp2))

/*
   Computes derivative factor by computing faceval[i] - faceval[i-1]:
   Let A = AM1 = A2, B = A0 = A1. Let Var at index i be f[i]. Then,
   dx*deriv = Af[-1] + Bf[0] + Bf[1] + Af[2] - (Af[-2] + Bf[-1] + Bf[0] + Af[1])
            = Af[-1] - Bf[-1] + Bf[1] - Af[1] + Af[2] - Af[-2]
            = (B-A)(f[1] - f[-1]) + A(f[2] - f[-2])
*/ 
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) ((A0 - AM1)*(Varp1 - Varm1) + AM1*(Varp2 - Varm2))

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
      metric_quantities *restrict metric);

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
      metric_quantities *restrict metric_derivs);

#endif // GRHAYLHDX_H_
