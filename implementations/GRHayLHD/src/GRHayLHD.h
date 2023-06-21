#ifndef GRHAYLHD_H_
#define GRHAYLHD_H_

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

// Computes derivative factor dx*deriv
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) ((A0 - AM1)*(Varp1 - Varm1) + AM1*(Varp2 - Varm2))

void GRHayLHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
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
      metric_quantities *restrict metric);

void GRHayLHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
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
      metric_quantities *restrict metric_derivs);

#endif // GRHAYLHD_H_
