#ifndef GRHAYL_IGH_H_
#define GRHAYL_IGH_H_

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLET.h"

//Interpolates to the +1/2 face
#define AM1 -0.0625
#define A0   0.5625
#define A1   0.5625
#define A2  -0.0625
#define COMPUTE_FCVAL(Varm1,Var,Varp1,Varp2) (AM1*(Varm1) + A0*(Var) + A1*(Varp1) + A2*(Varp2))

void GRHayL_IGH_interpolate_to_face_and_initialize_metric(
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
      metric_quantities *restrict metric);

#endif // GRHAYL_IGH_H_
