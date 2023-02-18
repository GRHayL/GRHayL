#include "stdlib.h"
#include "unit_tests.h"

void randomize_metric(
      double *restrict lapse,
      double *restrict gxx_ptr,
      double *restrict gxy_ptr,
      double *restrict gxz_ptr,
      double *restrict gyy_ptr,
      double *restrict gyz_ptr,
      double *restrict gzz_ptr,
      double *restrict betax,
      double *restrict betay,
      double *restrict betaz) {
  const double gyy = 1.0 + randf(0.0,1.0e-1);
  const double gzz = 1.0 + randf(0.0,1.0e-1);
  const double gxy = randf(-1.0e-1,1.0e-1);
  const double gxz = randf(-1.0e-1,1.0e-1);
  const double gyz = randf(-1.0e-1,1.0e-1);
  const double phi = randf(0.0,2.0);
  const double gxx = ( 12.0*phi - gxy * (gyz*gxz - gxy*gzz) - gxz * (gxy*gyz - gyy*gxz) )/(gyy*gzz - gyz*gyz);

  *gxx_ptr = gxx;
  *gxy_ptr = gxy;
  *gxz_ptr = gxz;
  *gyy_ptr = gyy;
  *gyz_ptr = gyz;
  *gzz_ptr = gzz;
  *lapse   = 1.0;
  *betax   = 0.0;
  *betay   = 0.0;
  *betaz   = 0.0;
}
