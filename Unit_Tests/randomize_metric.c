#include "stdlib.h"
#include "unit_tests.h"

void randomize_metric(metric_quantities *restrict metric) {
  double phi, gxx, gyy, gzz, gxy, gxz, gyz, lapse, betax, betay, betaz;
  gyy = 1.0 + randf(0.0,1.0e-1);
  gzz = 1.0 + randf(0.0,1.0e-1);
  gxy = randf(-1.0e-1,1.0e-1);
  gxz = randf(-1.0e-1,1.0e-1);
  gyz = randf(-1.0e-1,1.0e-1);
  phi = randf(0.0,2.0);
  gxx = ( 12.0*phi - gxy * (gyz*gxz - gxy*gzz) - gxz * (gxy*gyz - gyy*gxz) )/(gyy*gzz - gyz*gyz);
  lapse = 1.0;
  betax = 0.0;
  betay = 0.0;
  betaz = 0.0;

  initialize_metric(lapse,
                    gxx, gxy, gxz,
                    gyy, gyz, gzz,
                    betax, betay, betaz,
                    metric);
}
