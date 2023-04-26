#include "con2prim.h"

double compute_Bsq_from_Bup(
      const metric_quantities *restrict metric,
      const double *restrict Bup ) {

  double Bsq = 0.0;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      Bsq += metric->g4dn[i+1][j+1] * Bup[i] * Bup[j];

  return Bsq;
}
