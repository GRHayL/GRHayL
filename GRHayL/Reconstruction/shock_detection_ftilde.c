#include "reconstruction.h"

// Compute ftilde, which is used for flattening left and right face values
// DEPENDENCIES: P(MINUS2,MINUS1,PLUS_1,PLUS_2) and v^m(MINUS1,PLUS_1), where m=flux_dirn={1,2,3}={x,y,z}.
#define OMEGA1   0.75
#define OMEGA2  10.0
#define EPSILON2 0.33

double grhayl_shock_detection_ftilde(
      const double P[5],
      const double v_flux_dirn[5]) {

  double dP1 = P[PLUS_1] - P[MINUS1];
  double dP2 = P[PLUS_2] - P[MINUS2];

  // MODIFICATION TO STANDARD PPM:
  // Cure roundoff error issues when dP1==0 or dP2==0 to 15 or more significant digits.
  const double avg1=0.5*(P[PLUS_1] + P[MINUS1]);
  const double avg2=0.5*(P[PLUS_2] + P[MINUS2]);
  if(fabs(dP1)/avg1<1e-15) dP1=0.0; /* If this is triggered, there is NO shock */
  if(fabs(dP2)/avg2<1e-15) dP2=0.0; /* If this is triggered alone, there may be a shock. Otherwise if triggered with above, NO shock. */

  double dP1_over_dP2=1.0;
  if (dP2 != 0.0) dP1_over_dP2 = dP1/dP2;

  const double q1 = (dP1_over_dP2-OMEGA1)*OMEGA2;
  const double q2 = fabs(dP1)/MIN(P[PLUS_1], P[MINUS1]);

  // w==0 -> NOT inside a shock
  // w==1 -> inside a shock
  const double w = (q2 > EPSILON2 && q2*( (v_flux_dirn[MINUS1]) - (v_flux_dirn[PLUS_1]) ) > 0.0);

  return MIN(1.0, w*MAX(0.0,q1));
}
