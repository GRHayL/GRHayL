#include "ghl_reconstruction.h"

// Compute ftilde, which is used for flattening left and right face values
// DEPENDENCIES: P(MINUS2,MINUS1,PLUS_1,PLUS_2) and v^m(MINUS1,PLUS_1), where m=flux_dirn={1,2,3}={x,y,z}.

double ghl_shock_detection_ftilde(
      const ghl_parameters *restrict params,
      const double P[5],
      const double v_flux_dirn[5]) {

  double dP1 = P[PLUS_1] - P[MINUS1];
  double dP2 = P[PLUS_2] - P[MINUS2];

  const double avg1 = 0.5*(P[PLUS_1] + P[MINUS1]);
  const double avg2 = 0.5*(P[PLUS_2] + P[MINUS2]);

  // MODIFICATION TO STANDARD PPM:
  // Cure roundoff error issues when dP1==0 or dP2==0 to 15 or more significant digits.
  if(fabs(dP1)/avg1 < 1e-15) {
    /* If this is triggered, there is NO shock */
    dP1 = 0.0;
  }
  if(fabs(dP2)/avg2 < 1e-15) {
    /* If this is triggered, there may still be a shock */
    dP2 = 0.0;
  }

  double dP1_over_dP2 = 1.0;
  if (dP2 != 0.0) dP1_over_dP2 = dP1/dP2;

  const double q1 = (dP1_over_dP2 - params->ppm_flattening_omega1) * params->ppm_flattening_omega2;
  const double q2 = fabs(dP1)/MIN(P[PLUS_1], P[MINUS1]);

  // this if statement is equivalent to the w_j variable in the original Colella and Woodward paper
  if (q2 > params->ppm_flattening_epsilon && q2*( (v_flux_dirn[MINUS1]) - (v_flux_dirn[PLUS_1]) ) > 0.0) {
    // inside a shock
    return MIN(1.0, MAX(0.0, q1));
  } else {
    // NOT inside a shock
    return 0.0;
  }
}
