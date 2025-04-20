#include "ghl_reconstruction.h"

/**
 * @ingroup ppm_internal
 * @brief Computes flattening parameter needed by the PPM algorithm.
 *
 * @details
 * This function computes the \f$ \tilde{f} \f$ value needed for the
 * flattening procedure inside of the PPM reconstruction routine.
 * It aims to reduce oscillation near shocks where a higher value
 * corresponds to stronger flattening of the interpolation profile.
 *
 * This procedure is described in the appendix of \cite Colella_1984
 * from which it originates, but we provide an overview of the
 * implementation here.
 *
 * To begin, we compute differences and averages of the pressure across the cell:
 *
 * \f[
 * \begin{aligned}
 * \Delta P_1 &\equiv P_{i+1} - P_{i-1} \\
 * \Delta P_2 &\equiv P_{i+2} - P_{i-2} \\
 * A_1 &\equiv \frac{P_{i+1} + P_{i-1}}{2} \\
 * A_2 &\equiv \frac{P_{i+2} + P_{i-2}}{2}
 * \end{aligned}
 * \f]
 *
 * At this point, we introduce protection against numerical round-off issues by
 * setting the difference to zero if it is less than \f$ 10^{-15} \f$.
 *
 * We also define the ratio of the differences
 *
 * \f[
 * R_\Delta \equiv \frac{\Delta P_1}{\Delta P_2}
 * \f]
 *
 * This is set to 1 if \f$ \Delta P_2=0 \f$. Finally, we use these quantities
 * to determine if flattening needs to occur with the triggers
 *
 * \f[
 * q_2 > \epsilon
 * \f]
 *
 * \f[
 * q_2 \left( v_{i-1} - v_{i+1} \right) > 0
 * \f]
 *
 * where
 *
 * \f[
 * q_2 = \frac{\left| \Delta P_1 \right|}{\min{\left(P_{i+1}, P_{i-1}\right)}}
 * \f]
 *
 * If both of these conditions are met, then there is a shock. In this case,
 * the value of \f$ \tilde{f} \f$ is given by
 *
 * \f[
 * \tilde{f} = \min{\left[ 1, \max{\left( 0, q_1 \right)} \right]}
 * \f]
 *
 * where
 *
 * \f[
 * q_1 = \omega_2 \left(R_\Delta - \omega_1 \right)
 * \f]
 *
 * If either condition is false, then there is no shock and thus flattening
 * is unnecessary. Then, the function simply returns zero for \f$ \tilde{f} \f$.
 *
 * @param[in] params:      pointer to ghl_parameters struct
 *
 * @param[in] P:           1D array containing stencil for the pressure
 *
 * @param[in] v_flux_dirn: 1D array containing stencil for the fluid
 *                         velocity in the direction of the reconstruction
 *
 * @returns the computed \f$ \tilde{f} \f$ value
 */
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
  const double zero_cutoff = 1e-15;

  if(fabs(dP1)/avg1 < zero_cutoff) {
    // If this is triggered, there is no shock
    return 0.0;
  }
  if(fabs(dP2)/avg2 < zero_cutoff) {
    // If this is triggered, there may still be a shock
    dP2 = 0.0;
  }

  double dP1_over_dP2 = 1.0;
  if (dP2 != 0.0) {
    dP1_over_dP2 = dP1/dP2;
  }

  const double q1 = (dP1_over_dP2 - params->ppm_flattening_omega1)
                    * params->ppm_flattening_omega2;
  const double q2 = fabs(dP1)/MIN(P[PLUS_1], P[MINUS1]);

  // this if statement is equivalent to the w_j variable in the original Colella and Woodward paper
  if (q2 > params->ppm_flattening_epsilon && q2*( (v_flux_dirn[MINUS1]) - (v_flux_dirn[PLUS_1]) ) > 0.0) {
    // inside a shock
    return ghl_clamp(q1, 0.0, 1.0);
  } else {
    // Not inside a shock
    return 0.0;
  }
}
