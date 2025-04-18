#include "ghl_reconstruction.h"

/**
 * @ingroup ppm_internal
 * @brief Macro controlling the slope limit method. 2.0 = MC, 1 = minmod.
 *
 * @todo
 * I (S.Cupp) am not convinced that changing this to 1 will actually
 * recover a minmod method, though I do agree that 2 gives MC. One of the
 * following should probably be done:
 * - remove the (hard-coded) macro and just use the value in the code
 * - change the code to call the appropriate functions in PLM to make
 *   sure we know what it's doing
 */
#define SLOPE_LIMITER_COEFF 2.0

/**
 * @ingroup ppm_internal
 * @brief Applies slope limiter to a reconstructed variable.
 *
 * @details
 * This function computes the derivative at a point while
 * limiting the slope by comparing the first- and second-order
 * derivatives. The core algorithm is quite simple, though we
 * expand it here for more detailed explanations. This method
 * comes from equation 60 of \cite Marti_1996. We first compute
 * the first order derivative
 *
 * \f[
 * \begin{aligned}
 * \Delta U &= \frac{dU + dUp1}{2} \\
 *          &= \frac{U_{i} - U_{i-1} + U_{i+1} - U_{i}}{2} \\
 *          &= \frac{U_{i+1} - U_{i-1}}{2}
 * \end{aligned}
 * \f]
 *
 * Then, we determine the correct sign of the derivative using
 *
 * \f[
 * (0 < \Delta U) - (\Delta U < 0)
 * \f]
 *
 * This expression will give -1 if the slope is negative, +1 if
 * it is positive, and 0 if it is neither. Then, we compute the
 * magnitude of the slope with the expression
 *
 * \f[
 * \mathrm{min}\left[ |\Delta U| , 2 \mathrm{min}\left( |dU| , |dUp1| \right) \right]
 * \f]
 *
 * This operation takes the minimum of the 3 derivatives, thereby limiting the slope.
 *
 * @param[in] dU: value of the difference \f$ U_{i} - U_{i-1} \f$
 *
 * @param[in] dUp1: value of the difference \f$ U_{i+1} - U_{i} \f$
 *
 * @returns value of the limited slope at point \f$ i \f$
 */

// Eq. 60 \cite Marti_1996
// [note the factor of 2 missing in the |a_{j+1} - a_{j}| term].
// Recall that dU = U_{i} - U_{i-1}.
double ghl_slope_limit(
      const double dU,
      const double dUp1) {

  if(dU*dUp1 > 0.0) {
    //delta_m_U=0.5 * [ (u_(i+1)-u_i) + (u_i-u_(i-1)) ] = (u_(i+1) - u_(i-1))/2  <-- first derivative, second-order; this should happen most of the time (smooth flows)
    const double delta_m_U = 0.5*(dU + dUp1);
    // EXPLANATION OF BELOW LINE OF CODE.
    // In short, sign_delta_a_j = sign(delta_m_U) = (0.0 < delta_m_U) - (delta_m_U < 0.0).
    //    If delta_m_U>0, then (0.0 < delta_m_U)==1, and (delta_m_U < 0.0)==0, so sign_delta_a_j=+1
    //    If delta_m_U<0, then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==1, so sign_delta_a_j=-1
    //    If delta_m_U==0,then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==0, so sign_delta_a_j=0
    const int sign_delta_m_U = (0.0 < delta_m_U) - (delta_m_U < 0.0);
    //Decide whether to use 2nd order derivative or first-order derivative, limiting slope.
    return sign_delta_m_U*MIN(fabs(delta_m_U), SLOPE_LIMITER_COEFF*MIN(fabs(dUp1), fabs(dU)));
  }
  return 0.0;
}
