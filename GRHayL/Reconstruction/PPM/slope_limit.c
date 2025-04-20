#include "ghl_reconstruction.h"

/**
 * @ingroup recon_internal
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
 * @ingroup recon_internal
 * @brief Applies slope limiter to a reconstructed variable.
 *
 * @details
 * This function computes the derivative at a point while
 * limiting the slope by comparing the first- and second-order
 * derivatives. The core algorithm is quite simple, though we
 * expand it here for more detailed explanations. This method
 * comes from equation 60 of \cite Marti_1996. Note, however,
 * the factor of 2 missing in the
 * \f$ \left| a_{j+1} - a_{j} \right| \f$ term.
 *
 * We first compute the first order derivative
 *
 * \f[
 * \begin{aligned}
 * \Delta U &= \frac{dU + dUp1}{2} \\
 *          &= \frac{U_{i} - U_{i-1} + U_{i+1} - U_{i}}{2} \\
 *          &= \frac{U_{i+1} - U_{i-1}}{2}
 * \end{aligned}
 * \f]
 *
 * Then, we determine the sign and magnitude of the slope using
 *
 * \f[
 * \begin{align}
 * sgn(\mathrm{slope}) &= (0 < \Delta U) - (\Delta U < 0) \\
 * &=
 * \begin{cases}
 *   +1 & \mathrm{if\ } \Delta U > 0 \\
 *   -1 & \mathrm{if\ } \Delta U < 0 \\
 *    0 & \mathrm{otherwise}
 * \end{cases}
 * \end{align}
 * \f]
 *
 * \f[
 * \left|\mathrm{slope}\right| = \mathrm{min}\Bigl[ |\Delta U| ,\ 2 \mathrm{min}\left( |dU| ,\ |dUp1| \right) \Bigr]
 * \f]
 *
 * where the slope magnitude is limited by the second equation.
 *
 * @param[in] dU: value of the difference \f$ U_{i} - U_{i-1} \f$
 *
 * @param[in] dUp1: value of the difference \f$ U_{i+1} - U_{i} \f$
 *
 * @returns value of the limited slope at point \f$ i \f$
 */
double ghl_slope_limit(
      const double dU,
      const double dUp1) {

  if(dU*dUp1 > 0.0) {
    const double delta_m_U = 0.5*(dU + dUp1);
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
