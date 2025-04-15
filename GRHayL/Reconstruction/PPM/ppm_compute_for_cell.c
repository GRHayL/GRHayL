#include "ghl_reconstruction.h"

/**
 * @ingroup ppm
 * @brief Reconstructs variables at the points  
 * @sp10 \f$ Ur(i) = U \left(i-\frac{1}{2} + \epsilon \right) \f$  
 * @sp10 \f$ Ul(i) = U \left(i+\frac{1}{2} - \epsilon \right) \f$
 *
 * @details
 * This function computes the right and left face values of a cell for a
 * single variable using the piecewise parabolic method (PPM). This
 * process occurs in several stages. Before that, it is important to
 * note that this method assumes a uniform grid, and thus all the 'slopes'
 * here do not divide by the grid spacing because it cancels out in the end.
 * Hence, they are technically only differences between grid functions instead
 * of actual approximate slopes.
 *
 * First, we use the slope limiter to compute the limited slopes for various
 * stencil combinations. The @ref ghl_slope_limit function handles the individual
 * slope limits, and we compute the following arguments:
 *
 * \f[
 * \begin{align}
 * s_- &= \mathrm{slope\_limit} \left( U_{i-1} - U_{i-2}, U_i - U_{i-1} \right) \\
 * s_0 &= \mathrm{slope\_limit} \left( U_i - U_{i-1},     U_{i+1} - U_i \right) \\
 * s_+ &= \mathrm{slope\_limit} \left( U_{i+1} - U_i,     U_{i+2} - U_{i+1} \right)
 * \end{align}
 * \f]
 *
 * These slopes are then used to compute the face values for the cell using
 *
 * \f[
 * \begin{align}
 * U_r &= \frac{U_{i-1} + U_i}{2} + \frac{s_0 - s_+}{6} \\
 * U_l &= \frac{U_i + U_{i+1}}{2} + \frac{s_- - s_0}{6}
 * \end{align}
 * \f]
 *
 * This initial value is then flattened using the provided \f$ \tilde{f} \f$ values
 *
 * \f[
 * \begin{align}
 * U_r &= \tilde{f} U_i + U_r \left( 1 - \tilde{f} \right) \\
 * U_l &= \tilde{f} U_i + U_l \left( 1 - \tilde{f} \right)
 * \end{align}
 * \f]
 *
 * After flattening, the values are then monotonized. There are several cases here.
 *
 * If the face values and cell value are not monotonic (i.e.
 * \f$ (U_r - U_i) (U_i - U_l) \leq 0 \f$), then both face values
 * are set to the cell value
 *
 * \f[
 * U_r = U_l = U_i
 * \f]
 *
 * Otherwise, we compute the slope from the faces
 *
 * \f[
 * \Delta U = U_r - U_l
 * \f]
 *
 * and also multiply this by the slope of the cell value and the average of the faces
 *
 * \f[
 * U_\mathrm{temp} = \Delta U \left( U_i - \frac{U_r + U_l}{2} \right)
 * \f]
 *
 * Then, we check if \f$ U_\mathrm{temp} \f$ is too large with respect to
 * \f$ \Delta U \f$. The values default to
 *
 * \f[
 * \begin{align}
 * U_r &= U_r \\
 * U_l &= U_l
 * \end{align}
 * \f]
 *
 * However, if
 *
 * \f[
 * U_\mathrm{temp} > \frac{\left(\Delta U\right)^2}{6}
 * \f]
 *
 * then the left value is changed to
 *
 * \f[
 * U_l = 3U_i - 2U_r
 * \f]
 *
 * Similarly, if
 *
 * \f[
 * U_\mathrm{temp} < -\frac{\left(\Delta U\right)^2}{6}
 * \f]
 *
 * then the right value is changed to
 *
 * \f[
 * U_r = 3U_i - 2U_l
 * \f]
 *
 * @param[in] ftilde:  \f$ \tilde{f} \f$ for the flattening procedure
 *
 * @param[in] U:       1D array containing stencil of variable to reconstruct
 *
 * @param[out] Ur_ptr: pointer to a double; set to the value of the right face
 *
 * @param[out] Ul_ptr: pointer to a double; set to the value of the left face
 *
 * @returns void
 */
void ghl_ppm_compute_for_cell(
      const double ftilde,
      const double U[5],
      double *restrict Ur_ptr,
      double *restrict Ul_ptr) {

  const double U0 = U[PLUS_0];

  const double slope_limited_dU_m1 = ghl_slope_limit(U[MINUS1] - U[MINUS2], U0        - U[MINUS1]);
  const double slope_limited_dU_p0 = ghl_slope_limit(U0        - U[MINUS1], U[PLUS_1] - U0);
  const double slope_limited_dU_p1 = ghl_slope_limit(U[PLUS_1] - U0,        U[PLUS_2] - U[PLUS_1]);

  double Ur = 0.5*(U[PLUS_1] + U0) + (1.0/6.0)*(slope_limited_dU_p0 - slope_limited_dU_p1);
  double Ul = 0.5*(U0 + U[MINUS1]) + (1.0/6.0)*(slope_limited_dU_m1 - slope_limited_dU_p0);

  // First detect shocks / steep gradients
  // and flatten variables
  Ur = U0*ftilde + Ur*(1.0 - ftilde);
  Ul = U0*ftilde + Ul*(1.0 - ftilde);

  // Then monotonize all variables
  if ( (Ur - U0)*(U0 - Ul) <= 0.0) {
    *Ur_ptr = U0;
    *Ul_ptr = U0;
    return;
  }

  const double dU = Ur - Ul;
  const double Utmp = dU*( U0 - 0.5*(Ur + Ul) );

  if ( Utmp > (1.0/6.0)*(dU*dU)) {
    *Ur_ptr = Ur;
    *Ul_ptr = 3.0*U0 - 2.0*Ur;
  } else if ( Utmp < -(1.0/6.0)*(dU*dU)) {
    *Ur_ptr = 3.0*U0 - 2.0*Ul;
    *Ul_ptr = Ul;
  } else {
    *Ur_ptr = Ur;
    *Ul_ptr = Ul;
  }
}
