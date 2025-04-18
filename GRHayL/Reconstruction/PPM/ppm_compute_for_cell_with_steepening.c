#include "ghl_reconstruction.h"

/**
 * @ingroup ppm
 * @brief Reconstructs variables at the points  
 * @sp10 \f$ Ur(i) = U \left(i-\frac{1}{2} + \epsilon \right) \f$  
 * @sp10 \f$ Ul(i) = U \left(i+\frac{1}{2} - \epsilon \right) \f$
 *
 * @details
 * This function computes the right and left face values of a cell for a
 * single variable using the piecewise parabolic method (PPM). As this
 * function is nearly identical to @ref ghl_ppm_compute_for_cell, we do
 * not repeat the full discussion of the PPM method here. Instead, we
 * simply discuss the key difference between the functions. This function
 * applies an additional steepening algorithm to the given variable. This
 * is normally used for the density, as the density profile should be
 * narrowed in the presence of the contact discontinuity (see e.g. Appendix I
 * of \cite Marti_1996.
 * This procedure is implemented via the @ref ghl_steepen_var function.
 *
 * @param[in] params:   pointer to ghl_parameters struct
 *
 * @param[in] pressure: 1D array containing stencil for the pressure
 *
 * @param[in] Gamma_eff: value of \f$ \Gamma \f$ to be used during the steepening procedure
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
void ghl_ppm_compute_for_cell_with_steepening(
      const ghl_parameters *restrict params,
      const double pressure[5],
      const double Gamma_eff,
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

  // Apply steepening algorithm
  ghl_steepen_var(params, U, pressure, Gamma_eff, &Ur, &Ul);

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
