#include "ghl_reconstruction.h"

/**
 * @ingroup ppm
 * @brief Reconstructs variables at the points  
 * @sp10 \f$ Ur(i) = U \left(i-\frac{1}{2} + \epsilon \right) \f$  
 * @sp10 \f$ Ul(i) = U \left(i-\frac{1}{2} - \epsilon \right) \f$
 *
 * @details
 * This function computes face values of a variable using the piecewise
 * parabolic method (PPM). The returned face depends on the provided stencil.
 * For example, providing a data stencil from \f$ (i-3) \f$ to \f$ (i+2) \f$ of the cell
 * centers will return the right and left values for the face at \f$ i-\frac{1}{2} \f$.
 * This involves two calls to @ref ghl_ppm_compute_for_cell, and this function filters
 * the returned data to provide the desired return values. For more detail on
 * the implemented PPM routine, see documentation on the @ref ghl_ppm_compute_for_cell
 * function.
 *
 * The input `ftilde` can be computed by using @ref ghl_compute_ftilde, which uses
 * the same stencil range as this function and populates the `ftilde` array with
 * the needed data.
 *
 * @param[in] ftilde:    1D array containing \f$ \tilde{f} \f$ values for flattening procedure
 *
 * @param[in] var_data:  1D array containing stencil of variable to reconstruct
 *
 * @param[out] var_datar: pointer to a double; set to the value of the right side of the face
 *
 * @param[out] var_datal: pointer to a double; set to the value of the left side of the face
 *
 * @returns void
 */
void ghl_ppm_reconstruction(
      const double ftilde[2],
      const double var_data[6],
      double *restrict var_datar,
      double *restrict var_datal) {

  double tmpr, tmpl;

  ghl_ppm_compute_for_cell(ftilde[0], var_data, &tmpr, &tmpl);
  *var_datal = tmpr;

  ghl_ppm_compute_for_cell(ftilde[1], &var_data[1], &tmpr, &tmpl);
  *var_datar = tmpl;
}
