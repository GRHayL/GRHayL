#include "ghl_reconstruction.h"

/**
 * @ingroup ppm
 * @brief Computes flattening parameters needed by the PPM algorithm.
 *
 * @details
 * This function computes the \f$ \tilde{f} \f$ values needed for the
 * flattening procedure inside of the PPM reconstruction routine.
 * As with all the reconstruction functions, this function is designed
 * to facilitate the computation of the right and left reconstructed
 * variables for a face instead of a cell. As such, this function
 * computes two values for \f$ \tilde{f} \f$, one for each side of the face.
 * This 2-element array can then be directly passed to the PPM routine.
 *
 * This routine uses the @grhayl parameters ghl_parameters::ppm_flattening_epsilon,
 * ghl_parameters::ppm_flattening_omega1, and ghl_parameters::ppm_flattening_omega2.
 * The @ref ghl_initialize_params function sets these to default values matching
 * \cite Colella_1984, but they can be changed at any time.
 *
 * The core @ref ghl_shock_detection_ftilde routine for a single value uses a
 * 5-element stencil of the input values to compute \f$ \tilde{f} \f$, and the
 * method for computing \f$ \tilde{f} \f$ is discussed in detail in its documentation.
 *
 * @param[in] params:      pointer to ghl_parameters struct
 *
 * @param[in] pressure:    1D array containing stencil for the pressure
 *
 * @param[in] v_flux_dirn: 1D array containing stencil for the fluid
 *                         velocity in the direction of the reconstruction
 *
 * @param[out] ftilde:     1D array containing \f$ \tilde{f} \f$ values
 *
 * @returns void
 */
void ghl_compute_ftilde(
      const ghl_parameters *restrict params,
      const double pressure[6],
      const double v_flux_dirn[6],
      double ftilde[2]) {

  ftilde[0] = ghl_shock_detection_ftilde(params, pressure, v_flux_dirn);
  ftilde[1] = ghl_shock_detection_ftilde(params, &pressure[1], &v_flux_dirn[1]);
}
