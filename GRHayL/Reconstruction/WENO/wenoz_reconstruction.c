#include "ghl_reconstruction.h"

/**
 * @ingroup weno
 * @brief Reconstructs variables at the points  
 * @sp10 \f$ Ur(i) = U \left(i-\frac{1}{2} + \epsilon \right) \f$  
 * @sp10 \f$ Ul(i) = U \left(i-\frac{1}{2} - \epsilon \right) \f$
 *
 * @details
 * This function computes the right and left values of the left face of variable
 * \f$ U \f$ using the WENO-Z method. For example, to reconstruct the values
 * for the face \f$ i-\frac{1}{2} \f$, the array \f$ U \f$ should contain a
 * stencil centered around that face (i.e. values from \f$ i-2 \f$ to \f$ i+1 \f$).
 *
 * The method itself is described in the low-level function
 * @ref ghl_wenoz_reconstruction_for_cell.
 * 
 * @todo The implementer (Terrence Pierre Jacques) should review and approve
 *       this documentation.
 *
 * @param[in] U:   1D array containing values of variable \f$ U \f$
 *
 * @param[out] Ur: pointer to a double; set to the value of the right side of the face
 *
 * @param[out] Ul: pointer to a double; set to the value of the left side of the face
 *
 * @returns void
 */
void ghl_wenoz_reconstruction(
      const double U[6],
      double *restrict Ur,
      double *restrict Ul) {

  double tmpr, tmpl;        

  ghl_wenoz_reconstruction_for_cell(&U[0], &tmpr, &tmpl);
  *Ul = tmpr;

  ghl_wenoz_reconstruction_for_cell(&U[1], &tmpr, &tmpl);
  *Ur = tmpl;
 }
