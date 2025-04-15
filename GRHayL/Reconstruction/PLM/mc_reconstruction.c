#include "ghl_reconstruction.h"

/**
 * @ingroup plm
 * @brief Reconstructs variables at the points  
 * @sp10 \f$ Ur(i) = U \left(i-\frac{1}{2} + \epsilon \right) \f$  
 * @sp10 \f$ Ul(i) = U \left(i-\frac{1}{2} - \epsilon \right) \f$
 *
 * @details
 * This function computes the right and left values of the left face of variable \f$ U \f$
 * using the monotized-central method. For example, to reconstruct the values
 * for the face \f$ i-\frac{1}{2} \f$, the array \f$ U \f$ should contain a
 * stencil centered around that face (i.e. values from \f$ i-2 \f$ to \f$ i+1 \f$).
 * For the right side value `Ur`,
 *
 * \f[
 * U_R \equiv U_{i+\epsilon-1/2} = U_i - \frac{\Delta U}{2}
 * \f]
 *
 * where \f$ \Delta U \f$ is determined by repeated @ref ghl_minmod function calls.
 * These calls effectively select the smallest magnitude \f$ \Delta U \f$ from
 * \f$ \frac{1}{2}(U_{i+1} - U_{i-1}) \f$, \f$ 2(U_{i} - U_{i-1}) \f$, and
 * \f$ 2(U_{i+1}-U_{i}) \f$.
 *
 * Similarly, the left side value `Ul` chooses the smallest magnitude \f$ \Delta U \f$
 * from \f$ \frac{1}{2}(U_{i} - U_{i-2}) \f$, \f$ 2(U_{i-1} - U_{i-2}) \f$, and
 * \f$ (U_{i}-U_{i-1}) \f$ to compute
 *
 * \f[
 * U_L \equiv U_{i-\epsilon-1/2} = U_{i-1} + \frac{\Delta U}{2}
 * \f]
 *
 * If the two options of \f$ \Delta U \f$ are of different sign, then the method sets
 * \f$ \Delta U=0 \f$, as with @ref ghl_minmod_reconstruction.
 *
 * @param[in] U:   1D array containing values of variable \f$ U \f$
 *
 * @param[out] Ur: pointer to a double; set to the value of the right side of the face
 *
 * @param[out] Ul: pointer to a double; set to the value of the left side of the face
 *
 * @returns void
 */
void ghl_mc_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul) {

  const double DeltaU_0_m1 = U[2] - U[1];

  const double tmp_sigma_i   = ghl_minmod(DeltaU_0_m1, U[3] - U[2]);
  const double tmp_sigma_im1 = ghl_minmod(U[1] - U[0], DeltaU_0_m1);

  const double sigma_i   = ghl_minmod(0.5*(U[3] - U[1]), 2*tmp_sigma_i);
  const double sigma_im1 = ghl_minmod(0.5*(U[2] - U[0]), 2*tmp_sigma_im1);

  *Ur = U[2] - 0.5*sigma_i;
  *Ul = U[1] + 0.5*sigma_im1;
 }
