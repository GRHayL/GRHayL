#include "ghl_reconstruction.h"

/**
 * @ingroup plm
 * @brief Reconstructs variables at the points  
 * @sp10 \f$ Ur(i) = U \left(i-\frac{1}{2} + \epsilon \right) \f$  
 * @sp10 \f$ Ul(i) = U \left(i-\frac{1}{2} - \epsilon \right) \f$
 *
 * @details
 * This function computes the right and left values of the left face of variable \f$ U \f$
 * using the minmod method. For example, to reconstruct the values for the face
 * \f$ i-\frac{1}{2} \f$, the array \f$ U \f$ should contain a stencil centered
 * around that face (i.e. values from \f$ i-2 \f$ to \f$ i+1 \f$). This method
 * chooses the smaller \f$ \Delta U \f$ between \f$ U_i - U_{i-1} \f$ and
 * \f$ U_{i+1}-U \f$ to use to take a half-step backwards to the right side:
 *
 * \f[
 * U_R \equiv U_{i+\epsilon-1/2} = U_i - \frac{\Delta U}{2}
 * \f]
 *
 * Similarly, the smaller \f$ \Delta U \f$ between \f$ U_{i-1} - U_{i-2} \f$
 * and \f$ U_{i}-U_{i-1} \f$ is used to take a half-step forwards to the left side:
 *
 * \f[
 * U_L \equiv U_{i-\epsilon-1/2} = U_{i-1} + \frac{\Delta U}{2}
 * \f]
 *
 * If the two options of \f$ \Delta U \f$ are of different sign, then the method
 * sets \f$ \Delta U=0 \f$. The actual value selection is done via the @ref ghl_minmod
 * function.
 *
 * @param[in] U:   1D array containing values of variable \f$ U \f$
 *
 * @param[out] Ur: pointer to a double; set to the value of the right side of the face
 *
 * @param[out] Ul: pointer to a double; set to the value of the left side of the face
 *
 * @returns void
 */
void ghl_minmod_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul) {

  const double DeltaU_0_m1 = U[2]- U[1];

  const double sigma_i   = ghl_minmod(DeltaU_0_m1, U[3] - U[2]);
  const double sigma_im1 = ghl_minmod(U[1] - U[0], DeltaU_0_m1);

  *Ur = U[2] - 0.5*sigma_i;
  *Ul = U[1] + 0.5*sigma_im1;
 }
