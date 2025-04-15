#include "ghl_reconstruction.h"

/**
 * @ingroup plm
 * @brief Reconstructs variables at the points  
 * @sp10 \f$ Ur(i) = U \left(i-\frac{1}{2} + \epsilon \right) \f$  
 * @sp10 \f$ Ul(i) = U \left(i-\frac{1}{2} - \epsilon \right) \f$
 *
 * @details
 * This function computes the right and left values of the left face of variable \f$ U \f$
 * using the superbee method. For example, to reconstruct the values for the face
 * \f$ i-\frac{1}{2} \f$, the array \f$ U \f$ should contain a stencil centered
 * around that face (i.e. values from \f$ i-2 \f$ to \f$ i+1 \f$). For the right
 * side value `Ur`,
 *
 * \f[
 * U_R \equiv U_{i+\epsilon-1/2} = U_i - \frac{\Delta U}{2}
 * \f]
 *
 * where \f$ \Delta U \f$ is determined by using the @ref ghl_maxmod and @ref ghl_minmod
 * functions. Specifically, it is computed as
 *
 * \f[
 * \Delta U = \mathrm{maxmod}(\sigma_1,\ \sigma_2)
 * \f]
 *
 * where
 *
 * \f[
 * \sigma_1 = \mathrm{minmod}\left( U_{i} - U_{i-1},\ 2(U_{i+1} - U_{i}) \right)
 * \f]
 *
 * and
 *
 * \f[
 * \sigma_2 = \mathrm{minmod}\left( 2(U_{i} - U_{i-1}),\ U_{i+1} - U_{i} \right)
 * \f]
 *
 * Similarly, the left side value `Ul`
 *
 * \f[
 * U_L \equiv U_{i-\epsilon-1/2} = U_{i-1} + \frac{\Delta U}{2}
 * \f]
 *
 * uses
 *
 * \f[
 * \Delta U = \mathrm{maxmod}(\sigma_1,\ \sigma_2)
 * \f]
 *
 * where
 *
 * \f[
 * \sigma_1 = \mathrm{minmod}\left( U_{i-1} - U_{i-2},\ 2(U_{i} - U_{i-1}) \right)
 * \f]
 *
 * and
 *
 * \f[
 * \sigma_2 = \mathrm{minmod}\left( 2(U_{i-1} - U_{i-2}),\ U_{i} - U_{i-1} \right)
 * \f]
 *
 * @param[in] U:   1D array containing values of variable \f$ U \f$
 *
 * @param[out] Ur: pointer to a double; set to the value of the right side of the face
 *
 * @param[out] Ul: pointer to a double; set to the value of the left side of the face
 *
 * @returns void
 */
void ghl_superbee_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul) {

  const double Um1m2 = U[1] - U[0];
  const double U0m1  = U[2] - U[1];
  const double Up10  = U[3] - U[2];

  const double sigma_i   = ghl_maxmod(ghl_minmod(U0m1, 2.0*Up10),
                                      ghl_minmod(2.0*U0m1, Up10));

  const double sigma_im1 = ghl_maxmod(ghl_minmod(Um1m2, 2.0*U0m1),
                                      ghl_minmod(2.0*Um1m2, U0m1));

  *Ur = U[2] - 0.5*sigma_i;
  *Ul = U[1] + 0.5*sigma_im1;
 }
