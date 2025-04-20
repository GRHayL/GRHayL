#include "ghl_induction.h"

/**
 * @ingroup mag_gauge
 * @brief Compute RHS for \f$ \tilde{\Phi} = \psi^6\Phi \f$
 *
 * @todo
 * IllinoisGRMHD claimed that the `if` logic for the advection term
 * was to avoid cache misses. However, this seems like would have no
 * benefit for @grhayl given the explicit passing of stencils. This
 * should be revisited to check if a better implementation exists
 * for @grhayl.
 *
 * @details
 * This function computes
 *
 * \f[
 * \tilde{\Phi}^\mathrm{RHS} = -\partial_i \left( \alpha \sqrt{\gamma} A^i - \beta^i\tilde{\Phi} \right)
 *                             - \lambda \alpha \tilde{\Phi}
 * \f]
 *
 * using interpolated values of the metric and vector potential. These
 * quantities can be computed using the provided interpolators.
 * The vector potential contribution is just a standard finite difference
 * stencil
 *
 * \f[
 * -\frac{\alpha\sqrt{\gamma}}{dx^i} \partial_i A^i = \alpha\sqrt{\gamma} \left( \frac{A^1_{(i,j,k)} - A^1_{(i+1,j,k)}}{\Delta x^1}
 *                                                                             + \frac{A^2_{(i,j,k)} - A^2_{(i,j+1,k)}}{\Delta x^2}
 *                                                                             + \frac{A^3_{(i,j,k)} - A^3_{(i,j,k+1)}}{\Delta x^3} \right)
 *  \f]
 *
 *  while the shift term is instead upwinded in a similar manner to
 *  the usual BSSN treatment. For example,
 *
 *  \f[
 *  \partial_1\beta^1\tilde{\Phi} = \begin{cases} 
 *            \frac{1}{2\Delta x^1} \left( \beta^1_{(i-2,j,k)}\tilde{\Phi}_{(i-2,j,k)}
 *                                       - 4\beta^1_{(i-1,j,k)}\tilde{\Phi}_{(i-1,j,k)}
 *                                       + 3\beta^1_{(i,j,k)}\tilde{\Phi}_{(i,j,k)} \right) \ , & \beta^1_{(i,j,k)} < 0 \\
 *            \frac{1}{2\Delta x^1} \left( -\beta^1_{(i+2,j,k)}\tilde{\Phi}_{(i+2,j,k)}
 *                                       + 4\beta^1_{(i+1,j,k)}\tilde{\Phi}_{(i+1,j,k)}
 *                                       - 3\beta^1_{(i,j,k)}\tilde{\Phi}_{(i,j,k)} \right) \ , & \beta^1_{(i,j,k)} \ge 0
 * \end{cases}
 * \f]
 *
 * The final term arises from the generalized Lorenz gauge, and
 * can be directly computed without any stencils. The quantity
 * \f$ \lambda \f$ is the input quantity `Lorenz_damping_factor`.
 * To ensure consistency within the code, we recommend using
 * the stored value within an instance of the ghl_parameters struct.
 *
 * @param[in] dxi:                   1D array of the inverse grid spacing
 *
 * @param[in] Lorenz_damping_factor: sets the damping factor for the generalized Lorenz gauge term
 *
 * @param[in] alpha_interp:          value of the lapse at the location of \f$ \tilde{\Phi} \f$
 *
 * @param[in] shiftx_interp:         1D stencil array of \f$ \beta^x \f$ (on the same grid as \f$ \tilde{\Phi} \f$)
 *                                   in the x direction from \f$ i-2 \f$ to \f$ i+2 \f$
 *
 * @param[in] shifty_interp:         1D stencil array of \f$ \beta^y \f$ (on the same grid as \f$ \tilde{\Phi} \f$)
 *                                   in the y direction from \f$ i-2 \f$ to \f$ i+2 \f$
 *
 * @param[in] shiftz_interp:         1D stencil array of \f$ \beta^z \f$ (on the same grid as \f$ \tilde{\Phi} \f$)
 *                                   in the z direction from \f$ i-2 \f$ to \f$ i+2 \f$
 *
 * @param[in] sqrtg_Ai_stencil:      2D stencil array of \f$ \sqrt{-g} A_i \f$. The first dimension represents the
 *                                   spatial component \f$ i \f$, and the second dimension represents the stencil.
 *                                   The stencil goes from \f$ j \f$ to \f$ j+1 \f$ in the direction of the spatial component.
 *
 * @param[in] phitilde_stencil:      2D stencil array of \f$ \tilde{\Phi} \f$. The first dimension represents the
 *                                   direction of the stencil, and the second dimension represents the stencil itself.
 *                                   The stencil goes from \f$ i-2 \f$ to \f$ i+2 \f$.
 *
 * @returns value of \f$ \tilde{\Phi}^\mathrm{RHS} \f$
 */
double ghl_calculate_phitilde_rhs(
      const double dxi[3],
      const double Lorenz_damping_factor,
      const double alpha_interp,
      const double shiftx_interp[5],
      const double shifty_interp[5],
      const double shiftz_interp[5],
      const double sqrtg_Ai_stencil[3][2],
      const double phitilde_stencil[3][5]) {

  double phitilde_rhs = 0.0;

  // Here we compute [shift advection term] = \partial_j (\beta^j psi6phi)
  // \partial_x (\beta^x psi6phi) :
  if(shiftx_interp[2] < 0.0) {
    phitilde_rhs += 0.5*dxi[0]*(+   shiftx_interp[0]*phitilde_stencil[0][0]
                               -4.0*shiftx_interp[1]*phitilde_stencil[0][1]
                               +3.0*shiftx_interp[2]*phitilde_stencil[0][2]);
  } else {
    phitilde_rhs += 0.5*dxi[0]*(-   shiftx_interp[4]*phitilde_stencil[0][4]
                               +4.0*shiftx_interp[3]*phitilde_stencil[0][3]
                               -3.0*shiftx_interp[2]*phitilde_stencil[0][2]);
  }

  // \partial_y (\beta^y psi6phi) :
  if(shifty_interp[2] < 0.0) {
    phitilde_rhs += 0.5*dxi[1]*(+   shifty_interp[0]*phitilde_stencil[1][0]
                               -4.0*shifty_interp[1]*phitilde_stencil[1][1]
                               +3.0*shifty_interp[2]*phitilde_stencil[1][2]);
  } else {
    phitilde_rhs += 0.5*dxi[1]*(-   shifty_interp[4]*phitilde_stencil[1][4]
                               +4.0*shifty_interp[3]*phitilde_stencil[1][3]
                               -3.0*shifty_interp[2]*phitilde_stencil[1][2]);
  }

  // \partial_z (\beta^z psi6phi) :
  if(shiftz_interp[2] < 0.0) {
    phitilde_rhs += 0.5*dxi[2]*(+   shiftz_interp[0]*phitilde_stencil[2][0]
                               -4.0*shiftz_interp[1]*phitilde_stencil[2][1]
                               +3.0*shiftz_interp[2]*phitilde_stencil[2][2]);
  } else {
    phitilde_rhs += 0.5*dxi[2]*(-   shiftz_interp[4]*phitilde_stencil[2][4]
                               +4.0*shiftz_interp[3]*phitilde_stencil[2][3]
                               -3.0*shiftz_interp[2]*phitilde_stencil[2][2]);
  }

  // Next we add \partial_j (\alpha \sqrt{\gamma} A^j)
  // and add the damping factor for the generalized Lorenz gauge
  phitilde_rhs += dxi[0]*(sqrtg_Ai_stencil[0][0] - sqrtg_Ai_stencil[0][1])
                + dxi[1]*(sqrtg_Ai_stencil[1][0] - sqrtg_Ai_stencil[1][1])
                + dxi[2]*(sqrtg_Ai_stencil[2][0] - sqrtg_Ai_stencil[2][1])
                - Lorenz_damping_factor*alpha_interp*phitilde_stencil[0][2];

  return phitilde_rhs;
}
