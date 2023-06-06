#include "induction.h"

/* Function    : ghl_calculate_phitilde_rhs()
 * Description : compute RHS for \tilde{phi} = psi^6*Phi and the gauge contribution
 *               for A_i
 *
 * Inputs      : Lorenz_damping_factor - TODO: ?
 *             : vars                  - A_gauge_rhs_vars struct containing the
 *                                       stencils for interpolated variables;
 *                                       these can be computed using
 *                                       ghl_interpolate_for_A_gauge_rhs()
 *
 * Outputs     : vars                  - struct elements A_i_gauge_rhs and
 *                                       phitilde_rhs contain RHS values
 *
 */
double ghl_calculate_phitilde_rhs(
      const double dxi[3],
      const double Lorenz_damping_factor,
      const double alpha_interp,
      const double shiftx_interp[5],
      const double shifty_interp[5],
      const double shiftz_interp[5],
      const double Ai_stencil[3][2],
      const double phitilde_stencil[3][5]) {

  // \partial_t psi6phi = [shift advection term] + \partial_j (\alpha \sqrt{\gamma} A^j)
  // Here we compute [shift advection term] = \partial_j (\beta^j psi6phi)
  // TODO: IllinoisGRMHD used this logic in the hope of avoiding cache misses. However,
  // the necessity of passing stencils means that this has no benefit now. This
  // needs to be revisited.

  double phitilde_rhs = 0.0;

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
  // - \lambda * \alpha \psi^6 \Phi
  phitilde_rhs += dxi[0]*(Ai_stencil[0][0] - Ai_stencil[0][1])
                + dxi[1]*(Ai_stencil[1][0] - Ai_stencil[1][1])
                + dxi[2]*(Ai_stencil[2][0] - Ai_stencil[2][1])
                - Lorenz_damping_factor*alpha_interp*phitilde_stencil[0][2];

  return phitilde_rhs;
}
