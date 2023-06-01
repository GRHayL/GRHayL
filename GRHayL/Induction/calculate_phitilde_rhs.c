#include "induction.h"

/* Function    : grhayl_calculate_phitilde_rhs()
 * Description : compute RHS for \tilde{phi} = psi^6*Phi and the gauge contribution
 *               for A_i
 *
 * Inputs      : Lorenz_damping_factor - TODO: ?
 *             : vars                  - A_gauge_rhs_vars struct containing the
 *                                       stencils for interpolated variables;
 *                                       these can be computed using
 *                                       grhayl_interpolate_for_A_gauge_rhs()
 *
 * Outputs     : vars                  - struct elements A_i_gauge_rhs and
 *                                       phitilde_rhs contain RHS values
 *
 */
void grhayl_calculate_phitilde_rhs(
                                   const double Lorenz_damping_factor,
                                   A_gauge_rhs_vars *restrict vars) {

  // \partial_t psi6phi = [shift advection term] + \partial_j (\alpha \sqrt{\gamma} A^j)
  // Here we compute [shift advection term] = \partial_j (\beta^j psi6phi)
  // Cache misses are likely more expensive than branch mispredictions here,
  //       which is why we use if() statements and array lookups inside the if()'s.
  vars->phitilde_rhs = 0.0;

  // \partial_x (\beta^x psi6phi) :
  if(vars->shiftx_interp[2] < 0.0) {
    vars->phitilde_rhs += 0.5*vars->dxi[0]*(+   vars->shiftx_interp[0]*vars->phitildex[0]
                                           -4.0*vars->shiftx_interp[1]*vars->phitildex[1]
                                           +3.0*vars->shiftx_interp[2]*vars->phitildex[2]);
  } else {
    vars->phitilde_rhs += 0.5*vars->dxi[0]*(-   vars->shiftx_interp[4]*vars->phitildex[4]
                                           +4.0*vars->shiftx_interp[3]*vars->phitildex[3]
                                           -3.0*vars->shiftx_interp[2]*vars->phitildex[2]);
  }

  // \partial_y (\beta^y psi6phi) :
  if(vars->shifty_interp[2] < 0.0) {
    vars->phitilde_rhs += 0.5*vars->dxi[1]*(+   vars->shifty_interp[0]*vars->phitildey[0]
                                           -4.0*vars->shifty_interp[1]*vars->phitildey[1]
                                           +3.0*vars->shifty_interp[2]*vars->phitildey[2]);
  } else {
    vars->phitilde_rhs += 0.5*vars->dxi[1]*(-   vars->shifty_interp[4]*vars->phitildey[4]
                                           +4.0*vars->shifty_interp[3]*vars->phitildey[3]
                                           -3.0*vars->shifty_interp[2]*vars->phitildey[2]);
  }

  // \partial_z (\beta^z psi6phi) :
  if(vars->shiftz_interp[2] < 0.0) {
    vars->phitilde_rhs += 0.5*vars->dxi[2]*(+   vars->shiftz_interp[0]*vars->phitildez[0]
                                           -4.0*vars->shiftz_interp[1]*vars->phitildez[1]
                                           +3.0*vars->shiftz_interp[2]*vars->phitildez[2]);
  } else {
    vars->phitilde_rhs += 0.5*vars->dxi[2]*(-   vars->shiftz_interp[4]*vars->phitildez[4]
                                           +4.0*vars->shiftz_interp[3]*vars->phitildez[3]
                                           -3.0*vars->shiftz_interp[2]*vars->phitildez[2]);
  }

  // Next we add \partial_j (\alpha \sqrt{\gamma} A^j) to \partial_t psi6phi:
  vars->phitilde_rhs += vars->dxi[0]*(vars->alpha_sqrtg_Ax_interp[0] - vars->alpha_sqrtg_Ax_interp[1])
                      + vars->dxi[1]*(vars->alpha_sqrtg_Ay_interp[0] - vars->alpha_sqrtg_Ay_interp[1])
                      + vars->dxi[2]*(vars->alpha_sqrtg_Az_interp[0] - vars->alpha_sqrtg_Az_interp[1]);

  // *GENERALIZED* LORENZ GAUGE:
  // Finally, add damping factor to \partial_t psi6phi
  //subtract lambda * alpha psi^6 Phi
  vars->phitilde_rhs -= Lorenz_damping_factor*vars->alpha_interp*vars->phitildex[2];
}
