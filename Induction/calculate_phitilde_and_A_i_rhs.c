#include "induction_gem.h"

void calculate_phitilde_and_A_i_rhs(
                                   const double Lorenz_damping_factor,
                                   induction_gauge_rhs *restrict vars) {

  // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
  //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
  // - partial_i -> - (A_{i} - A_{i-1})/dX = (A_{i-1} - A_{i})/dX, for Ax
  vars->A_x_gauge_rhs = vars->dx[0]*(vars->alpha_Phi_minus_betaj_A_j_interp[1] - vars->alpha_Phi_minus_betaj_A_j_interp[0]);
  vars->A_y_gauge_rhs = vars->dx[1]*(vars->alpha_Phi_minus_betaj_A_j_interp[2] - vars->alpha_Phi_minus_betaj_A_j_interp[0]);
  vars->A_z_gauge_rhs = vars->dx[2]*(vars->alpha_Phi_minus_betaj_A_j_interp[3] - vars->alpha_Phi_minus_betaj_A_j_interp[0]);

  // \partial_t psi6phi = [shift advection term] + \partial_j (\alpha \sqrt{\gamma} A^j)
  // Here we compute [shift advection term] = \partial_j (\beta^j psi6phi)
  // Cache misses are likely more expensive than branch mispredictions here,
  //       which is why we use if() statements and array lookups inside the if()'s.
  vars->phitilde_rhs = 0.0;

  // \partial_x (\beta^x psi6phi) :
  if(vars->shiftx_interp[2] < 0.0) {
    vars->phitilde_rhs+=0.5*vars->dx[0]*(+    vars->shiftx_interp[0]*vars->phitildex[0]
                                         -4.0*vars->shiftx_interp[1]*vars->phitildex[1]
                                         +3.0*vars->shiftx_interp[2]*vars->phitildex[2]);
  } else {
    vars->phitilde_rhs+=0.5*vars->dx[0]*(-    vars->shiftx_interp[4]*vars->phitildex[4]
                                         +4.0*vars->shiftx_interp[3]*vars->phitildex[3]
                                         -3.0*vars->shiftx_interp[2]*vars->phitildex[2]);
  }

  // \partial_y (\beta^y psi6phi) :
  if(vars->shifty_interp[2] < 0.0) {
    vars->phitilde_rhs+=0.5*vars->dx[1]*(+    vars->shifty_interp[0]*vars->phitildey[0]
                                         -4.0*vars->shifty_interp[1]*vars->phitildey[1]
                                         +3.0*vars->shifty_interp[2]*vars->phitildey[2]);
  } else {
    vars->phitilde_rhs+=0.5*vars->dx[1]*(-    vars->shifty_interp[4]*vars->phitildey[4]
                                         +4.0*vars->shifty_interp[3]*vars->phitildey[3]
                                         -3.0*vars->shifty_interp[2]*vars->phitildey[2]);
  }

  // \partial_z (\beta^z psi6phi) :
  if(vars->shiftz_interp[2] < 0.0) {
    vars->phitilde_rhs+=0.5*vars->dx[2]*(+    vars->shiftz_interp[0]*vars->phitildez[0]
                                         -4.0*vars->shiftz_interp[1]*vars->phitildez[1]
                                         +3.0*vars->shiftz_interp[2]*vars->phitildez[2]);
  } else {
    vars->phitilde_rhs+=0.5*vars->dx[2]*(-    vars->shiftz_interp[4]*vars->phitildez[4]
                                         +4.0*vars->shiftz_interp[3]*vars->phitildez[3]
                                         -3.0*vars->shiftz_interp[2]*vars->phitildez[2]);
  }

  // Next we add \partial_j (\alpha \sqrt{\gamma} A^j) to \partial_t psi6phi:
  vars->phitilde_rhs+=vars->dx[0]*(vars->alpha_sqrtg_Ax_interp[0] - vars->alpha_sqrtg_Ax_interp[1])
                    + vars->dx[1]*(vars->alpha_sqrtg_Ay_interp[0] - vars->alpha_sqrtg_Ay_interp[1])
                    + vars->dx[2]*(vars->alpha_sqrtg_Az_interp[0] - vars->alpha_sqrtg_Az_interp[1]);

  // *GENERALIZED* LORENZ GAUGE:
  // Finally, add damping factor to \partial_t psi6phi
  //subtract lambda * alpha psi^6 Phi
  vars->phitilde_rhs+=-Lorenz_damping_factor*vars->alpha_interp*vars->phitildex[2];
}
