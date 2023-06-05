#ifndef INDUCTION_H_
#define INDUCTION_H_

#include "grhayl.h"

typedef struct em_quantities {
  double phitilde;
  double AD[3];
} em_quantities;

typedef struct HLL_2D_vars {
  double B1r, B1l;
  double B2r, B2l;
  double c1_min, c1_max;
  double c2_min, c2_max;
  double v1rr, v1rl, v1lr, v1ll;
  double v2rr, v2rl, v2lr, v2ll;
} HLL_2D_vars;

typedef struct induction_interp_vars {
  double alpha_interp;
  double alpha_Phi_minus_betaj_A_j_interp;
  double alpha_sqrtg_Ai_interp[3];
  double shifti_interp[3];
} induction_interp_vars;

//--------------------------------------------------

double grhayl_HLL_2D_flux(
      const double psi6,
      const HLL_2D_vars *restrict vars);

#ifdef __cplusplus
extern "C" {
#endif

void grhayl_interpolate_for_induction_rhs(
      const metric_quantities metric[2][2][2],
      const double psi_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      induction_interp_vars *restrict interp_vars);

double grhayl_calculate_phitilde_rhs(
      const double dxi[3],
      const double Lorenz_damping_factor,
      const double alpha_interp,
      const double shiftx_interp[5],
      const double shifty_interp[5],
      const double shiftz_interp[5],
      const double Ai_stencil[3][2],
      const double phitilde_stencil[3][5]);

#ifdef __cplusplus
}
#endif

#endif // INDUCTION_H_
