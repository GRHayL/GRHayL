#ifndef INDUCTION_H_
#define INDUCTION_H_

#include "ghl.h"

typedef struct HLL_2D_vars {
  double B1r, B1l;
  double B2r, B2l;
  double c1_min, c1_max;
  double c2_min, c2_max;
  double v1rr, v1rl, v1lr, v1ll;
  double v2rr, v2rl, v2lr, v2ll;
} HLL_2D_vars;

typedef struct induction_interp_vars {
  double alpha;
  double alpha_Phi_minus_betaj_A_j;
  double sqrtg_Ai[3];
  double betai[3];
} induction_interp_vars;

//------------------ Functions ---------------------

#ifdef __cplusplus
extern "C" {
#endif

double ghl_HLL_2D_flux_with_B(
      const double psi6,
      const HLL_2D_vars *restrict vars);

double ghl_HLL_2D_flux_with_Btilde(
      const HLL_2D_vars *restrict vars);

void ghl_interpolate_with_cell_centered_ADM(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      induction_interp_vars *restrict interp_vars);

void ghl_interpolate_with_cell_centered_BSSN(
      const ghl_metric_quantities metric[2][2][2],
      const double psi_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      induction_interp_vars *restrict interp_vars);

void ghl_interpolate_with_vertex_centered_ADM(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      induction_interp_vars *restrict interp_vars);

double ghl_calculate_phitilde_rhs(
      const double dxi[3],
      const double Lorenz_damping_factor,
      const double alpha_interp,
      const double shiftx_interp[5],
      const double shifty_interp[5],
      const double shiftz_interp[5],
      const double sqrtg_Ai_stencil[3][2],
      const double phitilde_stencil[3][5]);

#ifdef __cplusplus
}
#endif

//--------------------------------------------------

#endif // INDUCTION_H_
