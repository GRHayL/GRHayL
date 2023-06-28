#ifndef INDUCTION_H_
#define INDUCTION_H_

#include "ghl.h"

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
  double alpha;
  double alpha_Phi_minus_betaj_A_j;
  double sqrtg_Ai[3];
  double betai[3];
} induction_interp_vars;

//------------------ Functions ---------------------

#ifdef __cplusplus
extern "C" {
#endif

double ghl_HLL_2D_flux(
      const double psi6,
      const HLL_2D_vars *restrict vars);

void ghl_interpolate_with_cell_centered_ADM(
      const metric_quantities metric_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      induction_interp_vars *restrict interp_vars);

void ghl_interpolate_with_cell_centered_BSSN(
      const metric_quantities metric[2][2][2],
      const double psi_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      induction_interp_vars *restrict interp_vars);

void ghl_interpolate_with_vertex_centered_ADM(
      const metric_quantities metric_stencil[2][2][2],
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
      const double Ai_stencil[3][2],
      const double phitilde_stencil[3][5]);

//-------------- Helper Functions ------------------

void ghl_A_i_avg(
      const int size,
      const double Ax_stencil[size][size][size],
      const double Ay_stencil[size][size][size],
      const double Az_stencil[size][size][size],
      double A_to_phitilde[3],
      double A_to_Ax[3],
      double A_to_Ay[3],
      double A_to_Az[3]);

void ghl_BSSN_cell_interp(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      const double psi_stencil[size][size][size],
      metric_quantities *restrict metric_interp,
      double lapse_psi2_interp[3],
      double *restrict lapse_over_psi6_interp);

void ghl_ADM_cell_interp(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      metric_quantities *restrict metric_interp,
      double *restrict lapse_over_psi6_interp);

void ghl_ADM_vertex_interp(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      double gammaUU_interp[3][3]);

#ifdef __cplusplus
}
#endif

//--------------------------------------------------

#endif // INDUCTION_H_
