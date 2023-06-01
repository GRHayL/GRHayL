#ifndef INDUCTION_H_
#define INDUCTION_H_

#include "grhayl.h"

typedef struct em_quantities {
  double phitilde;
  double AD[3];
} em_quantities;

typedef struct A_no_gauge_vars {
  double psi6;
  double B1r, B1l;
  double B2r, B2l;
  double c1_min, c1_max;
  double c2_min, c2_max;
  double v1rr, v1rl, v1lr, v1ll;
  double v2rr, v2rl, v2lr, v2ll;
} A_no_gauge_vars;

typedef struct A_gauge_rhs_vars {
  double dxi[3];
  double alpha_interp;
  double alpha_Phi_minus_betaj_A_j_interp;
  double alpha_sqrtg_Ax_interp[2]; // [i,j,k], [i+1,j,  k  ]
  double alpha_sqrtg_Ay_interp[2]; // [i,j,k], [i,  j+1,k  ]
  double alpha_sqrtg_Az_interp[2]; // [i,j,k], [i,  j,  k+1]
  double phitildex[5], shiftx_interp[5]; // [i-2,j,  k  ], [i-1,j,  k  ], [i,j,k], [i+1,j,  k  ], [i+2,j,  k  ]
  double phitildey[5], shifty_interp[5]; // [i,  j-2,k  ], [i,  j-1,k  ], [i,j,k], [i,  j+1,k  ], [i,  j+2,k  ]
  double phitildez[5], shiftz_interp[5]; // [i,  j,  k-2], [i,  j,  k-1], [i,j,k], [i,  j,  k+1], [i,  j,  k+2]
  double phitilde_rhs, A_x_gauge_rhs, A_y_gauge_rhs, A_z_gauge_rhs;
} A_gauge_rhs_vars;

//--------------------------------------------------

double grhayl_HLL_flux(const A_no_gauge_vars *restrict vars);

#ifdef __cplusplus
extern "C" {
#endif

void grhayl_interpolate_for_A_gauge_rhs(
      const metric_quantities metric[2][2][2],
      const double psi_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      A_gauge_rhs_vars *restrict gauge_rhs_vars);

void grhayl_calculate_phitilde_rhs(
      const double Lorenz_damping_factor,
      A_gauge_rhs_vars *restrict vars);

#ifdef __cplusplus
}
#endif

#endif // INDUCTION_H_
