#ifndef INDUCTION_GEM_H_
#define INDUCTION_GEM_H_

#include "GRHayL.h"

typedef struct induction_lr {
  double B1r, B1l;
  double B2r, B2l;
  double c1_min, c1_max;
  double c2_min, c2_max;
  double v1rr, v1rl, v1lr, v1ll;
  double v2rr, v2rl, v2lr, v2ll;
} induction_lr;

typedef struct induction_gauge {
  double gupxx[2][2][2];
  double gupxy[2][2][2];
  double gupxz[2][2][2];
  double gupyy[2][2][2];
  double gupyz[2][2][2];
  double gupzz[2][2][2];
  double lapse[2][2][2];
  double psi[2][2][2];
  double shiftx[2][2][2];
  double shifty[2][2][2];
  double shiftz[2][2][2];
  double A_x[3][3][3];
  double A_y[3][3][3];
  double A_z[3][3][3];
  double phitilde;
} induction_gauge;

typedef struct induction_gauge_rhs {
  double dxi[3];
  double alpha_interp;
  double alpha_Phi_minus_betaj_A_j_interp[4]; // [i,j,k], [i,j-1,k], [i,j-1,k], [i,j,k-1]
  double alpha_sqrtg_Ax_interp[2]; // [i,j,k], [i+1,j,  k  ]
  double alpha_sqrtg_Ay_interp[2]; // [i,j,k], [i,  j+1,k  ]
  double alpha_sqrtg_Az_interp[2]; // [i,j,k], [i,  j,  k+1]
  double phitildex[5], shiftx_interp[5]; // [i-2,j,  k  ], [i-1,j,  k  ], [i,j,k], [i+1,j,  k  ], [i+2,j,  k  ]
  double phitildey[5], shifty_interp[5]; // [i,  j-2,k  ], [i,  j-1,k  ], [i,j,k], [i,  j+1,k  ], [i,  j+2,k  ]
  double phitildez[5], shiftz_interp[5]; // [i,  j,  k-2], [i,  j,  k-1], [i,j,k], [i,  j,  k+1], [i,  j,  k+2]
  double phitilde_rhs, A_x_gauge_rhs, A_y_gauge_rhs, A_z_gauge_rhs;
} induction_gauge_rhs;

//typedef struct A_to_B_stag {
//  double dxi, dyi, dzi;
//  double Bx_stagger, By_stagger, Bz_stagger;
//  double Ax_j[2], Ax_k[2];
//  double Ay_i[2], Ay_k[2];
//  double Az_i[2], Az_j[2];
//  double psi[4];
//} A_to_B_stag;
//
//typedef struct B_stag_to_B {
//  double Bx_stagger[2], By_stagger[2], Bz_stagger[2];
//  double Bx, By, Bz;
//} B_stag_to_B;

//--------------------------------------------------

double A_i_rhs_no_gauge_terms(
             const double psi6,
             induction_lr *restrict vars );

void interpolate_for_A_i_rhs(
             induction_gauge *restrict gauge_vars,
             induction_gauge_rhs *restrict gauge_rhs_vars );

void calculate_phitilde_and_A_i_rhs(
             const double Lorenz_damping_factor,
             induction_gauge_rhs *restrict vars);

//void compute_Bstagger_from_A(A_to_B_stag *restrict A_to_B);
//
//void compute_B_from_Bstagger(B_stag_to_B *restrict Bs_to_B);

#endif // INDUCTION_GEM_H_
