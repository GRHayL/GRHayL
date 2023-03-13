#ifndef GRHAYL_IGM_H_
#define GRHAYL_IGM_H_
#include "cctk.h"
#include "GRHayLET.h"

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in GRHayLET_evaluate_MHD_rhs.C.
static const int RHOB=0,PRESSURE=1,VX=2,VY=3,VZ=4,
  BX_CENTER=5,BY_CENTER=6,BZ_CENTER=7,BX_STAGGER=8,BY_STAGGER=9,BZ_STAGGER=10,
  VXR=11,VYR=12,VZR=13,VXL=14,VYL=15,VZL=16,MAXNUMVARS=17;  //<-- Be _sure_ to define MAXNUMVARS appropriately!

/******** Helper functions for the RHS calculations *************/

void reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const eos_parameters *restrict eos,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l);

void reconstruction_loop_no_rho_P(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const eos_parameters *restrict eos,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l);

// The const are commented out because C does not support implicit typecasting of types when
// they are more than 1 level removed from the top pointer. i.e. I can pass the argument with
// type "double *" for an argument expecting "const double *" because this is only 1 level
// down (pointer to double -> pointer to const double). It will not do
// pointer to pointer to double -> pointer to pointer to const double. I saw comments
// suggesting this may become part of the C23 standard, so I guess you can uncomment this
// in like 10 years.
void A_no_gauge_rhs(const cGH *restrict cctkGH, const int A_dir,
               /*const*/ double **out_prims_r,
               /*const*/ double **out_prims_l,
               const double *restrict phi_bssn,
               const double **cmin,
               const double **cmax,
               double *restrict A_rhs);

void phitilde_and_A_gauge_rhs(const cGH *cctkGH,
                   const double *restrict dX,
                   const double *restrict gupxx,
                   const double *restrict gupxy,
                   const double *restrict gupxz,
                   const double *restrict gupyy,
                   const double *restrict gupyz,
                   const double *restrict gupzz,
                   const double *restrict psi,
                   const double *restrict lapm1,
                   const double *restrict betax,
                   const double *restrict betay,
                   const double *restrict betaz,
                   const double *restrict Ax,
                   const double *restrict Ay,
                   const double *restrict Az,
                   const double *restrict phitilde,
                   double Lorenz_damping_factor,
                   double *restrict shiftx_interp,
                   double *restrict shifty_interp,
                   double *restrict shiftz_interp,
                   double *restrict alpha_interp,
                   double *restrict alpha_Phi_minus_betaj_A_j_interp,
                   double *restrict alpha_sqrtg_Ax_interp,
                   double *restrict alpha_sqrtg_Ay_interp,
                   double *restrict alpha_sqrtg_Az_interp,
                   double *restrict phitilde_rhs,
                   double *restrict Ax_rhs,
                   double *restrict Ay_rhs,
                   double *restrict Az_rhs);

/****************************************************************/
#endif // GRHAYL_IGM_H_
