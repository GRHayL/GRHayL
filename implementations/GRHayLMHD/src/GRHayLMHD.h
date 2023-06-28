#ifndef GRHAYLMHD_H_
#define GRHAYLMHD_H_

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in evaluate_MHD_rhs.C.
enum recon_indices{
      BX_STAGGER, BY_STAGGER, BZ_STAGGER,
      VXR, VYR, VZR, VXL,VYL, VZL, MAXNUMVARS};

//Interpolates to the +1/2 face
#define interpolate_to_face(Vm1,V0,Vp1,Vp2) (0.5625*(V0 + Vp1) - 0.0625*(Vm1+Vp2))
//Interpolates to the -1/2 face
#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

// Computes derivative factor dx*deriv
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) ((A0 - AM1)*(Varp1 - Varm1) + AM1*(Varp2 - Varm2))

void GRHayLMHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      metric_quantities *restrict metric);

void GRHayLMHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const CCTK_REAL dxi,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      metric_quantities *restrict metric_derivs);

void GRHayLMHD_set_symmetry_gzs_staggered(
      const cGH *cctkGH,
      const CCTK_REAL *X,
      const CCTK_REAL *Y,
      const CCTK_REAL *Z,
      CCTK_REAL *gridfunc,
      const CCTK_REAL *gridfunc_syms,
      const int stagger_x,  //TODO: unused
      const int stagger_y,  //TODO: unused
      const int stagger_z);

/******** Helper functions for the RHS calculations *************/

void GRHayLMHD_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const double *rho_b,
                         const double *pressure,
                         const double *v_flux,
                         const CCTK_REAL **in_prims,
                         CCTK_REAL **out_prims_r,
                         CCTK_REAL **out_prims_l);

void GRHayLMHD_calculate_MHD_dirn_rhs(
      const cGH *cctkGH,
      const int flux_dirn,
      const CCTK_REAL *restrict dX,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      const double *rho,
      const double *pressure,
      const double *vx,
      const double *vy,
      const double *vz,
      const CCTK_REAL **B_center,
      const double *B_stagger,
      /*const*/ CCTK_REAL **in_prims_r,
      /*const*/ CCTK_REAL **in_prims_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax,
      CCTK_REAL *restrict rho_star_flux,
      CCTK_REAL *restrict tau_flux,
      CCTK_REAL *restrict Stildex_flux,
      CCTK_REAL *restrict Stildey_flux,
      CCTK_REAL *restrict Stildez_flux,
      CCTK_REAL *restrict rho_star_rhs,
      CCTK_REAL *restrict tau_rhs,
      CCTK_REAL *restrict Stildex_rhs,
      CCTK_REAL *restrict Stildey_rhs,
      CCTK_REAL *restrict Stildez_rhs);

// The const are commented out because C does not support implicit typecasting of types when
// they are more than 1 level removed from the top pointer. i.e. I can pass the argument with
// type "CCTK_REAL *" for an argument expecting "const CCTK_REAL *" because this is only 1 level
// down (pointer to CCTK_REAL -> pointer to const CCTK_REAL). It will not do
// pointer to pointer to CCTK_REAL -> pointer to pointer to const CCTK_REAL. I saw comments
// suggesting this may become part of the C23 standard, so I guess you can uncomment this
// in like 10 years.
void GRHayLMHD_A_flux_rhs(
      const cGH *restrict cctkGH,
      const int A_dir,
      /*const*/ CCTK_REAL **out_prims_r,
      /*const*/ CCTK_REAL **out_prims_l,
      /*const*/ CCTK_REAL **cmin,
      /*const*/ CCTK_REAL **cmax,
      CCTK_REAL *restrict A_rhs);

/****************************************************************/
#endif // GRHAYLMHD_H_
