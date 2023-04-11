#ifndef GRHAYL_IGM_H_
#define GRHAYL_IGM_H_
#include "cctk.h"
#include "GRHayLET.h"

enum metric_indices{
      LAPSE, BETAX, BETAY, BETAZ,
      GXX, GXY, GXZ, GYY, GYZ, GZZ};

enum ext_curv{
      KXX, KXY, KXZ, KYY, KYZ, KZZ};

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in evaluate_MHD_rhs.C.
enum recon_indices{
      RHOB, PRESSURE, VX, VY, VZ,
      BX_CENTER, BY_CENTER, BZ_CENTER,
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

void GRHayL_IGM_convert_ADM_to_BSSN(const cGH *cctkGH,
      double *gxx, double *gxy, double *gxz,
      double *gyy, double *gyz, double *gzz,
      double *phi, double *psi,
      double *gtxx, double *gtxy, double *gtxz,
      double *gtyy, double *gtyz, double *gtzz,
      double *gtupxx, double *gtupxy, double *gtupxz,
      double *gtupyy, double *gtupyz, double *gtupzz);

void GRHayL_IGM_interpolate_to_face_and_initialize_metric(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dirn,
      const double *restrict lapse,
      const double *restrict betax,
      const double *restrict betay,
      const double *restrict betaz,
      const double *restrict gxx,
      const double *restrict gxy,
      const double *restrict gxz,
      const double *restrict gyy,
      const double *restrict gyz,
      const double *restrict gzz,
      metric_quantities *restrict metric);

void GRHayL_IGM_set_symmetry_gzs_staggered(
      const cGH *cctkGH,
      const double *X,
      const double *Y,
      const double *Z,
      double *gridfunc,
      const double *gridfunc_syms,
      const int stagger_x,  //TODO: unused
      const int stagger_y,  //TODO: unused
      const int stagger_z);

/******** Helper functions for the RHS calculations *************/

void GRHayL_IGM_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const eos_parameters *restrict eos,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l);

void GRHayL_IGM_reconstruction_loop_no_rho_P(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const eos_parameters *restrict eos,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l);

void GRHayL_IGM_calculate_tau_source_rhs(
      const cGH *cctkGH,
      const eos_parameters *restrict eos,
      const double **in_metric,
      const double **in_curv,
      const double **in_prims,
      double *restrict tau_rhs);

void GRHayL_IGM_compute_characteristic_speeds(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const eos_parameters *restrict eos,
      const double **in_metric,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      double *restrict cmin,
      double *restrict cmax);

void GRHayL_IGM_calculate_MHD_dirn_rhs(
      const cGH *cctkGH,
      const int flux_dirn,
      const double *restrict dX,
      const eos_parameters *restrict eos,
      const double **in_metric,
      const double **in_prims,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      const double *restrict cmin,
      const double *restrict cmax,
      double *restrict rho_star_flux,
      double *restrict tau_flux,
      double *restrict Stildex_flux,
      double *restrict Stildey_flux,
      double *restrict Stildez_flux,
      double *restrict rho_star_rhs,
      double *restrict tau_rhs,
      double *restrict Stildex_rhs,
      double *restrict Stildey_rhs,
      double *restrict Stildez_rhs);

// The const are commented out because C does not support implicit typecasting of types when
// they are more than 1 level removed from the top pointer. i.e. I can pass the argument with
// type "double *" for an argument expecting "const double *" because this is only 1 level
// down (pointer to double -> pointer to const double). It will not do
// pointer to pointer to double -> pointer to pointer to const double. I saw comments
// suggesting this may become part of the C23 standard, so I guess you can uncomment this
// in like 10 years.
void GRHayL_IGM_A_no_gauge_rhs(
      const cGH *restrict cctkGH,
      const int A_dir,
      /*const*/ double **out_prims_r,
      /*const*/ double **out_prims_l,
      const double *restrict phi_bssn,
      /*const*/ double **cmin,
      /*const*/ double **cmax,
      double *restrict A_rhs);

void GRHayL_IGM_phitilde_and_A_gauge_rhs(
      const cGH *cctkGH,
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
