#ifndef GRHAYL_IGH_H_
#define GRHAYL_IGH_H_

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
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
#define AM1 -0.0625
#define A0   0.5625
#define A1   0.5625
#define A2  -0.0625
#define COMPUTE_FCVAL(Varm1,Var,Varp1,Varp2) (AM1*(Varm1) + A0*(Var) + A1*(Varp1) + A2*(Varp2))

void GRHayL_IGH_interpolate_to_face_and_initialize_metric(
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

/******** Helper functions for the RHS calculations *************/

void GRHayL_IGH_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const eos_parameters *restrict eos,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l);

void GRHayL_IGH_reconstruction_loop_no_rho_P(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const eos_parameters *restrict eos,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l);

void GRHayL_IGH_calculate_tau_source_rhs(
      const cGH *cctkGH,
      const eos_parameters *restrict eos,
      const double **in_metric,
      const double **in_curv,
      const double **in_prims,
      double *restrict tau_rhs);

void GRHayL_IGH_compute_characteristic_speeds(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const eos_parameters *restrict eos,
      const double **in_metric,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      double *restrict cmin,
      double *restrict cmax);

void GRHayL_IGH_calculate_HD_dirn_rhs(
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

#endif // GRHAYL_IGH_H_
