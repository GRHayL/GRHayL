#ifndef GRHAYL_IGH_H_
#define GRHAYL_IGH_H_
#include "cctk.h"
#include "GRHayLET.h"

static const int LAPSE=0, BETAX=1,BETAY=2, BETAZ=3,
                 GXX=4, GXY=5, GXZ=6,
                 GYY=7, GYZ=8, GZZ=9;

static const int KXX=0, KXY=1, KXZ=2, KYY=3, KYZ=4, KZZ=5;

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in evaluate_MHD_rhs.C.
static const int RHOB=0,PRESSURE=1,VX=2,VY=3,VZ=4,MAXNUMVARS=5;  //<-- Be _sure_ to define MAXNUMVARS appropriately!

void GRHayL_convert_ADM_to_BSSN(const cGH *cctkGH,
      double *gxx, double *gxy, double *gxz,
      double *gyy, double *gyz, double *gzz,
      double *phi, double *psi,
      double *gtxx, double *gtxy, double *gtxz,
      double *gtyy, double *gtyz, double *gtzz,
      double *gtupxx, double *gtupxy, double *gtupxz,
      double *gtupyy, double *gtupyz, double *gtupzz);

void interpolate_to_face_and_initialize_metric(
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

void calculate_tau_source_rhs(
      const cGH *cctkGH,
      const eos_parameters *restrict eos,
      const double **in_metric,
      const double **in_curv,
      const double **in_prims,
      double *restrict tau_rhs);

void compute_characteristic_speeds(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const eos_parameters *restrict eos,
      const double **in_metric,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      double *restrict cmin,
      double *restrict cmax);

void calculate_HD_dirn_rhs(
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
