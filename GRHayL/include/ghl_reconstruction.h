#ifndef GHL_RECONSTRUCTION_H_
#define GHL_RECONSTRUCTION_H_

#include "ghl.h"

// Integer constants to keep track of stencil.
enum reconstruction_stencil {
  MINUS2, MINUS1,
  PLUS_0, PLUS_1,
  PLUS_2
};

#ifdef __cplusplus
extern "C" {
#endif

// PPM functions
void ghl_ppm_reconstruction(
      const double ftilde[2],
      const double var_data[6],
      double *restrict var_datar,
      double *restrict var_datal);

void ghl_ppm_reconstruction_with_steepening(
      const ghl_parameters *restrict params,
      const double pressure[6],
      const double Gamma_eff,
      const double ftilde[2],
      const double var_data[6],
      double *restrict var_datar,
      double *restrict var_datal);

void ghl_ppm_compute_for_cell(
      const double ftilde,
      const double U[5],
      double *restrict Ur_ptr,
      double *restrict Ul_ptr);

void ghl_ppm_compute_for_cell_with_steepening(
      const ghl_parameters *restrict params,
      const double pressure[5],
      const double Gamma_eff,
      const double ftilde,
      const double U[5],
      double *restrict Ur_ptr,
      double *restrict Ul_ptr);

void ghl_compute_ftilde(
      const ghl_parameters *restrict params,
      const double pressure[6],
      const double v_flux_dirn[6],
      double ftilde[2]);

double ghl_shock_detection_ftilde(
      const ghl_parameters *restrict params,
      const double P[5],
      const double v_flux_dirn[5]);

void ghl_steepen_var(
      const ghl_parameters *restrict params,
      const double rho[5],
      const double P[5],
      const double Gamma_eff,
      double *restrict rhor,
      double *restrict rhol);

double ghl_slope_limit(
      const double dU,
      const double dUp1);

// PLM functions
double ghl_minmod(
      const double a,
      const double b);

double ghl_maxmod(
      const double a,
      const double b);

void ghl_minmod_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul);

void ghl_mc_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul);

void ghl_superbee_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul);

void ghl_wenoz_reconstruction(
      const double U[6],
      double *restrict Ur,
      double *restrict Ul);

void ghl_wenoz_reconstruction_right_left_faces(
      const double U[5],
      double *restrict Ur,
      double *restrict Ul);

#ifdef __cplusplus
}
#endif

#endif // GHL_RECONSTRUCTION_H_
