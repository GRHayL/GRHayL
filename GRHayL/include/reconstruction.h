#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

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
void ghl_ppm(
      const ghl_parameters *restrict params,
      const double rho[6],
      const double pressure[6],
      const double var_data[][6],
      const int num_vars,
      const double v_flux_dirn[6],
      const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
      double *restrict rhor,
      double *restrict rhol,
      double *restrict pressr,
      double *restrict pressl,
      double *restrict var_datar,
      double *restrict var_datal);

void ghl_ppm_no_rho_P(
      const ghl_parameters *restrict params,
      const double pressure[6],
      const double var_data[][6],
      const int num_vars,
      const double v_flux_dirn[6],
      const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
      double *restrict var_datar,
      double *restrict var_datal);

double ghl_slope_limit(
      const double dU,
      const double dUp1);

void ghl_compute_UrUl_onevar(
      const double U[5],
      double *restrict Ur,
      double *restrict Ul);

double ghl_shock_detection_ftilde(
      const ghl_parameters *restrict params,
      const double P[5],
      const double v_flux_dirn[5]);

void ghl_steepen_rhor_rhol(
      const ghl_parameters *restrict params,
      const double rho[5],
      const double P[5],
      const double Gamma_eff,
      double *restrict rhor,
      double *restrict rhol);

void ghl_flatten_and_monotonize_Ur_and_Ul(
      const double U,
      const double ftilde,
      double *restrict Ur,
      double *restrict Ul);

// TVD functions
double ghl_minmod(
      const double a,
      const double b);

double ghl_maxmod(
      const double a,
      const double b);

void ghl_minmod_reconstruction(
      const double U_m2,
      const double U_m1,
      const double U,
      const double U_p1,
      double *restrict Ur,
      double *restrict Ul);

void ghl_mc_reconstruction(
      const double U_m2,
      const double U_m1,
      const double U,
      const double U_p1,
      double *restrict Ur,
      double *restrict Ul);

void ghl_superbee_reconstruction(
      const double U_m2,
      const double U_m1,
      const double U,
      const double U_p1,
      double *restrict Ur,
      double *restrict Ul);

#ifdef __cplusplus
}
#endif

#endif // RECONSTRUCTION_H_
