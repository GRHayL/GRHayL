#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

#include "GRHayL.h"

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )

// Integer constants to keep track of stencil.
enum reconstruction_stencil {
  MINUS2, MINUS1,
  PLUS_0, PLUS_1,
  PLUS_2
};

void simple_ppm(
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

void simple_ppm_no_rho_P(
      const double pressure[6],
      const double var_data[][6],
      const int num_vars,
      const double v_flux_dirn[6],
      const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
      double *restrict var_datar,
      double *restrict var_datal);

double slope_limit(
      const double dU,
      const double dUp1);

void compute_UrUl_onevar(const double U[5], double *restrict Ur, double *restrict Ul);

double shock_detection_ftilde(const double P[5], const double v_flux_dirn[5]);
void steepen_rhor_rhol(const double rho[5], const double P[5], const double Gamma_eff,
                              double *restrict rhor, double *restrict rhol);
void flatten_and_monotonize_Ur_and_Ul(const double U, const double ftilde, double *restrict Ur, double *restrict Ul);
#endif // RECONSTRUCTION_H_
