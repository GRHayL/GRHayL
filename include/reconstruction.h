#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

#include "GRHayL.h"

void simple_ppm(const double rho[6], const double pressure[6], const double var_data[][6],
                const int num_vars, const double v_flux_dirn[6],
                const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
                double *restrict rhor, double *restrict rhol,
                double *restrict pressr, double *restrict pressl,
                double *restrict var_datar, double *restrict var_datal);

void simple_ppm_no_rho_P(const double pressure[6], const double var_data[][6],
                const int num_vars, const double v_flux_dirn[6],
                const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
                double *restrict var_datar, double *restrict var_datal);

#endif // RECONSTRUCTION_H_
