#include "ghl.h"
#include "ghl_radiation.h"

void calc_density_ratios(const double mb,
                         const double rho,
                         const double eps,
                         const ghl_metric_quantities *metric,
                         const double n_e,
                         const double n_ae,
                         const double n_x,
                         const double J_e,
                         const double J_ae,
                         const double J_x,
                         double *number_ratio_e,
                         double *number_ratio_ae,
                         double *number_ratio_x,
                         double *energy_ratio_e,
                         double *energy_ratio_ae,
                         double *energy_ratio_x){
  
  const double nb = rho/mb;
  *number_ratio_e = n_e/nb/metric->sqrt_detgamma;
  *number_ratio_ae = n_ae/nb/metric->sqrt_detgamma;
  *number_ratio_x = n_x/nb/metric->sqrt_detgamma;

  const double rho_e_gas = rho*(1+eps);
  const double rho_e_e = J_e/metric->sqrt_detgamma;
  const double rho_e_ae = J_ae/metric->sqrt_detgamma;
  const double rho_e_x = J_x/metric->sqrt_detgamma;
  const double rho_e_tot = rho_e_gas + rho_e_e + rho_e_ae + rho_e_x;
  *energy_ratio_e = rho_e_e/rho_e_tot;
  *energy_ratio_ae = rho_e_ae/rho_e_tot;
  *energy_ratio_x = rho_e_x/rho_e_tot;

}
