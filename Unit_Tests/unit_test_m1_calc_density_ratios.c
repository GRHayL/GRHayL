#include "ghl_unit_tests.h"
#include "ghl.h"
#include "ghl_radiation.h"

int main(){

  double mb = 1;
  double rho = .5;
  double eps = 1;
  ghl_metric_quantities metric;
  double n_e = 1;
  double n_ae = 2;
  double n_x = 3;
  double J_e = 2;
  double J_ae = 3;
  double J_x = 4;
  double nr_e = 0;
  double nr_ae = 0;
  double nr_x = 0;
  double er_e = 0;
  double er_ae = 0;
  double er_x = 0;

  ghl_initialize_metric(1,0,0,0,
                          1,0,0,
                            1,0,
                              1,
                        &metric);

  calc_density_ratios(mb, rho, eps, &metric, n_e, n_ae, n_x, J_e, J_ae, J_x, &nr_e, &nr_ae, &nr_x, &er_e, &er_ae, &er_x);

  printf("Number Density ratio (Electrons neu.):      %e\n", nr_e);
  printf("Number Density ratio (Anti-Electrons neu.): %e\n", nr_ae);
  printf("Number Density ratio (Heavy neu.):          %e\n\n", nr_x);
  printf("Energy Density ratio (Electrons neu.):      %e\n", er_e);
  printf("Energy Density ratio (Anti-Electrons neu.): %e\n", er_ae);
  printf("Energy Density ratio (Heavy neu.):          %e\n\n", er_x);

  return 1;
}
