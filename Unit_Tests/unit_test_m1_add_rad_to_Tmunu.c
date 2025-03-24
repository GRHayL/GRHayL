#include "ghl_unit_tests.h"
#include "ghl_radiation.h"

int main(){

  ghl_stress_energy Tmunu_test;
  ghl_initialize_stress_energy(1, 0, 0, 0,
                                  1, 0, 0,
                                     1, 0,
                                        1,
                               &Tmunu_test);

  ghl_stress_energy Tmunu_rad_test;
  ghl_initialize_stress_energy(1, 0, 0, 0,
                                  1, 0, 0,
                                     1, 0,
                                        1,
                                &Tmunu_rad_test);

  ghl_metric_quantities metric;
  ghl_initialize_metric(1, 0, 0, 0,
                           1, 0, 0,
                              1, 0,
                                 1,
                            &metric);

  add_rad_to_Tmunu(&Tmunu_test, &Tmunu_rad_test, &metric);

  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 4; j++){
      printf("Tmunu: index1 = %d, index2 = %d, Value = %e\n", i, j, Tmunu_test.T4[i][j]);
    }
  }



}
