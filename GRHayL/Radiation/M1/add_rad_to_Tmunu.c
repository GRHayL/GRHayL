#include "ghl.h"
#include "ghl_radiation.h"

void add_rad_to_Tmunu(ghl_stress_energy *T4DD,
                      const ghl_stress_energy *rT4DD,
                      const ghl_metric_quantities *metric){

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      T4DD->T4[i][j] += rT4DD->T4[i][j]/metric->sqrt_detgamma;
    }
  }

}
