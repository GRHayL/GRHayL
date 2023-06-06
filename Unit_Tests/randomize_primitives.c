#include "stdlib.h"
#include "unit_tests.h"

void randomize_primitives(
      const eos_parameters *restrict eos,
      const double rho,
      const double press,
      double *restrict eps,
      double *restrict vx,
      double *restrict vy,
      double *restrict vz,
      double *restrict Bx,
      double *restrict By,
      double *restrict Bz) {

  // Fix Y_e
  //const double Ye_test = 0.1;

  // Fix W
  const double W_test  = 2.0;

  // Fix log10(Pmag/P)
  const double logPmoP = -5.0;

  // First, set the velocities
  // Velocity magnitude
  const double v = sqrt(1.0-1.0/(W_test*W_test));
  *vx = v*randf(0.0,1.0);
  const double vx2 = (*vx)*(*vx);
  *vy = sqrt(v*v - vx2)*randf(0.0,1.0);
  const double vy2 = (*vy)*(*vy);
  *vz = sqrt(v*v - vx2 - vy2);

  // Set eps
  double P_cold = 0.0;
  double eps_cold = 0.0;
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, rho, &P_cold, &eps_cold);
  *eps = eps_cold + (press-P_cold)/(eos->Gamma_th-1.0)/rho;

  // Now the magnetic fields. We'll set them aligned
  // with the velocities, for simplicity.
  const double Bhatx = *vx/v;
  const double Bhaty = *vy/v;
  const double Bhatz = *vz/v;
  const double B     = sqrt(2.0*pow(10.0,logPmoP)*press);
  *Bx    = -Bhatx * B;
  *By    = -Bhaty * B;
  *Bz    = -Bhatz * B;

}
