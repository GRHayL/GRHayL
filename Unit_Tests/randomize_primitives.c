#include "stdlib.h"
#include "ghl_unit_tests.h"

void ghl_randomize_primitives(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double press,
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
