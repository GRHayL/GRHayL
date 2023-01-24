#include "stdlib.h"
#include "unit_tests.h"

void initial_random_data( const double xrho, const double xpress,
                          metric_quantities *restrict metric,
                          primitive_quantities *restrict prims ) {

  // Fix Y_e
  //const double Ye_test = 0.1;

  // Fix W
  const double W_test  = 2.0;

  // Fix log10(Pmag/P)
  const double logPmoP = -5.0;

  // First, set the velocities
  // Velocity magnitude
  const double v = sqrt(1.0-1.0/(W_test*W_test));
  const double vx = v*randf(0.0,1.0);
  const double vy = sqrt(v*v - vx*vx)*randf(0.0,1.0);
  const double vz = sqrt(v*v - vx*vx - vy*vy);

  // Now the magnetic fields. We'll set them aligned
  // with the velocities, for simplicity.
  const double Bhatx = vx/v;
  const double Bhaty = vy/v;
  const double Bhatz = vz/v;
  const double B     = sqrt(2.0*pow(10.0,logPmoP)*xpress);
  const double Bx    = -Bhatx * B;
  const double By    = -Bhaty * B;
  const double Bz    = -Bhatz * B;

  randomize_metric(metric);

  double poison = 1e200;
  initialize_primitives(xrho, xpress, poison,
                        vx, vy, vz, Bx, By, Bz,
                        poison, poison, poison, // entropy, Y_e=xye, temp=xtemp
                        prims);
}
