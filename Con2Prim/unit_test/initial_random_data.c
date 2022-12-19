#include "stdlib.h"
#include "../con2prim_gem.h"

inline double randf(double low,double high) {
    return (rand()/(double)(RAND_MAX))*(high-low)+low;
}

void initial_random_data( const double xrho, const double xpress,
                          const bool random_metric,
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

  // For the metric, we must ensure that the metric is somewhat reasonable (i.e
  // positive-definite) and therefore solve for the final component.
  double phi, gxx, gyy, gzz, gxy, gxz, gyz, lapse, betax, betay, betaz;
  if(random_metric) {
    gyy = 1.0 + randf(0.0,1.0e-1);
    gzz = 1.0 + randf(0.0,1.0e-1);
    gxy = randf(-1.0e-1,1.0e-1);
    gxz = randf(-1.0e-1,1.0e-1);
    gyz = randf(-1.0e-1,1.0e-1);
    phi = randf(0.0,2.0);
    gxx = ( 12.0*phi - gxy * (gyz*gxz - gxy*gzz) - gxz * (gxy*gyz - gyy*gxz) )/(gyy*gzz - gyz*gyz);
    lapse = 1.0; //randf(1.0e-10,1.0);
    betax = 0.0; //v*randf(-1.0,1.0);
    betay = 0.0; //sqrt(v*v - betax*betax)*randf(0.0,1.0);
    betaz = 0.0; //sqrt(v*v - betax*betax - betay*betay);
  } else {
    gxx = 1.0;
    gyy = 1.0;
    gzz = 1.0;
    gxy = 0.0;
    gxz = 0.0;
    gyz = 0.0;
    lapse = 1.0;
    betax = 0.0;
    betay = 0.0;
    betaz = 0.0;
  }

  double poison = 1e200;
  // Store the metric randomized values into the structs
  initialize_metric(lapse,
                    gxx, gxy, gxz,
                    gyy, gyz, gzz,
                    betax, betay, betaz,
                    metric);

  initialize_primitives(xrho, xpress, poison,
                        vx, vy, vz, Bx, By, Bz,
                        poison, poison, poison, // entropy, Y_e=xye, temp=xtemp
                        prims);
}
