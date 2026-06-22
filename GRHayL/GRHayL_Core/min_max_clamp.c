#include <math.h>

int ghl_imin(int a, int b) {
  return a < b ? a : b;
}

int ghl_imax(int a, int b) {
  return a > b ? a : b;
}

int ghl_iclamp(int x, int x_min, int x_max) {
  return ghl_imin(ghl_imax(x, x_min), x_max);
}

double ghl_clamp(double x, double x_min, double x_max) {
  return fmin(fmax(x, x_min), x_max);
}
