#include <math.h>
#include <stdlib.h>

#include "ghl_unit_tests.h"

typedef void (*ghl_reconstruction_test_function)(
      const double *restrict U,
      double *restrict Ur,
      double *restrict Ul);

static void ghl_weno7z_reconstruction_wrapper(
      const double *restrict U,
      double *restrict Ur,
      double *restrict Ul) {
  ghl_weno7z_reconstruction(U, Ur, Ul);
}

static void ghl_weno9z_reconstruction_wrapper(
      const double *restrict U,
      double *restrict Ur,
      double *restrict Ul) {
  ghl_weno9z_reconstruction(U, Ur, Ul);
}

static double sin_cell_average(
      const double x_left,
      const double dx) {
  return (cos(x_left) - cos(x_left + dx)) / dx;
}

static void compute_interface_errors(
      ghl_reconstruction_test_function reconstruct,
      const int nghosts,
      const int arraylength,
      double *restrict error_r,
      double *restrict error_l) {
  const double xmin = -1.0;
  const double xmax = 1.0;
  const double dx = (xmax - xmin) / arraylength;

  double *cell_avg = (double *)malloc(sizeof(double) * arraylength);
  if(cell_avg == NULL)
    ghl_error("Could not allocate cell-average array for convergence test.\n");

  for(int index = 0; index < arraylength; ++index) {
    const double x_left = xmin + index * dx;
    cell_avg[index] = sin_cell_average(x_left, dx);
  }

  *error_r = 0.0;
  *error_l = 0.0;

  int sample_count = 0;
  for(int index = nghosts; index < arraylength - (nghosts - 1); ++index) {
    double var_r, var_l;
    reconstruct(&cell_avg[index - nghosts], &var_r, &var_l);

    const double x_face = xmin + index * dx;
    const double exact = sin(x_face);

    *error_r += fabs(var_r - exact);
    *error_l += fabs(var_l - exact);
    ++sample_count;
  }

  *error_r /= sample_count;
  *error_l /= sample_count;
  free(cell_avg);
}

static void test_scheme_convergence(
      const char *restrict name,
      ghl_reconstruction_test_function reconstruct,
      const int nghosts,
      const int *restrict resolutions,
      const int resolution_count,
      const double min_rate) {
  double error_r[4];
  double error_l[4];

  for(int i = 0; i < resolution_count; ++i) {
    compute_interface_errors(
        reconstruct,
        nghosts,
        resolutions[i],
        &error_r[i],
        &error_l[i]);

    ghl_info("%s convergence N=%d : L1(Ur)=%.16e L1(Ul)=%.16e\n",
             name, resolutions[i], error_r[i], error_l[i]);
  }

  for(int i = 1; i < resolution_count; ++i) {
    const double rate_r = log(error_r[i - 1] / error_r[i]) / log(2.0);
    const double rate_l = log(error_l[i - 1] / error_l[i]) / log(2.0);

    ghl_info("%s observed order %d->%d : Ur=%.6f Ul=%.6f\n",
             name, resolutions[i - 1], resolutions[i], rate_r, rate_l);

    if(rate_r < min_rate || rate_l < min_rate) {
      ghl_error("%s convergence test failed.\n"
                "  Observed orders from N=%d to N=%d were Ur=%.6f and Ul=%.6f.\n"
                "  Minimum accepted order is %.6f.\n",
                name, resolutions[i - 1], resolutions[i], rate_r, rate_l, min_rate);
    }
  }
}

int main(int argc, char **argv) {
  (void)argc;
  (void)argv;

  const int weno7z_resolutions[] = {20, 40, 80};
  const int weno9z_resolutions[] = {10, 20, 40};

  test_scheme_convergence(
      "WENO7Z",
      ghl_weno7z_reconstruction_wrapper,
      4,
      weno7z_resolutions,
      3,
      6.8);

  test_scheme_convergence(
      "WENO9Z",
      ghl_weno9z_reconstruction_wrapper,
      5,
      weno9z_resolutions,
      3,
      8.5);

  ghl_info("WENO7Z/WENO9Z convergence test has passed!\n");
  return 0;
}
