#include "unit_tests.h"

int main(int argc, char **argv) {

  // First test: reset_prims_to_atmosphere() function
  // This test requires no input data

  // Second test: GRHayL_enforce_detgtij_and_initialize_metric() function

  // We simply generate a large number of random metrics to send through
  // the function.

  const int arraylength = 100;

  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  for(int i=0; i<arraylength; i++)
    ghl_randomize_metric(
        &lapse[i], &betax[i], &betay[i], &betaz[i],
        &gxx[i], &gxy[i], &gxz[i],
        &gyy[i], &gyz[i], &gzz[i]);

  FILE* infile = fopen("grhayL_core_test_suite_input.bin", "wb");
  check_file_was_successfully_open(infile, "grhayL_core_test_suite_input.bin");

  fwrite(&arraylength, sizeof(int), 1, infile);

  fwrite(lapse, sizeof(double), arraylength, infile);
  fwrite(betax, sizeof(double), arraylength, infile);
  fwrite(betay, sizeof(double), arraylength, infile);
  fwrite(betaz, sizeof(double), arraylength, infile);
  fwrite(gxx,   sizeof(double), arraylength, infile);
  fwrite(gxy,   sizeof(double), arraylength, infile);
  fwrite(gxz,   sizeof(double), arraylength, infile);
  fwrite(gyy,   sizeof(double), arraylength, infile);
  fwrite(gyz,   sizeof(double), arraylength, infile);
  fwrite(gzz,   sizeof(double), arraylength, infile);

  fclose(infile);
}
