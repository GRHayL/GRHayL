#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  const int arraylength = 100000;
  const int NGHOSTS     = 3;

  double *var = (double*) malloc(sizeof(double)*arraylength);
  double *var_r = (double*) malloc(sizeof(double)*arraylength);
  double *var_l = (double*) malloc(sizeof(double)*arraylength);

  // Initialize random data.
  for(int index=0; index<arraylength/2; index++) {
    var[index] = randf(-1.0,1.0);
  }
  for(int index=arraylength/2; index<arraylength; index++) {
    var[index] = randf(1e-14, 1e-1);
  }

  FILE* outfile = fopen_with_check("WENOZ_reconstruction_input.bin", "wb");
  fwrite(&arraylength, sizeof(int), 1, outfile);
  fwrite(var, sizeof(double), arraylength, outfile);
  fclose(outfile);

  double *varpert = (double*) malloc(sizeof(double)*arraylength);
  for(int index=0; index<arraylength; index++) {
    varpert[index] = var[index] * (1 + randf(-1.0,1.0)*1.0e-14);
  }

  outfile = fopen_with_check("WENOZ_reconstruction_output.bin", "wb");
  FILE *outpert = fopen("WENOZ_reconstruction_output_pert.bin", "wb");

  for(int index=NGHOSTS; index<arraylength-2; index++) {
    ghl_wenoz_reconstruction(&var[index-3], &var_r[index], &var_l[index]);
  }

  fwrite(var_r, sizeof(double), arraylength, outfile);
  fwrite(var_l, sizeof(double), arraylength, outfile);

  for(int index=NGHOSTS; index<arraylength-2; index++) {
    ghl_wenoz_reconstruction(&varpert[index-3], &var_r[index], &var_l[index]);
  }
  
  fwrite(var_r, sizeof(double), arraylength, outpert);
  fwrite(var_l, sizeof(double), arraylength, outpert);

  fclose(outfile);
  fclose(outpert);
}
