#include "unit_tests.h"

int main(int argc, char **argv) {

  const int arraylength = 100000;

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

  FILE* outfile = fopen_with_check("PLM_reconstruction_input.bin", "wb");
  fwrite(&arraylength, sizeof(int), 1, outfile);
  fwrite(var, sizeof(double), arraylength, outfile);
  fclose(outfile);

  double *varpert = (double*) malloc(sizeof(double)*arraylength);
  for(int index=0; index<arraylength; index++) {
    varpert[index] = var[index] * (1 + randf(-1.0,1.0)*1.0e-14);
  }

  outfile = fopen_with_check("PLM_reconstruction_output.bin", "wb");
  FILE *outpert = fopen("PLM_reconstruction_output_pert.bin", "wb");

  for(int method=0; method<3; method++) {
    void(*ghl_reconstruction)(
          const double,
          const double,
          const double,
          const double,
          double *restrict,
          double *restrict);

    switch(method) {
      case 0:
        ghl_reconstruction = &ghl_minmod_reconstruction;
        break;
      case 1:
        ghl_reconstruction = &ghl_mc_reconstruction;
        break;
      case 2:
        ghl_reconstruction = &ghl_superbee_reconstruction;
        break;
    }
    for(int index=0; index<arraylength; index++) {
      ghl_reconstruction(var[index-2], var[index-1], var[index], var[index+1], &var_r[index], &var_l[index]);
    }
    fwrite(var_r, sizeof(double), arraylength, outfile);
    fwrite(var_l, sizeof(double), arraylength, outfile);

    for(int index=0; index<arraylength; index++) {
      ghl_reconstruction(varpert[index-2], varpert[index-1], varpert[index], varpert[index+1], &var_r[index], &var_l[index]);
    }
    fwrite(var_r, sizeof(double), arraylength, outpert);
    fwrite(var_l, sizeof(double), arraylength, outpert);
  }
  fclose(outfile);
  fclose(outpert);
}
