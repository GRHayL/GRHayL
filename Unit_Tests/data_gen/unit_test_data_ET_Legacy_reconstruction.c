#include "unit_tests.h"

int main(int argc, char **argv) {
  const double poison = 1e200;

  const double W_max = 10.0;
  const int neos = 1;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;
  const double Gamma_th = 2.0;

  eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
                                             poison, poison, poison,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  const int dirlength = 20;
  const int arraylength = dirlength*dirlength*dirlength;

  double *rho = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);

  // Initialize random data. Note that for this test,
  // we needn't worry too much with physical reasonableness.
  // Discontinuities or v>1 might give physically strange
  // output, but as long as cmin+cmax != 0 the function will
  // produce numbers.
  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        rho[index] = randf(0.0,100.0);

        double P_cold;
        ghl_hybrid_compute_P_cold(&eos, rho[index], &P_cold);
        press[index] = randf(0.5*P_cold, 100*P_cold);

        vx[index] = randf(-1.0,1.0);
        vy[index] = randf(-1.0,1.0);
        vz[index] = randf(-1.0,1.0);
  }

  FILE* outfile = fopen_with_check("ET_Legacy_reconstruction_input.bin", "wb");
  fwrite(&dirlength, sizeof(int), 1, outfile);
  fwrite(rho  , sizeof(double), arraylength, outfile);
  fwrite(press, sizeof(double), arraylength, outfile);
  fwrite(vx   , sizeof(double), arraylength, outfile);
  fwrite(vy   , sizeof(double), arraylength, outfile);
  fwrite(vz   , sizeof(double), arraylength, outfile);
  fclose(outfile);

  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        rho[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        press[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        vx[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        vy[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        vz[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
  }

  outfile = fopen_with_check("ET_Legacy_reconstruction_input_pert.bin", "wb");
  fwrite(rho  , sizeof(double), arraylength, outfile);
  fwrite(press, sizeof(double), arraylength, outfile);
  fwrite(vx   , sizeof(double), arraylength, outfile);
  fwrite(vy   , sizeof(double), arraylength, outfile);
  fwrite(vz   , sizeof(double), arraylength, outfile);
  fclose(outfile);
}
