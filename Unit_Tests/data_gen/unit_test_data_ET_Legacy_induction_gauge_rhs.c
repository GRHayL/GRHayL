#include "unit_tests.h"

int main(int argc, char **argv) {


  const int dirlength     = 21;
  const int arraylength = dirlength*dirlength*dirlength;

  const int gauss_center = 10;
  const double dX[3] = {0.1, 0.1, 0.1};

  double *gupxx = (double*) malloc(sizeof(double)*arraylength);
  double *gupxy = (double*) malloc(sizeof(double)*arraylength);
  double *gupxz = (double*) malloc(sizeof(double)*arraylength);
  double *gupyy = (double*) malloc(sizeof(double)*arraylength);
  double *gupyz = (double*) malloc(sizeof(double)*arraylength);
  double *gupzz = (double*) malloc(sizeof(double)*arraylength);

  double *psi = (double*) malloc(sizeof(double)*arraylength);
  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde = (double*) malloc(sizeof(double)*arraylength);
  double *Ax = (double*) malloc(sizeof(double)*arraylength);
  double *Ay = (double*) malloc(sizeof(double)*arraylength);
  double *Az = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betaz_interp = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_Phi_minus_betaj_A_j_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Ax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Ay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Az_interp = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Az_rhs = (double*) malloc(sizeof(double)*arraylength);

  // Data which should be written before it is used is poisoned to validate behavior.
  // RHSs for A are set to 0 because they are assumed to already contain the zero-gauge
  // contribution to the RHS before entering this function.
  const double poison = 1e200;
  // This cannot be in parallel because the randomized quantities need to happen in the
  // same order every time.
  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        metric_quantities ADM_metric;
        double agxx, agxy, agxz, agyy, agyz, agzz;
        ghl_randomize_metric(
              &lapse[index], &betax[index], &betay[index], &betaz[index],
              &agxx, &agxy, &agxz,
              &agyy, &agyz, &agzz);

        // For most tests, we just set the shift to zero to simplify testing;
        // however, this code has branches which explicitly depend on the shift,
        // so we need to actually make it non-zero.

        betax[index] = randf(-1e-3,1e-3);
        betay[index] = randf(-1e-3,1e-3);
        betaz[index] = randf(-1e-3,1e-3);

        ghl_initialize_metric(
              lapse[index], betax[index], betay[index], betaz[index],
              agxx, agxy, agxz, agyy, agyz, agzz,
              &ADM_metric);

        ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        psi[index]   = sqrt(metric_aux.psi2);
        gupxx[index] = metric_aux.psi4*ADM_metric.gammaUU[0][0];
        gupxy[index] = metric_aux.psi4*ADM_metric.gammaUU[0][1];
        gupxz[index] = metric_aux.psi4*ADM_metric.gammaUU[0][2];
        gupyy[index] = metric_aux.psi4*ADM_metric.gammaUU[1][1];
        gupyz[index] = metric_aux.psi4*ADM_metric.gammaUU[1][2];
        gupzz[index] = metric_aux.psi4*ADM_metric.gammaUU[2][2];

        const int x = abs(i-gauss_center)*dX[0];
        const int y = abs(j-gauss_center)*dX[1];
        const int z = abs(k-gauss_center)*dX[2];
        const double r2 = x*x + y*y + z*z;

        phitilde[index] = exp(-r2/(2.0*1.0));;
        Ax[index]       = exp(-r2/(2.0*4.0));
        Ay[index]       = exp(-r2/(2.0*5.0));;
        Az[index]       = exp(-r2/(2.0*2.0));;

        alpha_interp[index]  = poison;
        betax_interp[index] = poison;
        betay_interp[index] = poison;
        betaz_interp[index] = poison;

        alpha_Phi_minus_betaj_A_j_interp[index] = poison;
        alpha_sqrtg_Ax_interp[index]            = poison;
        alpha_sqrtg_Ay_interp[index]            = poison;
        alpha_sqrtg_Az_interp[index]            = poison;

        Ax_rhs[index]       = 0.0;
        Ay_rhs[index]       = 0.0;
        Az_rhs[index]       = 0.0;
        phitilde_rhs[index] = poison;
  }

  FILE* outfile = fopen("ET_Legacy_induction_gauge_rhs_input.bin", "wb");
  check_file_was_successfully_open(outfile, "unit_test_ET_Legacy_induction_gauge_rhs_input.bin");

  fwrite(&dirlength, sizeof(int), 1, outfile);
  fwrite(gupxx, sizeof(double), arraylength, outfile);
  fwrite(gupxy, sizeof(double), arraylength, outfile);
  fwrite(gupxz, sizeof(double), arraylength, outfile);
  fwrite(gupyy, sizeof(double), arraylength, outfile);
  fwrite(gupyz, sizeof(double), arraylength, outfile);
  fwrite(gupzz, sizeof(double), arraylength, outfile);

  fwrite(psi,   sizeof(double), arraylength, outfile);
  fwrite(lapse, sizeof(double), arraylength, outfile);
  fwrite(betax, sizeof(double), arraylength, outfile);
  fwrite(betay, sizeof(double), arraylength, outfile);
  fwrite(betaz, sizeof(double), arraylength, outfile);

  fwrite(phitilde, sizeof(double), arraylength, outfile);
  fwrite(Ax,       sizeof(double), arraylength, outfile);
  fwrite(Ay,       sizeof(double), arraylength, outfile);
  fwrite(Az,       sizeof(double), arraylength, outfile);

  fclose(outfile);

  for(int index=0; index<arraylength; index++) {
    
    gupxx[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupxy[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupxz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupyy[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupyz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupzz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    lapse[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    psi[index]      *= (1.0 + randf(-1,1)*1.0e-14);
    betax[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    betay[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    betaz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    phitilde[index] *= (1.0 + randf(-1,1)*1.0e-14);
    Ax[index]       *= (1.0 + randf(-1,1)*1.0e-14);
    Ay[index]       *= (1.0 + randf(-1,1)*1.0e-14);
    Az[index]       *= (1.0 + randf(-1,1)*1.0e-14);
  }

  outfile = fopen("ET_Legacy_induction_gauge_rhs_input_pert.bin", "wb");
  check_file_was_successfully_open(outfile, "unit_test_ET_Legacy_induction_gauge_rhs_input_pert.bin");

  fwrite(gupxx, sizeof(double), arraylength, outfile);
  fwrite(gupxy, sizeof(double), arraylength, outfile);
  fwrite(gupxz, sizeof(double), arraylength, outfile);
  fwrite(gupyy, sizeof(double), arraylength, outfile);
  fwrite(gupyz, sizeof(double), arraylength, outfile);
  fwrite(gupzz, sizeof(double), arraylength, outfile);

  fwrite(psi,   sizeof(double), arraylength, outfile);
  fwrite(lapse, sizeof(double), arraylength, outfile);
  fwrite(betax, sizeof(double), arraylength, outfile);
  fwrite(betay, sizeof(double), arraylength, outfile);
  fwrite(betaz, sizeof(double), arraylength, outfile);

  fwrite(phitilde, sizeof(double), arraylength, outfile);
  fwrite(Ax,       sizeof(double), arraylength, outfile);
  fwrite(Ay,       sizeof(double), arraylength, outfile);
  fwrite(Az,       sizeof(double), arraylength, outfile);

  fclose(outfile);
}
