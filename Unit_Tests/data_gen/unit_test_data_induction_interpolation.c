#include "unit_tests.h"

int main(int argc, char **argv) {
  const int dirlength     = 21;
  const int arraylength = dirlength*dirlength*dirlength;

  const int gauss_center = 10;
  const double dX[3] = {0.1, 0.1, 0.1};

  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *psi    = (double*) malloc(sizeof(double)*arraylength);
  double *gtupxx = (double*) malloc(sizeof(double)*arraylength);
  double *gtupxy = (double*) malloc(sizeof(double)*arraylength);
  double *gtupxz = (double*) malloc(sizeof(double)*arraylength);
  double *gtupyy = (double*) malloc(sizeof(double)*arraylength);
  double *gtupyz = (double*) malloc(sizeof(double)*arraylength);
  double *gtupzz = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde = (double*) malloc(sizeof(double)*arraylength);
  double *Ax = (double*) malloc(sizeof(double)*arraylength);
  double *Ay = (double*) malloc(sizeof(double)*arraylength);
  double *Az = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betaz_interp = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_Phi_minus_betaj_A_j_interp = (double*) malloc(sizeof(double)*arraylength);
  double *sqrtg_Ax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *sqrtg_Ay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *sqrtg_Az_interp = (double*) malloc(sizeof(double)*arraylength);

  // Data which should be written before it is used is poisoned to validate behavior.
  // RHSs for A are set to 0 because they are assumed to already contain the zero-gauge
  // contribution to the RHS before entering this function.
  const double poison = 1e200;
  // This cannot be in parallel because the randomized quantities need to happen in the
  // same order every time.
  for(int k=0; k<dirlength; k++) {
    for(int j=0; j<dirlength; j++) {
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        ghl_metric_quantities ADM_metric;
        ghl_randomize_metric(
              &lapse[index], &betax[index], &betay[index], &betaz[index],
              &gxx[index], &gxy[index], &gxz[index],
              &gyy[index], &gyz[index], &gzz[index]);

        // For most tests, we just set the shift to zero to simplify testing;
        // however, this code has branches which explicitly depend on the shift,
        // so we need to actually make it non-zero.

        betax[index] = randf(-1e-3,1e-3);
        betay[index] = randf(-1e-3,1e-3);
        betaz[index] = randf(-1e-3,1e-3);

        ghl_initialize_metric(
              lapse[index], betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        psi[index]    = sqrt(metric_aux.psi2);
        gtupxx[index] = metric_aux.psi4*ADM_metric.gammaUU[0][0];
        gtupxy[index] = metric_aux.psi4*ADM_metric.gammaUU[0][1];
        gtupxz[index] = metric_aux.psi4*ADM_metric.gammaUU[0][2];
        gtupyy[index] = metric_aux.psi4*ADM_metric.gammaUU[1][1];
        gtupyz[index] = metric_aux.psi4*ADM_metric.gammaUU[1][2];
        gtupzz[index] = metric_aux.psi4*ADM_metric.gammaUU[2][2];

        const int x = abs(i-gauss_center)*dX[0];
        const int y = abs(j-gauss_center)*dX[1];
        const int z = abs(k-gauss_center)*dX[2];
        const double r2 = x*x + y*y + z*z;

        phitilde[index] = exp(-r2/(2.0*1.0));;
        Ax[index]       = exp(-r2/(2.0*4.0));
        Ay[index]       = exp(-r2/(2.0*5.0));;
        Az[index]       = exp(-r2/(2.0*2.0));;

        alpha_interp[index] = poison;
        betax_interp[index] = poison;
        betay_interp[index] = poison;
        betaz_interp[index] = poison;

        alpha_Phi_minus_betaj_A_j_interp[index] = poison;
        sqrtg_Ax_interp[index] = poison;
        sqrtg_Ay_interp[index] = poison;
        sqrtg_Az_interp[index] = poison;
      }
    }
  }

  // First, I write out the variables that are needed for all the codes.
  FILE* outfile = fopen_with_check("induction_interpolation_input.bin", "wb");
  fwrite(&dirlength, sizeof(int), 1, outfile);
  fwrite(lapse, sizeof(double), arraylength, outfile);
  fwrite(betax, sizeof(double), arraylength, outfile);
  fwrite(betay, sizeof(double), arraylength, outfile);
  fwrite(betaz, sizeof(double), arraylength, outfile);

  fwrite(phitilde, sizeof(double), arraylength, outfile);
  fwrite(Ax,       sizeof(double), arraylength, outfile);
  fwrite(Ay,       sizeof(double), arraylength, outfile);
  fwrite(Az,       sizeof(double), arraylength, outfile);
  fclose(outfile);

  outfile = fopen_with_check("induction_interpolation_ADM_input.bin", "wb");
  fwrite(gxx, sizeof(double), arraylength, outfile);
  fwrite(gxy, sizeof(double), arraylength, outfile);
  fwrite(gxz, sizeof(double), arraylength, outfile);
  fwrite(gyy, sizeof(double), arraylength, outfile);
  fwrite(gyz, sizeof(double), arraylength, outfile);
  fwrite(gzz, sizeof(double), arraylength, outfile);
  fclose(outfile);

  outfile = fopen_with_check("induction_interpolation_BSSN_input.bin", "wb");
  fwrite(psi,    sizeof(double), arraylength, outfile);
  fwrite(gtupxx, sizeof(double), arraylength, outfile);
  fwrite(gtupxy, sizeof(double), arraylength, outfile);
  fwrite(gtupxz, sizeof(double), arraylength, outfile);
  fwrite(gtupyy, sizeof(double), arraylength, outfile);
  fwrite(gtupyz, sizeof(double), arraylength, outfile);
  fwrite(gtupzz, sizeof(double), arraylength, outfile);
  fclose(outfile);

  ghl_test_compute_ccc_ADM(
        dirlength, lapse, betax, betay, betaz,
        gxx, gxy, gxz, gyy, gyz, gzz,
        phitilde, Ax, Ay, Az,
        alpha_interp,
        betax_interp,
        betay_interp,
        betaz_interp,
        alpha_Phi_minus_betaj_A_j_interp,
        sqrtg_Ax_interp,
        sqrtg_Ay_interp,
        sqrtg_Az_interp);

  outfile = fopen_with_check("induction_interpolation_ccc_ADM_output.bin", "wb");
  fwrite(alpha_interp, sizeof(double), arraylength, outfile);
  fwrite(betax_interp, sizeof(double), arraylength, outfile);
  fwrite(betay_interp, sizeof(double), arraylength, outfile);
  fwrite(betaz_interp, sizeof(double), arraylength, outfile);
  fwrite(alpha_Phi_minus_betaj_A_j_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ax_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ay_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Az_interp, sizeof(double), arraylength, outfile);
  fclose(outfile);

  ghl_test_compute_ccc_BSSN(
        dirlength, lapse, betax, betay, betaz, psi,
        gtupxx, gtupxy, gtupxz, gtupyy, gtupyz, gtupzz,
        phitilde, Ax, Ay, Az,
        alpha_interp,
        betax_interp,
        betay_interp,
        betaz_interp,
        alpha_Phi_minus_betaj_A_j_interp,
        sqrtg_Ax_interp,
        sqrtg_Ay_interp,
        sqrtg_Az_interp);

  outfile = fopen_with_check("induction_interpolation_ccc_BSSN_output.bin", "wb");
  fwrite(alpha_interp, sizeof(double), arraylength, outfile);
  fwrite(betax_interp, sizeof(double), arraylength, outfile);
  fwrite(betay_interp, sizeof(double), arraylength, outfile);
  fwrite(betaz_interp, sizeof(double), arraylength, outfile);
  fwrite(alpha_Phi_minus_betaj_A_j_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ax_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ay_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Az_interp, sizeof(double), arraylength, outfile);
  fclose(outfile);

  // For the vertex-centered tests, we can just use the same data, but
  // they'll represent different coordinate locations relative to the
  // EM variables.
  ghl_test_compute_vvv_ADM(
        dirlength, lapse, betax, betay, betaz,
        gxx, gxy, gxz, gyy, gyz, gzz,
        phitilde, Ax, Ay, Az,
        alpha_Phi_minus_betaj_A_j_interp,
        sqrtg_Ax_interp,
        sqrtg_Ay_interp,
        sqrtg_Az_interp);

  outfile = fopen_with_check("induction_interpolation_vvv_ADM_output.bin", "wb");
  fwrite(alpha_Phi_minus_betaj_A_j_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ax_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ay_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Az_interp, sizeof(double), arraylength, outfile);
  fclose(outfile);

  // Perturb inputs
  for(int index=0; index<arraylength; index++) {
    lapse[index] *= (1.0 + randf(-1,1)*1.0e-14);
    betax[index] *= (1.0 + randf(-1,1)*1.0e-14);
    betay[index] *= (1.0 + randf(-1,1)*1.0e-14);
    betaz[index] *= (1.0 + randf(-1,1)*1.0e-14);

    gxx[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gxy[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gxz[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gyy[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gyz[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gzz[index] *= (1.0 + randf(-1,1)*1.0e-14);

    gtupxx[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gtupxy[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gtupxz[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gtupyy[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gtupyz[index] *= (1.0 + randf(-1,1)*1.0e-14);
    gtupzz[index] *= (1.0 + randf(-1,1)*1.0e-14);
    psi[index]    *= (1.0 + randf(-1,1)*1.0e-14);

    phitilde[index] *= (1.0 + randf(-1,1)*1.0e-14);
    Ax[index]       *= (1.0 + randf(-1,1)*1.0e-14);
    Ay[index]       *= (1.0 + randf(-1,1)*1.0e-14);
    Az[index]       *= (1.0 + randf(-1,1)*1.0e-14);
  }

  ghl_test_compute_ccc_ADM(
        dirlength, lapse, betax, betay, betaz,
        gxx, gxy, gxz, gyy, gyz, gzz,
        phitilde, Ax, Ay, Az,
        alpha_interp,
        betax_interp,
        betay_interp,
        betaz_interp,
        alpha_Phi_minus_betaj_A_j_interp,
        sqrtg_Ax_interp,
        sqrtg_Ay_interp,
        sqrtg_Az_interp);

  outfile = fopen_with_check("induction_interpolation_ccc_ADM_output_pert.bin", "wb");
  fwrite(alpha_interp, sizeof(double), arraylength, outfile);
  fwrite(betax_interp, sizeof(double), arraylength, outfile);
  fwrite(betay_interp, sizeof(double), arraylength, outfile);
  fwrite(betaz_interp, sizeof(double), arraylength, outfile);
  fwrite(alpha_Phi_minus_betaj_A_j_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ax_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ay_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Az_interp, sizeof(double), arraylength, outfile);
  fclose(outfile);

  ghl_test_compute_ccc_BSSN(
        dirlength, lapse, betax, betay, betaz, psi,
        gtupxx, gtupxy, gtupxz, gtupyy, gtupyz, gtupzz,
        phitilde, Ax, Ay, Az,
        alpha_interp,
        betax_interp,
        betay_interp,
        betaz_interp,
        alpha_Phi_minus_betaj_A_j_interp,
        sqrtg_Ax_interp,
        sqrtg_Ay_interp,
        sqrtg_Az_interp);

  outfile = fopen_with_check("induction_interpolation_ccc_BSSN_output_pert.bin", "wb");
  fwrite(alpha_interp, sizeof(double), arraylength, outfile);
  fwrite(betax_interp, sizeof(double), arraylength, outfile);
  fwrite(betay_interp, sizeof(double), arraylength, outfile);
  fwrite(betaz_interp, sizeof(double), arraylength, outfile);
  fwrite(alpha_Phi_minus_betaj_A_j_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ax_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ay_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Az_interp, sizeof(double), arraylength, outfile);
  fclose(outfile);

  ghl_test_compute_vvv_ADM(
        dirlength, lapse, betax, betay, betaz,
        gxx, gxy, gxz, gyy, gyz, gzz,
        phitilde, Ax, Ay, Az,
        alpha_Phi_minus_betaj_A_j_interp,
        sqrtg_Ax_interp,
        sqrtg_Ay_interp,
        sqrtg_Az_interp);

  outfile = fopen_with_check("induction_interpolation_vvv_ADM_output_pert.bin", "wb");
  fwrite(alpha_Phi_minus_betaj_A_j_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ax_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Ay_interp, sizeof(double), arraylength, outfile);
  fwrite(sqrtg_Az_interp, sizeof(double), arraylength, outfile);
  fclose(outfile);

  free(lapse); free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);

  free(psi);
  free(gtupxx); free(gtupxy); free(gtupxz);
  free(gtupyy); free(gtupyz); free(gtupzz);

  free(phitilde); free(Ax); free(Ay); free(Az);

  free(alpha_interp); free(betax_interp); free(betay_interp); free(betaz_interp);
  free(alpha_Phi_minus_betaj_A_j_interp);
  free(sqrtg_Ax_interp); free(sqrtg_Ay_interp); free(sqrtg_Az_interp);
}
