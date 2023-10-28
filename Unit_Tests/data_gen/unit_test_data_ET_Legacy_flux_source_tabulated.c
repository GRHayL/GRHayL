#include "unit_tests.h"

static void randomize_rho_Ye_T_P_S(
    ghl_eos_parameters *restrict eos,
    double *restrict rho,
    double *restrict Y_e,
    double *restrict temperature,
    double *restrict press,
    double *restrict entropy) {

  const double rhoL = exp(randf(log(eos->table_rho_min), log(eos->table_rho_max)));
  const double temperatureL = exp(randf(log(eos->table_T_min), log(eos->table_T_max)));
  const double Y_eL = randf(eos->table_Y_e_min, eos->table_Y_e_max);

  double pressL, epsL, entropyL;
  ghl_tabulated_compute_P_eps_S_from_T(eos, rhoL, Y_eL, temperatureL, &pressL, &epsL, &entropyL);

  *rho = rhoL;
  *Y_e = Y_eL;
  *temperature = temperatureL;
  *press = pressL;
  *entropy = entropyL;
}

int main(int argc, char **argv) {

  if(argc != 2) {
    ghl_info("Usage: %s <EOS table path>\n", argv[0]);
    exit(1);
  }

  const char *tablepath = argv[1];
  const double rho_b_atm = 1e-12;
  const double rho_b_min = -1;
  const double rho_b_max = -1;
  const double Y_e_atm = 0.5;
  const double Y_e_min = -1;
  const double Y_e_max = -1;
  const double T_atm = 1e-2;
  const double T_min = -1;
  const double T_max = -1;

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
      tablepath, rho_b_atm, rho_b_min, rho_b_max, Y_e_atm, Y_e_min, Y_e_max, T_atm, T_min, T_max,
      &eos);

  const int dirlength = 20;
  const int arraylength = dirlength * dirlength * dirlength;

  // Allocate memory for metric
  double *gxx = (double *)malloc(sizeof(double) * arraylength);
  double *gxy = (double *)malloc(sizeof(double) * arraylength);
  double *gxz = (double *)malloc(sizeof(double) * arraylength);
  double *gyy = (double *)malloc(sizeof(double) * arraylength);
  double *gyz = (double *)malloc(sizeof(double) * arraylength);
  double *gzz = (double *)malloc(sizeof(double) * arraylength);
  double *lapse = (double *)malloc(sizeof(double) * arraylength);
  double *betax = (double *)malloc(sizeof(double) * arraylength);
  double *betay = (double *)malloc(sizeof(double) * arraylength);
  double *betaz = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for extrinsic curvature
  double *kxx = (double *)malloc(sizeof(double) * arraylength);
  double *kxy = (double *)malloc(sizeof(double) * arraylength);
  double *kxz = (double *)malloc(sizeof(double) * arraylength);
  double *kyy = (double *)malloc(sizeof(double) * arraylength);
  double *kyz = (double *)malloc(sizeof(double) * arraylength);
  double *kzz = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for primitives
  double *rho = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e = (double *)malloc(sizeof(double) * arraylength);
  double *temperature = (double *)malloc(sizeof(double) * arraylength);
  double *press = (double *)malloc(sizeof(double) * arraylength);
  double *entropy = (double *)malloc(sizeof(double) * arraylength);
  double *vx = (double *)malloc(sizeof(double) * arraylength);
  double *vy = (double *)malloc(sizeof(double) * arraylength);
  double *vz = (double *)malloc(sizeof(double) * arraylength);
  double *Bx = (double *)malloc(sizeof(double) * arraylength);
  double *By = (double *)malloc(sizeof(double) * arraylength);
  double *Bz = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for right face
  double *rho_r = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e_r = (double *)malloc(sizeof(double) * arraylength);
  double *temperature_r = (double *)malloc(sizeof(double) * arraylength);
  double *press_r = (double *)malloc(sizeof(double) * arraylength);
  double *entropy_r = (double *)malloc(sizeof(double) * arraylength);
  double *vx_r = (double *)malloc(sizeof(double) * arraylength);
  double *vy_r = (double *)malloc(sizeof(double) * arraylength);
  double *vz_r = (double *)malloc(sizeof(double) * arraylength);
  double *Bx_r = (double *)malloc(sizeof(double) * arraylength);
  double *By_r = (double *)malloc(sizeof(double) * arraylength);
  double *Bz_r = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for left face
  double *rho_l = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e_l = (double *)malloc(sizeof(double) * arraylength);
  double *temperature_l = (double *)malloc(sizeof(double) * arraylength);
  double *press_l = (double *)malloc(sizeof(double) * arraylength);
  double *entropy_l = (double *)malloc(sizeof(double) * arraylength);
  double *vx_l = (double *)malloc(sizeof(double) * arraylength);
  double *vy_l = (double *)malloc(sizeof(double) * arraylength);
  double *vz_l = (double *)malloc(sizeof(double) * arraylength);
  double *Bx_l = (double *)malloc(sizeof(double) * arraylength);
  double *By_l = (double *)malloc(sizeof(double) * arraylength);
  double *Bz_l = (double *)malloc(sizeof(double) * arraylength);

  // Initialize random data. Note that for this test,
  // we needn't worry too much with physical reasonableness.
  for(int k = 0; k < dirlength; k++) {
    for(int j = 0; j < dirlength; j++) {
      for(int i = 0; i < dirlength; i++) {
        const int index = indexf(dirlength, i, j, k);

        ghl_randomize_metric(
            &lapse[index], &betax[index], &betay[index], &betaz[index], &gxx[index], &gxy[index],
            &gxz[index], &gyy[index], &gyz[index], &gzz[index]);

        gxx[index] = 1.23;  /*randf(scale);*/
        gxy[index] = 0.03;  /*randf(scale);*/
        gxz[index] = 0.07;  /*randf(scale);*/
        gyy[index] = 2.33;  /*randf(scale);*/
        gyz[index] = 0.041; /*randf(scale);*/
        gzz[index] = 3.11;  /*randf(scale);*/

        kxx[index] = randf(0, 1.0);
        kxy[index] = randf(0, 1.0);
        kxz[index] = randf(0, 1.0);
        kyy[index] = randf(0, 1.0);
        kyz[index] = randf(0, 1.0);
        kzz[index] = randf(0, 1.0);

        // clang-format off
        randomize_rho_Ye_T_P_S(&eos, &rho  [index], &Y_e  [index], &temperature  [index], &press  [index], &entropy  [index]);
        randomize_rho_Ye_T_P_S(&eos, &rho_r[index], &Y_e_r[index], &temperature_r[index], &press_r[index], &entropy_r[index]);
        randomize_rho_Ye_T_P_S(&eos, &rho_l[index], &Y_e_l[index], &temperature_l[index], &press_l[index], &entropy_l[index]);
        // clang-format on

        ghl_randomize_primitives(
            &eos, rho[index], press[index], &vx[index], &vy[index], &vz[index], &Bx[index],
            &By[index], &Bz[index]);

        ghl_randomize_primitives(
            &eos, rho_r[index], press_r[index], &vx_r[index], &vy_r[index], &vz_r[index],
            &Bx_r[index], &By_r[index], &Bz_r[index]);

        ghl_randomize_primitives(
            &eos, rho_l[index], press_l[index], &vx_l[index], &vy_l[index], &vz_l[index],
            &Bx_l[index], &By_l[index], &Bz_l[index]);
      }
    }
  }

  FILE *outfile = fopen_with_check("ET_Legacy_flux_source_input.bin", "wb");
  fwrite(&dirlength, sizeof(int), 1, outfile);
  fwrite(gxx, sizeof(double), arraylength, outfile);
  fwrite(gxy, sizeof(double), arraylength, outfile);
  fwrite(gxz, sizeof(double), arraylength, outfile);
  fwrite(gyy, sizeof(double), arraylength, outfile);
  fwrite(gyz, sizeof(double), arraylength, outfile);
  fwrite(gzz, sizeof(double), arraylength, outfile);
  fwrite(lapse, sizeof(double), arraylength, outfile);
  fwrite(betax, sizeof(double), arraylength, outfile);
  fwrite(betay, sizeof(double), arraylength, outfile);
  fwrite(betaz, sizeof(double), arraylength, outfile);
  fwrite(kxx, sizeof(double), arraylength, outfile);
  fwrite(kxy, sizeof(double), arraylength, outfile);
  fwrite(kxz, sizeof(double), arraylength, outfile);
  fwrite(kyy, sizeof(double), arraylength, outfile);
  fwrite(kyz, sizeof(double), arraylength, outfile);
  fwrite(kzz, sizeof(double), arraylength, outfile);

  fwrite(rho, sizeof(double), arraylength, outfile);
  fwrite(press, sizeof(double), arraylength, outfile);
  fwrite(vx, sizeof(double), arraylength, outfile);
  fwrite(vy, sizeof(double), arraylength, outfile);
  fwrite(vz, sizeof(double), arraylength, outfile);
  fwrite(Bx, sizeof(double), arraylength, outfile);
  fwrite(By, sizeof(double), arraylength, outfile);
  fwrite(Bz, sizeof(double), arraylength, outfile);
  fwrite(Y_e, sizeof(double), arraylength, outfile);
  fwrite(entropy, sizeof(double), arraylength, outfile);

  fwrite(rho_r, sizeof(double), arraylength, outfile);
  fwrite(press_r, sizeof(double), arraylength, outfile);
  fwrite(vx_r, sizeof(double), arraylength, outfile);
  fwrite(vy_r, sizeof(double), arraylength, outfile);
  fwrite(vz_r, sizeof(double), arraylength, outfile);
  fwrite(Bx_r, sizeof(double), arraylength, outfile);
  fwrite(By_r, sizeof(double), arraylength, outfile);
  fwrite(Bz_r, sizeof(double), arraylength, outfile);
  fwrite(Y_e_r, sizeof(double), arraylength, outfile);
  fwrite(entropy_r, sizeof(double), arraylength, outfile);

  fwrite(rho_l, sizeof(double), arraylength, outfile);
  fwrite(press_l, sizeof(double), arraylength, outfile);
  fwrite(vx_l, sizeof(double), arraylength, outfile);
  fwrite(vy_l, sizeof(double), arraylength, outfile);
  fwrite(vz_l, sizeof(double), arraylength, outfile);
  fwrite(Bx_l, sizeof(double), arraylength, outfile);
  fwrite(By_l, sizeof(double), arraylength, outfile);
  fwrite(Bz_l, sizeof(double), arraylength, outfile);
  fwrite(Y_e_l, sizeof(double), arraylength, outfile);
  fwrite(entropy_l, sizeof(double), arraylength, outfile);

  for(int k = 0; k < dirlength; k++) {
    for(int j = 0; j < dirlength; j++) {
      for(int i = 0; i < dirlength; i++) {
        const int index = indexf(dirlength, i, j, k);

        double dummy = 0;
        rho[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Y_e[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        temperature[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vx[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vy[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vz[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bx[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        By[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bz[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        ghl_tabulated_compute_P_eps_S_from_T(
            &eos, rho[index], Y_e[index], temperature[index], &press[index], &dummy,
            &entropy[index]);

        rho_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Y_e_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        temperature_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vx_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vy_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vz_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bx_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        By_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bz_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        ghl_tabulated_compute_P_eps_S_from_T(
            &eos, rho_r[index], Y_e_r[index], temperature_r[index], &press_r[index], &dummy,
            &entropy_r[index]);

        rho_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Y_e_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        temperature_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vx_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vy_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vz_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bx_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        By_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bz_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        ghl_tabulated_compute_P_eps_S_from_T(
            &eos, rho_l[index], Y_e_l[index], temperature_l[index], &press_l[index], &dummy,
            &entropy_l[index]);
      }
    }
  }

  FILE *outpert = fopen_with_check("ET_Legacy_flux_source_input_pert.bin", "wb");
  fwrite(gxx, sizeof(double), arraylength, outpert);
  fwrite(gxy, sizeof(double), arraylength, outpert);
  fwrite(gxz, sizeof(double), arraylength, outpert);
  fwrite(gyy, sizeof(double), arraylength, outpert);
  fwrite(gyz, sizeof(double), arraylength, outpert);
  fwrite(gzz, sizeof(double), arraylength, outpert);
  fwrite(lapse, sizeof(double), arraylength, outpert);
  fwrite(betax, sizeof(double), arraylength, outpert);
  fwrite(betay, sizeof(double), arraylength, outpert);
  fwrite(betaz, sizeof(double), arraylength, outpert);
  fwrite(kxx, sizeof(double), arraylength, outpert);
  fwrite(kxy, sizeof(double), arraylength, outpert);
  fwrite(kxz, sizeof(double), arraylength, outpert);
  fwrite(kyy, sizeof(double), arraylength, outpert);
  fwrite(kyz, sizeof(double), arraylength, outpert);
  fwrite(kzz, sizeof(double), arraylength, outpert);

  fwrite(rho, sizeof(double), arraylength, outpert);
  fwrite(press, sizeof(double), arraylength, outpert);
  fwrite(vx, sizeof(double), arraylength, outpert);
  fwrite(vy, sizeof(double), arraylength, outpert);
  fwrite(vz, sizeof(double), arraylength, outpert);
  fwrite(Bx, sizeof(double), arraylength, outpert);
  fwrite(By, sizeof(double), arraylength, outpert);
  fwrite(Bz, sizeof(double), arraylength, outpert);
  fwrite(Y_e, sizeof(double), arraylength, outpert);
  fwrite(entropy, sizeof(double), arraylength, outpert);

  fwrite(rho_r, sizeof(double), arraylength, outpert);
  fwrite(press_r, sizeof(double), arraylength, outpert);
  fwrite(vx_r, sizeof(double), arraylength, outpert);
  fwrite(vy_r, sizeof(double), arraylength, outpert);
  fwrite(vz_r, sizeof(double), arraylength, outpert);
  fwrite(Bx_r, sizeof(double), arraylength, outpert);
  fwrite(By_r, sizeof(double), arraylength, outpert);
  fwrite(Bz_r, sizeof(double), arraylength, outpert);
  fwrite(Y_e_r, sizeof(double), arraylength, outpert);
  fwrite(entropy_r, sizeof(double), arraylength, outpert);

  fwrite(rho_l, sizeof(double), arraylength, outpert);
  fwrite(press_l, sizeof(double), arraylength, outpert);
  fwrite(vx_l, sizeof(double), arraylength, outpert);
  fwrite(vy_l, sizeof(double), arraylength, outpert);
  fwrite(vz_l, sizeof(double), arraylength, outpert);
  fwrite(Bx_l, sizeof(double), arraylength, outpert);
  fwrite(By_l, sizeof(double), arraylength, outpert);
  fwrite(Bz_l, sizeof(double), arraylength, outpert);
  fwrite(Y_e_l, sizeof(double), arraylength, outpert);
  fwrite(entropy_l, sizeof(double), arraylength, outpert);

  // Set r/l faces for next direction
  for(int k = 0; k < dirlength; k++) {
    for(int j = 0; j < dirlength; j++) {
      for(int i = 0; i < dirlength; i++) {
        const int index = indexf(dirlength, i, j, k);

        // clang-format off
        randomize_rho_Ye_T_P_S(&eos, &rho_r[index], &Y_e_r[index], &temperature_r[index], &press_r[index], &entropy_r[index]);
        randomize_rho_Ye_T_P_S(&eos, &rho_l[index], &Y_e_l[index], &temperature_l[index], &press_l[index], &entropy_l[index]);
        // clang-format on

        ghl_randomize_primitives(
            &eos, rho_r[index], press_r[index], &vx_r[index], &vy_r[index], &vz_r[index],
            &Bx_r[index], &By_r[index], &Bz_r[index]);

        ghl_randomize_primitives(
            &eos, rho_l[index], press_l[index], &vx_l[index], &vy_l[index], &vz_l[index],
            &Bx_l[index], &By_l[index], &Bz_l[index]);
      }
    }
  }

  fwrite(rho_r, sizeof(double), arraylength, outfile);
  fwrite(press_r, sizeof(double), arraylength, outfile);
  fwrite(vx_r, sizeof(double), arraylength, outfile);
  fwrite(vy_r, sizeof(double), arraylength, outfile);
  fwrite(vz_r, sizeof(double), arraylength, outfile);
  fwrite(Bx_r, sizeof(double), arraylength, outfile);
  fwrite(By_r, sizeof(double), arraylength, outfile);
  fwrite(Bz_r, sizeof(double), arraylength, outfile);
  fwrite(Y_e_r, sizeof(double), arraylength, outfile);
  fwrite(entropy_r, sizeof(double), arraylength, outfile);

  fwrite(rho_l, sizeof(double), arraylength, outfile);
  fwrite(press_l, sizeof(double), arraylength, outfile);
  fwrite(vx_l, sizeof(double), arraylength, outfile);
  fwrite(vy_l, sizeof(double), arraylength, outfile);
  fwrite(vz_l, sizeof(double), arraylength, outfile);
  fwrite(Bx_l, sizeof(double), arraylength, outfile);
  fwrite(By_l, sizeof(double), arraylength, outfile);
  fwrite(Bz_l, sizeof(double), arraylength, outfile);
  fwrite(Y_e_l, sizeof(double), arraylength, outfile);
  fwrite(entropy_l, sizeof(double), arraylength, outfile);

  for(int k = 0; k < dirlength; k++) {
    for(int j = 0; j < dirlength; j++) {
      for(int i = 0; i < dirlength; i++) {
        const int index = indexf(dirlength, i, j, k);

        double dummy = 0;

        rho_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Y_e_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        temperature_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vx_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vy_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vz_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bx_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        By_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bz_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        ghl_tabulated_compute_P_eps_S_from_T(
            &eos, rho_r[index], Y_e_r[index], temperature_r[index], &press_r[index], &dummy,
            &entropy_r[index]);

        rho_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Y_e_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        temperature_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vx_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vy_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vz_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bx_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        By_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bz_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        ghl_tabulated_compute_P_eps_S_from_T(
            &eos, rho_l[index], Y_e_l[index], temperature_l[index], &press_l[index], &dummy,
            &entropy_l[index]);
      }
    }
  }

  fwrite(rho_r, sizeof(double), arraylength, outpert);
  fwrite(press_r, sizeof(double), arraylength, outpert);
  fwrite(vx_r, sizeof(double), arraylength, outpert);
  fwrite(vy_r, sizeof(double), arraylength, outpert);
  fwrite(vz_r, sizeof(double), arraylength, outpert);
  fwrite(Bx_r, sizeof(double), arraylength, outpert);
  fwrite(By_r, sizeof(double), arraylength, outpert);
  fwrite(Bz_r, sizeof(double), arraylength, outpert);
  fwrite(Y_e_r, sizeof(double), arraylength, outpert);
  fwrite(entropy_r, sizeof(double), arraylength, outpert);

  fwrite(rho_l, sizeof(double), arraylength, outpert);
  fwrite(press_l, sizeof(double), arraylength, outpert);
  fwrite(vx_l, sizeof(double), arraylength, outpert);
  fwrite(vy_l, sizeof(double), arraylength, outpert);
  fwrite(vz_l, sizeof(double), arraylength, outpert);
  fwrite(Bx_l, sizeof(double), arraylength, outpert);
  fwrite(By_l, sizeof(double), arraylength, outpert);
  fwrite(Bz_l, sizeof(double), arraylength, outpert);
  fwrite(Y_e_l, sizeof(double), arraylength, outpert);
  fwrite(entropy_l, sizeof(double), arraylength, outpert);

  // Set r/l faces for next direction
  for(int k = 0; k < dirlength; k++) {
    for(int j = 0; j < dirlength; j++) {
      for(int i = 0; i < dirlength; i++) {
        const int index = indexf(dirlength, i, j, k);

        // clang-format off
        randomize_rho_Ye_T_P_S(&eos, &rho_r[index], &Y_e_r[index], &temperature_r[index], &press_r[index], &entropy_r[index]);
        randomize_rho_Ye_T_P_S(&eos, &rho_l[index], &Y_e_l[index], &temperature_l[index], &press_l[index], &entropy_l[index]);
        // clang-format on

        ghl_randomize_primitives(
            &eos, rho_r[index], press_r[index], &vx_r[index], &vy_r[index], &vz_r[index],
            &Bx_r[index], &By_r[index], &Bz_r[index]);

        ghl_randomize_primitives(
            &eos, rho_l[index], press_l[index], &vx_l[index], &vy_l[index], &vz_l[index],
            &Bx_l[index], &By_l[index], &Bz_l[index]);
      }
    }
  }

  fwrite(rho_r, sizeof(double), arraylength, outfile);
  fwrite(press_r, sizeof(double), arraylength, outfile);
  fwrite(vx_r, sizeof(double), arraylength, outfile);
  fwrite(vy_r, sizeof(double), arraylength, outfile);
  fwrite(vz_r, sizeof(double), arraylength, outfile);
  fwrite(Bx_r, sizeof(double), arraylength, outfile);
  fwrite(By_r, sizeof(double), arraylength, outfile);
  fwrite(Bz_r, sizeof(double), arraylength, outfile);
  fwrite(Y_e_r, sizeof(double), arraylength, outfile);
  fwrite(entropy_r, sizeof(double), arraylength, outfile);

  fwrite(rho_l, sizeof(double), arraylength, outfile);
  fwrite(press_l, sizeof(double), arraylength, outfile);
  fwrite(vx_l, sizeof(double), arraylength, outfile);
  fwrite(vy_l, sizeof(double), arraylength, outfile);
  fwrite(vz_l, sizeof(double), arraylength, outfile);
  fwrite(Bx_l, sizeof(double), arraylength, outfile);
  fwrite(By_l, sizeof(double), arraylength, outfile);
  fwrite(Bz_l, sizeof(double), arraylength, outfile);
  fwrite(Y_e_l, sizeof(double), arraylength, outfile);
  fwrite(entropy_l, sizeof(double), arraylength, outfile);

  for(int k = 0; k < dirlength; k++) {
    for(int j = 0; j < dirlength; j++) {
      for(int i = 0; i < dirlength; i++) {
        const int index = indexf(dirlength, i, j, k);

        double dummy = 0;

        rho_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Y_e_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        temperature_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vx_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vy_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vz_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bx_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        By_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bz_r[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        ghl_tabulated_compute_P_eps_S_from_T(
            &eos, rho_r[index], Y_e_r[index], temperature_r[index], &press_r[index], &dummy,
            &entropy_r[index]);

        rho_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Y_e_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        temperature_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vx_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vy_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        vz_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bx_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        By_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        Bz_l[index] *= 1 + randf(-1.0, 1.0) * 1.0e-14;
        ghl_tabulated_compute_P_eps_S_from_T(
            &eos, rho_l[index], Y_e_l[index], temperature_l[index], &press_l[index], &dummy,
            &entropy_l[index]);
      }
    }
  }

  fwrite(rho_r, sizeof(double), arraylength, outpert);
  fwrite(press_r, sizeof(double), arraylength, outpert);
  fwrite(vx_r, sizeof(double), arraylength, outpert);
  fwrite(vy_r, sizeof(double), arraylength, outpert);
  fwrite(vz_r, sizeof(double), arraylength, outpert);
  fwrite(Bx_r, sizeof(double), arraylength, outpert);
  fwrite(By_r, sizeof(double), arraylength, outpert);
  fwrite(Bz_r, sizeof(double), arraylength, outpert);
  fwrite(Y_e_r, sizeof(double), arraylength, outpert);
  fwrite(entropy_r, sizeof(double), arraylength, outpert);

  fwrite(rho_l, sizeof(double), arraylength, outpert);
  fwrite(press_l, sizeof(double), arraylength, outpert);
  fwrite(vx_l, sizeof(double), arraylength, outpert);
  fwrite(vy_l, sizeof(double), arraylength, outpert);
  fwrite(vz_l, sizeof(double), arraylength, outpert);
  fwrite(Bx_l, sizeof(double), arraylength, outpert);
  fwrite(By_l, sizeof(double), arraylength, outpert);
  fwrite(Bz_l, sizeof(double), arraylength, outpert);
  fwrite(Y_e_l, sizeof(double), arraylength, outpert);
  fwrite(entropy_l, sizeof(double), arraylength, outpert);

  fclose(outfile);
  fclose(outpert);
}
