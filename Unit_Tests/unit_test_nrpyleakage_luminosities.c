#include "neutrinos.h"
#include "unit_tests.h"

void
generate_test_data(const ghl_eos_parameters *restrict eos) {

  srand(100);

  const int npoints = 1024;

  for(int perturb=0;perturb<=1;perturb++) {
    char filename[64];
    if( perturb )
      sprintf(filename, "nrpyleakage_luminosities_perturbed.bin");
    else
      sprintf(filename, "nrpyleakage_luminosities_unperturbed.bin");

    FILE *fp = fopen_with_check(filename, "wb");
    fwrite(&npoints, sizeof(int), 1, fp);
    for(int n=0;n<npoints;n++) {

      // Get random metric values
      double alpha;
      __attribute__((unused)) double betax, betay, betaz;
      double gammaxx, gammaxy, gammaxz, gammayy, gammayz, gammazz;
      ghl_randomize_metric(
            &alpha, &betax, &betay, &betaz,
            &gammaxx, &gammaxy, &gammaxz,
            &gammayy, &gammayz, &gammazz);

      // Get random primitive values
      double rho = pow(10, randf(log10(eos->rho_min), log10(eos->rho_max)));
      double Y_e = randf(eos->Y_e_min, eos->Y_e_max);
      double T   = pow(10, randf(log10(eos->T_min), log10(eos->T_max)));
      double W   = randf(1, 10);

      // Get random optical depths (not sure these are reasonable values)
      ghl_neutrino_optical_depths tau;
      tau.nue [0] = randf(1, 1000);
      tau.nue [1] = randf(1, 1000);
      tau.anue[0] = randf(1, 1000);
      tau.anue[1] = randf(1, 1000);
      tau.nux [0] = randf(1, 1000);
      tau.nux [1] = randf(1, 1000);

      if( perturb ) {
        alpha       *= (1+randf(-1,1)*1e-14);
        gammaxx     *= (1+randf(-1,1)*1e-14);
        gammaxy     *= (1+randf(-1,1)*1e-14);
        gammaxz     *= (1+randf(-1,1)*1e-14);
        gammayy     *= (1+randf(-1,1)*1e-14);
        gammayz     *= (1+randf(-1,1)*1e-14);
        gammazz     *= (1+randf(-1,1)*1e-14);
        rho         *= (1+randf(-1,1)*1e-14);
        Y_e         *= (1+randf(-1,1)*1e-14);
        T           *= (1+randf(-1,1)*1e-14);
        W           *= (1+randf(-1,1)*1e-14);
        tau.nue [0] *= (1+randf(-1,1)*1e-14);
        tau.nue [1] *= (1+randf(-1,1)*1e-14);
        tau.anue[0] *= (1+randf(-1,1)*1e-14);
        tau.anue[1] *= (1+randf(-1,1)*1e-14);
        tau.nux [0] *= (1+randf(-1,1)*1e-14);
        tau.nux [1] *= (1+randf(-1,1)*1e-14);
      }

      // Compute luminosities
      ghl_neutrino_luminosities lum;
      NRPyLeakage_compute_ghl_neutrino_luminosities(eos, alpha,
                                                gammaxx, gammaxy, gammaxz,
                                                gammayy, gammayz, gammazz,
                                                rho, Y_e, T, W,
                                                &tau, &lum);

      // Output to file
      if( !perturb ) {
        fwrite(&alpha  , sizeof(double)                 , 1, fp);
        fwrite(&gammaxx, sizeof(double)                 , 1, fp);
        fwrite(&gammaxy, sizeof(double)                 , 1, fp);
        fwrite(&gammaxz, sizeof(double)                 , 1, fp);
        fwrite(&gammayy, sizeof(double)                 , 1, fp);
        fwrite(&gammayz, sizeof(double)                 , 1, fp);
        fwrite(&gammazz, sizeof(double)                 , 1, fp);
        fwrite(&rho    , sizeof(double)                 , 1, fp);
        fwrite(&Y_e    , sizeof(double)                 , 1, fp);
        fwrite(&T      , sizeof(double)                 , 1, fp);
        fwrite(&W      , sizeof(double)                 , 1, fp);
        fwrite(&tau    , sizeof(ghl_neutrino_optical_depths), 1, fp);
      }
      fwrite(&lum      , sizeof(ghl_neutrino_luminosities)  , 1, fp);
    }
    fclose(fp);
  }
}

void
run_unit_test(const ghl_eos_parameters *restrict eos) {

  int n1, n2;

  FILE *fp_unpert = fopen_with_check("nrpyleakage_luminosities_unperturbed.bin", "rb");
  FILE *fp_pert   = fopen_with_check("nrpyleakage_luminosities_perturbed.bin", "rb");

  int err = 0;
  err += fread(&n1, sizeof(int), 1, fp_unpert);
  err += fread(&n2, sizeof(int), 1, fp_pert  );
  if( err != 2 || n1 != n2 ) {
    fclose(fp_unpert);
    fclose(fp_pert);
    ghl_error("Problem reading number of points from file (err: %d, n1: %d, n2: %d)\n",
                 err, n1, n2);
  }

  const int npoints=n1;
  for(int n=0;n<npoints;n++) {

    // Read metric and primitive quantities from the unperturbed data file
    double alpha;
    double gammaxx, gammaxy, gammaxz, gammayy, gammayz, gammazz;
    double rho, Y_e, T, W;
    ghl_neutrino_optical_depths tau;

    err  = 0;
    err += fread(&alpha  , sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammaxx, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammaxy, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammaxz, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammayy, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammayz, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammazz, sizeof(double)                 , 1, fp_unpert);
    err += fread(&rho    , sizeof(double)                 , 1, fp_unpert);
    err += fread(&Y_e    , sizeof(double)                 , 1, fp_unpert);
    err += fread(&T      , sizeof(double)                 , 1, fp_unpert);
    err += fread(&W      , sizeof(double)                 , 1, fp_unpert);
    err += fread(&tau    , sizeof(ghl_neutrino_optical_depths), 1, fp_unpert);

    if( err != 12 ) {
      fclose(fp_unpert); fclose(fp_pert);
      ghl_error("Failed to read inputs from unperturbed data file\n");
    }

    // Compute luminosities
    ghl_neutrino_luminosities lum;
    NRPyLeakage_compute_ghl_neutrino_luminosities(eos, alpha,
                                              gammaxx, gammaxy, gammaxz,
                                              gammayy, gammayz, gammazz,
                                              rho, Y_e, T, W,
                                              &tau, &lum);

    // Now read luminosities from unperturbed and perturbed data files
    ghl_neutrino_luminosities lum_trusted, lum_pert;
    if( 1 != fread(&lum_trusted, sizeof(ghl_neutrino_luminosities), 1, fp_unpert) ) {
      fclose(fp_unpert); fclose(fp_pert);
      ghl_error("Failed to read luminosities from unperturbed data file\n");
    }

    if( 1 != fread(&lum_pert, sizeof(ghl_neutrino_luminosities), 1, fp_pert) ) {
      fclose(fp_unpert); fclose(fp_pert);
      ghl_error("Failed to read luminosities from perturbed data file\n");
    }

    validate(lum_trusted.nue , lum.nue , lum_pert.nue );
    validate(lum_trusted.anue, lum.anue, lum_pert.anue);
    validate(lum_trusted.nux , lum.nux , lum_pert.nux );
  }
  fclose(fp_unpert);
  fclose(fp_pert);
}

#include "nrpyleakage_main.h"
