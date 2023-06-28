#include "unit_tests.h"

int main(int argc, char **argv) {
  FILE* infile = fopen_with_check("induction_interpolation_input.bin","rb");

  int dirlength;
  int key = fread(&dirlength, sizeof(int), 1, infile);
  if( key != 1 || dirlength < 1 )
    ghl_error("An error has occured with reading the grid size. "
                 "Please check that Noble2D_initial_data.bin"
                 "is up-to-date with current test version.\n");
  const int arraylength = dirlength*dirlength*dirlength;

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

  // Read in data from file to ensure portability
  key  = fread(lapse, sizeof(double), arraylength, infile);
  key += fread(betax, sizeof(double), arraylength, infile);
  key += fread(betay, sizeof(double), arraylength, infile);
  key += fread(betaz, sizeof(double), arraylength, infile);

  key += fread(phitilde, sizeof(double), arraylength, infile);
  key += fread(Ax,       sizeof(double), arraylength, infile);
  key += fread(Ay,       sizeof(double), arraylength, infile);
  key += fread(Az,       sizeof(double), arraylength, infile);
  fclose(infile);
  if(key != arraylength*8)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen_with_check("induction_interpolation_ADM_input.bin","rb");
  key  = fread(gxx, sizeof(double), arraylength, infile);
  key += fread(gxy, sizeof(double), arraylength, infile);
  key += fread(gxz, sizeof(double), arraylength, infile);
  key += fread(gyy, sizeof(double), arraylength, infile);
  key += fread(gyz, sizeof(double), arraylength, infile);
  key += fread(gzz, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*6)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Data which should be written before it is used is poisoned to validate behavior.
  // RHSs for A are set to 0 because they are assumed to already contain the zero-gauge
  // contribution to the RHS before entering this function.
  const double poison = 1e300;

  for(int k=0; k<dirlength; k++) {
    for(int j=0; j<dirlength; j++) {
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        alpha_interp[index] = poison;
        betax_interp[index] = poison;
        betay_interp[index] = poison;
        betaz_interp[index] = poison;

        alpha_Phi_minus_betaj_A_j_interp[index] = poison;
        sqrtg_Ax_interp[index]                  = poison;
        sqrtg_Ay_interp[index]                  = poison;
        sqrtg_Az_interp[index]                  = poison;
      }
    }
  }

  ghl_test_compute_ccc_ADM(
        dirlength, lapse, betax, betay, betaz,
        gxx, gxy, gxz, gyy, gyz, gzz,
        phitilde, Ax, Ay, Az,
        alpha_interp, betax_interp,
        betay_interp, betaz_interp,
        alpha_Phi_minus_betaj_A_j_interp,
        sqrtg_Ax_interp, sqrtg_Ay_interp,
        sqrtg_Az_interp);

  infile = fopen_with_check("induction_interpolation_ccc_ADM_output.bin", "rb");

  double *alpha_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *betax_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *betay_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *betaz_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *phitilde_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_trusted       = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_trusted       = (double*) malloc(sizeof(double)*arraylength);
  double *Az_trusted       = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(alpha_trusted,    sizeof(double), arraylength, infile);
  key += fread(betax_trusted,    sizeof(double), arraylength, infile);
  key += fread(betay_trusted,    sizeof(double), arraylength, infile);
  key += fread(betaz_trusted,    sizeof(double), arraylength, infile);
  key += fread(phitilde_trusted, sizeof(double), arraylength, infile);
  key += fread(Ax_trusted,       sizeof(double), arraylength, infile);
  key += fread(Ay_trusted,       sizeof(double), arraylength, infile);
  key += fread(Az_trusted,       sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*8)
    ghl_error("An error has occured with reading in trusted data. Please check that comparison data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen_with_check("induction_interpolation_ccc_ADM_output_pert.bin", "rb");

  double *alpha_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *betax_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *betay_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *betaz_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *phitilde_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_pert       = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_pert       = (double*) malloc(sizeof(double)*arraylength);
  double *Az_pert       = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(alpha_pert,    sizeof(double), arraylength, infile);
  key += fread(betax_pert,    sizeof(double), arraylength, infile);
  key += fread(betay_pert,    sizeof(double), arraylength, infile);
  key += fread(betaz_pert,    sizeof(double), arraylength, infile);
  key += fread(phitilde_pert, sizeof(double), arraylength, infile);
  key += fread(Ax_pert,       sizeof(double), arraylength, infile);
  key += fread(Ay_pert,       sizeof(double), arraylength, infile);
  key += fread(Az_pert,       sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*8)
    ghl_error("An error has occured with reading in perturbed data. Please check that comparison data\n"
                 "is up-to-date with current test version.\n");

  for(int k=1; k<dirlength-1; k++) {
    for(int j=1; j<dirlength-1; j++) {
      for(int i=1; i<dirlength-1; i++) {
        const int index = indexf(dirlength,i,j,k);

        if( validate(alpha_trusted[index], alpha_interp[index], alpha_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of alpha.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", alpha_trusted[index], alpha_interp[index], alpha_pert[index],
                                                   relative_error(alpha_trusted[index], alpha_interp[index]),
                                                   relative_error(alpha_trusted[index], alpha_pert[index]));

        if( validate(betax_trusted[index], betax_interp[index], betax_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of betax.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", betax_trusted[index], betax_interp[index], betax_pert[index],
                                                   relative_error(betax_trusted[index], betax_interp[index]),
                                                   relative_error(betax_trusted[index], betax_pert[index]));

        if( validate(betay_trusted[index], betay_interp[index], betay_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of betay.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", betay_trusted[index], betay_interp[index], betay_pert[index],
                                                   relative_error(betay_trusted[index], betay_interp[index]),
                                                   relative_error(betay_trusted[index], betay_pert[index]));

        if( validate(betaz_trusted[index], betaz_interp[index], betaz_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of betaz.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", betaz_trusted[index], betaz_interp[index], betaz_pert[index],
                                                   relative_error(betaz_trusted[index], betaz_interp[index]),
                                                   relative_error(betaz_trusted[index], betaz_pert[index]));

        if( validate(phitilde_trusted[index], alpha_Phi_minus_betaj_A_j_interp[index], phitilde_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of alpha_Phi_minus_betaj_A_j.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", phitilde_trusted[index], alpha_Phi_minus_betaj_A_j_interp[index], phitilde_pert[index],
                                                   relative_error(phitilde_trusted[index], alpha_Phi_minus_betaj_A_j_interp[index]),
                                                   relative_error(phitilde_trusted[index], phitilde_pert[index]));

        if( validate(Ax_trusted[index], sqrtg_Ax_interp[index], Ax_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of Ax.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", Ax_trusted[index], sqrtg_Ax_interp[index], Ax_pert[index],
                                                   relative_error(Ax_trusted[index], sqrtg_Ax_interp[index]),
                                                   relative_error(Ax_trusted[index], Ax_pert[index]));

        if( validate(Ay_trusted[index], sqrtg_Ay_interp[index], Ay_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of Ay.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", Ay_trusted[index], sqrtg_Ay_interp[index], Ay_pert[index],
                                                   relative_error(Ay_trusted[index], sqrtg_Ay_interp[index]),
                                                   relative_error(Ay_trusted[index], Ay_pert[index]));

        if( validate(Az_trusted[index], sqrtg_Az_interp[index], Az_pert[index]) )
          ghl_error("Test unit_test_induction_interpolation_ccc_ADM has failed for interpolation of Az.\n"
                       "  trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", Az_trusted[index], sqrtg_Az_interp[index], Az_pert[index],
                                                   relative_error(Az_trusted[index], sqrtg_Az_interp[index]),
                                                   relative_error(Az_trusted[index], Az_pert[index]));
      }
    }
  }
  ghl_info("Interpolation test for cell-centered ADM inputs has passed!\n");
  free(lapse); free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(Ax); free(Ay); free(Az); free(phitilde);
  free(alpha_interp);
  free(betax_interp); free(betay_interp); free(betaz_interp);
  free(sqrtg_Ax_interp); free(sqrtg_Ay_interp); free(sqrtg_Az_interp);
  free(alpha_Phi_minus_betaj_A_j_interp);
  free(alpha_trusted);
  free(betax_trusted); free(betay_trusted); free(betaz_trusted);
  free(Ax_trusted); free(Ay_trusted); free(Az_trusted);
  free(phitilde_trusted);
  free(alpha_pert);
  free(betax_pert); free(betay_pert); free(betaz_pert);
  free(Ax_pert); free(Ay_pert); free(Az_pert);
  free(phitilde_pert);
}
