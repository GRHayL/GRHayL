#include "unit_tests.h"

int main(int argc, char **argv) {
  const double poison = 1e200;

  const double W_max = 10.0;
  const int neos = 1;
  const double rho_ppoly_in[1] = {0.0};
  const double Gamma_ppoly_in[1] = {2.0};
  const double k_ppoly0 = 1.0;
  const double Gamma_th = 2.0;

  eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
                                             poison, poison, poison,
                                             neos, rho_ppoly_in, Gamma_ppoly_in,
                                             k_ppoly0, Gamma_th, &eos);

  FILE* infile = fopen("superbee_reconstruction_input.bin", "rb");
  check_file_was_successfully_open(infile, "superbee_reconstruction_input.bin");

  int dirlength;
  int key = fread(&dirlength, sizeof(int), 1, infile);
  const int arraylength = dirlength*dirlength*dirlength;

  double *rho = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);

  double *rhor_trusted   = (double*) malloc(sizeof(double)*arraylength);
  double *rhol_trusted   = (double*) malloc(sizeof(double)*arraylength);
  double *pressr_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *pressl_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vxr_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *vxl_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *vyr_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *vyl_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *vzr_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *vzl_trusted    = (double*) malloc(sizeof(double)*arraylength);

  double *rhor_pert   = (double*) malloc(sizeof(double)*arraylength);
  double *rhol_pert   = (double*) malloc(sizeof(double)*arraylength);
  double *pressr_pert = (double*) malloc(sizeof(double)*arraylength);
  double *pressl_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vxr_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vxl_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vyr_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vyl_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vzr_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vzl_pert    = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(rho  , sizeof(double), arraylength, infile);
  key += fread(press, sizeof(double), arraylength, infile);
  key += fread(vx   , sizeof(double), arraylength, infile);
  key += fread(vy   , sizeof(double), arraylength, infile);
  key += fread(vz   , sizeof(double), arraylength, infile);
  fclose(infile);

  if(key != arraylength*5)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen("superbee_reconstruction_output.bin","rb");
  check_file_was_successfully_open(infile, "superbee_reconstruction_output.bin");
  key  = fread(rhor_trusted  , sizeof(double), arraylength, infile);
  key += fread(rhol_trusted  , sizeof(double), arraylength, infile);
  key += fread(pressr_trusted, sizeof(double), arraylength, infile);
  key += fread(pressl_trusted, sizeof(double), arraylength, infile);
  key += fread(vxr_trusted   , sizeof(double), arraylength, infile);
  key += fread(vxl_trusted   , sizeof(double), arraylength, infile);
  key += fread(vyr_trusted   , sizeof(double), arraylength, infile);
  key += fread(vyl_trusted   , sizeof(double), arraylength, infile);
  key += fread(vzr_trusted   , sizeof(double), arraylength, infile);
  key += fread(vzl_trusted   , sizeof(double), arraylength, infile);
  fclose(infile);

  if(key != arraylength*10)
    ghl_error("An error has occured with reading in comparison data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  FILE* inpert = fopen("superbee_reconstruction_output_pert.bin","rb");
  check_file_was_successfully_open(inpert, "superbee_reconstruction_output_pert.bin");

  key  = fread(rhor_pert  , sizeof(double), arraylength, inpert);
  key += fread(rhol_pert  , sizeof(double), arraylength, inpert);
  key += fread(pressr_pert, sizeof(double), arraylength, inpert);
  key += fread(pressl_pert, sizeof(double), arraylength, inpert);
  key += fread(vxr_pert   , sizeof(double), arraylength, inpert);
  key += fread(vxl_pert   , sizeof(double), arraylength, inpert);
  key += fread(vyr_pert   , sizeof(double), arraylength, inpert);
  key += fread(vyl_pert   , sizeof(double), arraylength, inpert);
  key += fread(vzr_pert   , sizeof(double), arraylength, inpert);
  key += fread(vzl_pert   , sizeof(double), arraylength, inpert);
  fclose(inpert);

  if(key != arraylength*10)
    ghl_error("An error has occured with reading in comparison data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // These are set up to match the loops in the ET version of IllinoisGRMHD.
  const int num_vars = 5;
  double *gfs[5] = {rho, press, vx, vy, vz};
  
  double gfs_r[5];
  double gfs_l[5]; 

  for(int k=2; k<dirlength-1; k++)
    for(int j=2; j<dirlength-1; j++)
      for(int i=2; i<dirlength-1; i++) {
        const int index = indexf(dirlength, i, j, k);
        double var_data[num_vars][4];

        for(int ind=0; ind<4; ind++) {
          const int stencil = indexf(dirlength, i, j+(ind-2), k);
          for(int var=0; var<num_vars; var++) {
            var_data[var][ind] = gfs[var][stencil];
          }
        }        
        
    for(int var=0; var<num_vars; var++) {
        ghl_superbee_reconstruction( var_data[var][0], //U_m2
                                        var_data[var][1], //U_m1
                                        var_data[var][2], //U
                                        var_data[var][3], //U_p1
                                        &gfs_r[var],
                                        &gfs_l[var]);
    }

      if( validate(rhor_trusted[index], gfs_r[0], rhor_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable rho_r.\n"
                       "  rho_r trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", rhor_trusted[index], gfs_r[0], rhor_pert[index],
                                                   relative_error(rhor_trusted[index], gfs_r[0]),
                                                   relative_error(rhor_trusted[index], rhor_pert[index]));
        if( validate(pressr_trusted[index], gfs_r[1], pressr_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable press_r.\n"
                       "  press_r trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", pressr_trusted[index], gfs_r[1], pressr_pert[index],
                                                   relative_error(pressr_trusted[index], gfs_r[1]),
                                                   relative_error(pressr_trusted[index], pressr_pert[index]));
        if( validate(vxr_trusted[index], gfs_r[2], vxr_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable vx_r.\n"
                       "  vx_r trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", vxr_trusted[index], gfs_r[2], vxr_pert[index],
                                                   relative_error(vxr_trusted[index], gfs_r[2]),
                                                   relative_error(vxr_trusted[index], vxr_pert[index]));
        if( validate(vyr_trusted[index], gfs_r[3], vyr_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable vy_r.\n"
                       "  vy_r trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", vyr_trusted[index], gfs_r[3], vyr_pert[index],
                                                   relative_error(vyr_trusted[index], gfs_r[3]),
                                                   relative_error(vyr_trusted[index], vyr_pert[index]));
        if( validate(vzr_trusted[index], gfs_r[4], vzr_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable vz_r.\n"
                       "  vz_r trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", vzr_trusted[index], gfs_r[4], vzr_pert[index],
                                                   relative_error(vzr_trusted[index], gfs_r[4]),
                                                   relative_error(vzr_trusted[index], vzr_pert[index]));

        if( validate(rhol_trusted[index], gfs_l[0], rhol_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable rho_l.\n"
                       "  rho_l trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", rhol_trusted[index], gfs_l[0], rhol_pert[index],
                                                   relative_error(rhol_trusted[index], gfs_l[0]),
                                                   relative_error(rhol_trusted[index], rhol_pert[index]));
        if( validate(pressr_trusted[index], gfs_l[1], pressl_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable press_l.\n"
                       "  press_l trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", pressl_trusted[index], gfs_l[1], pressl_pert[index],
                                                   relative_error(pressl_trusted[index], gfs_l[1]),
                                                   relative_error(pressl_trusted[index], pressl_pert[index]));
        if( validate(vxl_trusted[index], gfs_l[2], vxl_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable vx_l.\n"
                       "  vx_l trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", vxl_trusted[index], gfs_l[2], vxl_pert[index],
                                                   relative_error(vxl_trusted[index], gfs_l[2]),
                                                   relative_error(vxl_trusted[index], vxl_pert[index]));
        if( validate(vyl_trusted[index], gfs_l[3], vyl_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable vy_l.\n"
                       "  vy_l trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", vyl_trusted[index], gfs_l[3], vyl_pert[index],
                                                   relative_error(vyl_trusted[index], gfs_l[3]),
                                                   relative_error(vyl_trusted[index], vyl_pert[index]));
        if( validate(vzl_trusted[index], gfs_l[4], vzl_pert[index]) )
          ghl_error("Test unit_test_superbee_reconstruction has failed for variable vz_l.\n"
                       "  vz_l trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", vzl_trusted[index], gfs_l[4], vzl_pert[index],
                                                   relative_error(vzl_trusted[index], gfs_l[4]),
                                                   relative_error(vzl_trusted[index], vzl_pert[index]));
  }

  ghl_info("superbee reconstruction test has passed!\n");
  free(rho); free(press);
  free(vx); free(vy); free(vz);

  free(rhor_trusted); free(rhol_trusted);
  free(pressr_trusted); free(pressl_trusted);
  free(vxr_trusted); free(vxl_trusted); free(vyr_trusted);
  free(vyl_trusted); free(vzr_trusted); free(vzl_trusted);

  free(rhor_pert); free(rhol_pert);
  free(pressr_pert); free(pressl_pert);
  free(vxr_pert); free(vxl_pert); free(vyr_pert);
  free(vyl_pert); free(vzr_pert); free(vzl_pert);
}
