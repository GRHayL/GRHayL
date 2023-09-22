#include "unit_tests.h"

int main(int argc, char **argv) {
  FILE* infile = fopen_with_check("PLM_reconstruction_input.bin", "rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);

  double *var = (double*) malloc(sizeof(double)*arraylength);

  double *varr_trusted   = (double*) malloc(sizeof(double)*arraylength);
  double *varl_trusted   = (double*) malloc(sizeof(double)*arraylength);
  double *varr_pert   = (double*) malloc(sizeof(double)*arraylength);
  double *varl_pert   = (double*) malloc(sizeof(double)*arraylength);

  key = fread(var, sizeof(double), arraylength, infile);
  fclose(infile);

  if(key != arraylength)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen_with_check("PLM_reconstruction_output.bin","rb");
  FILE* inpert = fopen_with_check("PLM_reconstruction_output_pert.bin","rb");

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

    key  = fread(varr_trusted  , sizeof(double), arraylength, infile);
    key += fread(varl_trusted  , sizeof(double), arraylength, infile);

    if(key != arraylength*2)
      ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                   "is up-to-date with current test version.\n");


    key  = fread(varr_pert  , sizeof(double), arraylength, inpert);
    key += fread(varl_pert  , sizeof(double), arraylength, inpert);

    if(key != arraylength*2)
      ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                   "is up-to-date with current test version.\n");

    // These are set up to match the loops in the ET version of IllinoisGRMHD.
    for(int index=2; index<arraylength-1; index++) {
      double var_r, var_l;
      ghl_reconstruction(var[index-2], var[index-1], var[index], var[index+1], &var_r, &var_l);

      char mname[20];
      if(method==0) {
        sprintf(mname, "minmod");
      } else if(method==1) {
        sprintf(mname, "mc");
      } else {
        sprintf(mname, "superbee");
      }
      if( ghl_pert_test_fail(varr_trusted[index], var_r, varr_pert[index]) )
          ghl_error("Test unit_test_reconstruction has failed for method %.20s.\n"
                       "  Right face trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", mname, varr_trusted[index], var_r, varr_pert[index],
                                                   relative_error(varr_trusted[index], var_r),
                                                   relative_error(varr_trusted[index], varr_pert[index]));

      if( ghl_pert_test_fail(varl_trusted[index], var_l, varl_pert[index]) )
        ghl_error("Test unit_test_superbee_reconstruction has failed for method %.20s.\n"
                     "  Left face trusted %.14e computed %.14e perturbed %.14e\n"
                     "  rel.err. %.14e %.14e\n", mname, varl_trusted[index], var_l, varl_pert[index],
                                                 relative_error(varl_trusted[index], var_l),
                                                 relative_error(varl_trusted[index], varl_pert[index]));
    }
  }
  fclose(infile);
  fclose(inpert);

  ghl_info("PLM reconstruction test has passed!\n");
  free(var);
  free(varr_trusted); free(varl_trusted);
  free(varr_pert); free(varl_pert);
}
