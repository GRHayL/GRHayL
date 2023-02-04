#include "unit_tests.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

static double eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in);

int main(int argc, char **argv) {
  const double poison = 1e200;


  const int eos_type = 0;
  const double W_max = 10.0;

  eos_parameters eos;
  initialize_general_eos(eos_type, W_max,
             poison, poison, poison,
             &eos);

  const int neos = 1;
  const double rho_ppoly_in[1] = {0.0};
  const double Gamma_ppoly_in[1] = {2.0};
  const double k_ppoly0 = 1.0;
  const double Gamma_th = 2.0;

  initialize_hybrid_functions(&eos);
  initialize_hybrid_eos(
             neos, rho_ppoly_in,
             Gamma_ppoly_in, k_ppoly0,
             Gamma_th, &eos);

  const int dirlength = 20;
  const int arraylength = dirlength*dirlength*dirlength;

  double *rho = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);

  double *rhor_trusted   = (double*) malloc(sizeof(double)*arraylength);
  double *pressr_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vxr_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *vzr_trusted    = (double*) malloc(sizeof(double)*arraylength);

  double *rhol_trusted   = (double*) malloc(sizeof(double)*arraylength);
  double *pressl_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vxl_trusted    = (double*) malloc(sizeof(double)*arraylength);
  double *vzl_trusted    = (double*) malloc(sizeof(double)*arraylength);

  double *rhor_pert   = (double*) malloc(sizeof(double)*arraylength);
  double *pressr_pert = (double*) malloc(sizeof(double)*arraylength);
  double *rhol_pert   = (double*) malloc(sizeof(double)*arraylength);
  double *pressl_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vxr_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vzr_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vxl_pert    = (double*) malloc(sizeof(double)*arraylength);
  double *vzl_pert = (double*) malloc(sizeof(double)*arraylength);

  FILE* infile;
  infile = fopen("simple_ppm_initial_data.bin", "rb");
  check_file_was_successfully_open(infile, "simple_ppm_initial_data.bin");
  int key = fread(rho  , sizeof(double), arraylength, infile);
  key    += fread(press, sizeof(double), arraylength, infile);
  key    += fread(vx   , sizeof(double), arraylength, infile);
  key    += fread(vz   , sizeof(double), arraylength, infile);
  fclose(infile);
  if(key != arraylength*4)
    grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen("simple_ppm.bin","rb");
  check_file_was_successfully_open(infile, "simple_ppm.bin");
  FILE* inpert = fopen("simple_ppm_pert.bin","rb");
  check_file_was_successfully_open(infile, "simple_ppm_pert.bin");

  key  = fread(rhor_trusted  , sizeof(double), arraylength, infile);
  key += fread(pressr_trusted, sizeof(double), arraylength, infile);
  key += fread(vxr_trusted   , sizeof(double), arraylength, infile);
  key += fread(vzr_trusted   , sizeof(double), arraylength, infile);

  key += fread(rhol_trusted  , sizeof(double), arraylength, infile);
  key += fread(pressl_trusted, sizeof(double), arraylength, infile);
  key += fread(vxl_trusted   , sizeof(double), arraylength, infile);
  key += fread(vzl_trusted   , sizeof(double), arraylength, infile);

  key += fread(rhor_pert  , sizeof(double), arraylength, inpert);
  key += fread(pressr_pert, sizeof(double), arraylength, inpert);
  key += fread(vxr_pert   , sizeof(double), arraylength, inpert);
  key += fread(vzr_pert   , sizeof(double), arraylength, inpert);

  key += fread(rhol_pert  , sizeof(double), arraylength, inpert);
  key += fread(pressl_pert, sizeof(double), arraylength, inpert);
  key += fread(vxl_pert   , sizeof(double), arraylength, inpert);
  key += fread(vzl_pert   , sizeof(double), arraylength, inpert);
  if(key != arraylength*16)
    grhayl_error("An error has occured with reading in comparison data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  const int num_vars = 2;

#pragma omp parallel for
  // we are setting the reconstruction direction to x
  for(int k=1; k<dirlength; k++)
    for(int j=1; j<dirlength; j++)
      for(int i=3; i<dirlength-2; i++) {
        const int index = indexf(dirlength,i,j,k);

        double rhor, rhol, pressr, pressl;
        double rho_stencil[6], press_stencil[6], v_flux_dir[6];
        double var_data[num_vars][6], var_datar[num_vars], var_datal[num_vars];

        for(int ind=0; ind<6; ind++) {
          const int stencil  = indexf(dirlength, i+ind-3, j, k); // PPM needs indices from -3 to +2
          v_flux_dir[ind]    = vx[stencil]; // Could be smaller; doesn't use full stencil
          rho_stencil[ind]   = rho[stencil];
          press_stencil[ind] = press[stencil];
          var_data[0][ind]   = vx[stencil];
          var_data[1][ind]   = vz[stencil];
        }
        const double Gamma_eff = eos_Gamma_eff(&eos, rho[index], press[index]);

        simple_ppm(rho_stencil, press_stencil, var_data,
                   num_vars, v_flux_dir, Gamma_eff,
                   &rhor, &rhol, &pressr, &pressl,
                   var_datar, var_datal);

        if( validate(rhor_trusted[index], rhor, rhor_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable rho_r at index (%d,%d,%d).\n", i, j, k);
        if( validate(pressr_trusted[index], pressr, pressr_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable press_r at index (%d,%d,%d).\n", i, j, k);
        if( validate(vxr_trusted[index], var_datar[0], vxr_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vx_r at index (%d,%d,%d).\n", i, j, k);
        if( validate(vzr_trusted[index], var_datar[1], vzr_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vz_r at index (%d,%d,%d).\n", i, j, k);

        if( validate(rhol_trusted[index], rhol, rhol_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable rho_l at index (%d,%d,%d).\n", i, j, k);
        if( validate(pressl_trusted[index], pressl, pressl_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable press_l at index (%d,%d,%d).\n", i, j, k);
        if( validate(vxl_trusted[index], var_datal[0], vxl_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vx_l at index (%d,%d,%d).\n", i, j, k);
        if( validate(vzl_trusted[index], var_datal[1], vzl_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vz_l at index (%d,%d,%d).\n", i, j, k);
  }

  key  = fread(vxr_trusted   , sizeof(double), arraylength, infile);
  key += fread(vzr_trusted   , sizeof(double), arraylength, infile);
  key += fread(vxl_trusted   , sizeof(double), arraylength, infile);
  key += fread(vzl_trusted   , sizeof(double), arraylength, infile);

  key += fread(vxr_pert   , sizeof(double), arraylength, inpert);
  key += fread(vzr_pert   , sizeof(double), arraylength, inpert);
  key += fread(vxl_pert   , sizeof(double), arraylength, inpert);
  key += fread(vzl_pert   , sizeof(double), arraylength, inpert);
  fclose(infile);
  fclose(inpert);
  if(key != arraylength*8)
    grhayl_error("An error has occured with reading in comparison data. Please check that data\n"
                 "is up-to-date with current test version.\n");

#pragma omp parallel for
  // we are setting the reconstruction direction to z
  for(int k=3; k<dirlength-2; k++)
    for(int j=1; j<dirlength; j++)
      for(int i=1; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        double press_stencil[6], v_flux_dir[6];
        double var_data[num_vars][6], var_datar[num_vars], var_datal[num_vars];

        for(int ind=0; ind<6; ind++) {
          const int stencil  = indexf(dirlength, i+ind-3, j, k); // PPM needs indices from -3 to +2
          v_flux_dir[ind]    = vz[stencil]; // Could be smaller; doesn't use full stencil
          press_stencil[ind] = press[stencil];
          var_data[0][ind]   = vx[stencil];
          var_data[1][ind]   = vz[stencil];
        }
        const double Gamma_eff = eos_Gamma_eff(&eos, rho[index], press[index]);
        simple_ppm_no_rho_P(press_stencil, var_data,
                   num_vars, v_flux_dir, Gamma_eff,
                   var_datar, var_datal);

        if( validate(vxr_trusted[index], var_datar[0], vxr_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vx_r at index (%d,%d,%d).\n", i, j, k);
        if( validate(vzr_trusted[index], var_datar[1], vzr_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vz_r at index (%d,%d,%d).\n", i, j, k);

        if( validate(vxl_trusted[index], var_datal[0], vxl_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vx_l at index (%d,%d,%d).\n", i, j, k);
        if( validate(vzl_trusted[index], var_datal[1], vzl_pert[index]) )
          grhayl_error("Test unit_test_simple_ppm has failed for variable vz_l at index (%d,%d,%d).\n", i, j, k);
  }
}

double eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in) {
  double K, Gamma;
  eos->hybrid_get_K_and_Gamma(eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return eos->Gamma_th + (Gamma - eos->Gamma_th)*P_cold/press_in;
}
