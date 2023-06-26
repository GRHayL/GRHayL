#include "unit_tests.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

static double eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in);

int main(int argc, char **argv) {
  const double poison = 1e300;

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

  FILE* infile = fopen_with_check("ET_Legacy_reconstruction_input.bin", "rb");

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

  infile = fopen_with_check("ET_Legacy_reconstruction_output.bin","rb");

  FILE* inpert = fopen_with_check("ET_Legacy_reconstruction_output_pert.bin","rb");

  for(int flux_dirn = 0; flux_dirn<3; flux_dirn++) {
    const int num_vars = 3;

    double *vflux[3] = {vx, vy, vz};
    const int xdir = (flux_dirn==0);
    const int ydir = (flux_dirn==1);
    const int zdir = (flux_dirn==2);

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

    if(key != arraylength*10)
      ghl_error("An error has occured with reading in comparison data. Please check that data\n"
                   "is up-to-date with current test version.\n");

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

    if(key != arraylength*10)
      ghl_error("An error has occured with reading in comparison data. Please check that data\n"
                   "is up-to-date with current test version.\n");

    // These are set up to match the loops in the ET version of IllinoisGRMHD.
    const int imin = 3;
    const int jmin = 3;
    const int kmin = 3;
    const int imax = dirlength - 3 - !xdir;
    const int jmax = dirlength - 3 - !ydir;
    const int kmax = dirlength - 3 - !zdir;

    for(int k=kmin; k<kmax; k++)
      for(int j=jmin; j<jmax; j++)
        for(int i=imin; i<imax; i++) {
          const int index = indexf(dirlength,i,j,k);

          double rhor, rhol, pressr, pressl;
          double rho_stencil[6], press_stencil[6], v_flux_dir[6];
          double var_data[num_vars][6], var_datar[num_vars], var_datal[num_vars];

          for(int ind=0; ind<6; ind++) {
            const int stencil  = indexf(dirlength, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3)); // PPM needs indices from -3 to +2
            v_flux_dir[ind]    = vflux[flux_dirn][stencil]; // Could be smaller; doesn't use full stencil
            rho_stencil[ind]   = rho[stencil];
            press_stencil[ind] = press[stencil];
            var_data[0][ind]   = vx[stencil];
            var_data[1][ind]   = vy[stencil];
            var_data[2][ind]   = vz[stencil];
          }
          const double Gamma_eff = eos_Gamma_eff(&eos, rho[index], press[index]);

          ghl_ppm(rho_stencil, press_stencil, var_data,
                     num_vars, v_flux_dir, Gamma_eff,
                     &rhor, &rhol, &pressr, &pressl,
                     var_datar, var_datal);

          if( validate(rhor_trusted[index], rhor, rhor_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable rho_r.\n"
                         "  rho_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", rhor_trusted[index], rhor, rhor_pert[index],
                                                     relative_error(rhor_trusted[index], rhor),
                                                     relative_error(rhor_trusted[index], rhor_pert[index]));
          if( validate(pressr_trusted[index], pressr, pressr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable press_r.\n"
                         "  press_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", pressr_trusted[index], pressr, pressr_pert[index],
                                                     relative_error(pressr_trusted[index], pressr),
                                                     relative_error(pressr_trusted[index], pressr_pert[index]));
          if( validate(vxr_trusted[index], var_datar[0], vxr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_r.\n"
                         "  vx_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxr_trusted[index], var_datar[0], vxr_pert[index],
                                                     relative_error(vxr_trusted[index], var_datar[0]),
                                                     relative_error(vxr_trusted[index], vxr_pert[index]));
          if( validate(vyr_trusted[index], var_datar[1], vyr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vy_r.\n"
                         "  vy_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vyr_trusted[index], var_datar[1], vyr_pert[index],
                                                     relative_error(vyr_trusted[index], var_datar[1]),
                                                     relative_error(vyr_trusted[index], vyr_pert[index]));
          if( validate(vzr_trusted[index], var_datar[2], vzr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vz_r.\n"
                         "  vz_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vzr_trusted[index], var_datar[2], vzr_pert[index],
                                                     relative_error(vzr_trusted[index], var_datar[2]),
                                                     relative_error(vzr_trusted[index], vzr_pert[index]));

          if( validate(rhol_trusted[index], rhol, rhol_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable rho_l.\n"
                         "  rho_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", rhol_trusted[index], rhol, rhol_pert[index],
                                                     relative_error(rhol_trusted[index], rhol),
                                                     relative_error(rhol_trusted[index], rhol_pert[index]));
          if( validate(pressl_trusted[index], pressl, pressl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable press_l.\n"
                         "  press_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", pressl_trusted[index], pressl, pressl_pert[index],
                                                     relative_error(pressl_trusted[index], pressl),
                                                     relative_error(pressl_trusted[index], pressl_pert[index]));
          if( validate(vxl_trusted[index], var_datal[0], vxl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_l.\n"
                         "  vx_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxl_trusted[index], var_datal[0], vxl_pert[index],
                                                     relative_error(vxl_trusted[index], var_datal[0]),
                                                     relative_error(vxl_trusted[index], vxl_pert[index]));
          if( validate(vyl_trusted[index], var_datal[1], vyl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vy_l.\n"
                         "  vy_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vyl_trusted[index], var_datal[1], vyl_pert[index],
                                                     relative_error(vyl_trusted[index], var_datal[1]),
                                                     relative_error(vyl_trusted[index], vyl_pert[index]));
          if( validate(vzl_trusted[index], var_datal[2], vzl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vz_l.\n"
                         "  vz_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vzl_trusted[index], var_datal[2], vzl_pert[index],
                                                     relative_error(vzl_trusted[index], var_datal[2]),
                                                     relative_error(vzl_trusted[index], vzl_pert[index]));
    }

    const int num_vars2 = 2;

    for(int k=kmin; k<kmax; k++)
      for(int j=jmin; j<jmax; j++)
        for(int i=imin; i<imax; i++) {
          const int index = indexf(dirlength,i,j,k);

          double press_stencil[6], v_flux_dir[6];
          double var_data[num_vars2][6], var_datar[num_vars2], var_datal[num_vars2];

          for(int ind=0; ind<6; ind++) {
            const int stencil  = indexf(dirlength, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3)); // PPM needs indices from -3 to +2
            v_flux_dir[ind]    = vflux[flux_dirn][stencil]; // Could be smaller; doesn't use full stencil
            press_stencil[ind] = press[stencil];
            var_data[0][ind]   = vx[stencil];
            var_data[1][ind]   = vz[stencil];
          }
          const double Gamma_eff = eos_Gamma_eff(&eos, rho[index], press[index]);
          ghl_ppm_no_rho_P(press_stencil, var_data,
                     num_vars2, v_flux_dir, Gamma_eff,
                     var_datar, var_datal);

          if( validate(vxr_trusted[index], var_datar[0], vxr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_r.\n"
                         "  vx_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxr_trusted[index], var_datar[0], vxr_pert[index],
                                                     relative_error(vxr_trusted[index], var_datar[0]),
                                                     relative_error(vxr_trusted[index], vxr_pert[index]));
          if( validate(vzr_trusted[index], var_datar[1], vzr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vz_r.\n"
                         "  vz_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vzr_trusted[index], var_datar[1], vzr_pert[index],
                                                     relative_error(vzr_trusted[index], var_datar[1]),
                                                     relative_error(vzr_trusted[index], vzr_pert[index]));

          if( validate(vxl_trusted[index], var_datal[0], vxl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_l.\n"
                         "  vx_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxl_trusted[index], var_datal[0], vxl_pert[index],
                                                     relative_error(vxl_trusted[index], var_datal[0]),
                                                     relative_error(vxl_trusted[index], vxl_pert[index]));
          if( validate(vzl_trusted[index], var_datal[1], vzl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vz_l.\n"
                         "  vz_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vzl_trusted[index], var_datal[1], vzl_pert[index],
                                                     relative_error(vzl_trusted[index], var_datal[1]),
                                                     relative_error(vzl_trusted[index], vzl_pert[index]));
    }
  }
  fclose(infile);
  fclose(inpert);
  ghl_info("ET_Legacy reconstruction test has passed!\n");
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

double eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in) {
  double K, Gamma;
  ghl_hybrid_get_K_and_Gamma(eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return eos->Gamma_th + (Gamma - eos->Gamma_th)*P_cold/press_in;
}
