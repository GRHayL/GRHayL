#include "ghl_unit_tests.h"
#include <limits.h>
#include <stdint.h>
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

static double eos_Gamma_eff(const ghl_eos_parameters *restrict eos, const double rho_in, const double press_in);

int main(int argc, char **argv) {
  const double dummy_density = 1.0;

  const ghl_con2prim_id_t None = ghl_con2prim_id_None;
  const ghl_con2prim_id_t backups[3] = {None, None, None};

  // None of these parameters actually matter. We are only using
  // the default values set in the initialize function for PPM.
  ghl_parameters params;
  ghl_initialize_params(None, backups, false, false, true, 0, 10, 0.0, &params);

  const int neos = 1;
  const double rho_ppoly_in[1] = {0.0};
  const double Gamma_ppoly_in[1] = {2.0};
  const double k_ppoly0 = 1.0;
  const double Gamma_th = 2.0;

  ghl_eos_parameters eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        dummy_density, dummy_density, dummy_density,
        neos, rho_ppoly_in, Gamma_ppoly_in,
        k_ppoly0, Gamma_th, &eos);

  FILE* infile = fopen_with_check("ET_Legacy_reconstruction_input.bin", "rb");

  int dirlength;
  if(fread(&dirlength, sizeof(int), 1, infile) != 1) {
    ghl_error(
          "An error has occurred with reading the direction length. Please check that "
          "data\n"
          "is up-to-date with current test version.\n");
  }
  if(dirlength < 8) {
    ghl_error(
          "The ET Legacy reconstruction data must have direction length at least 8.\n");
  }

  const size_t dir_count = (size_t)dirlength;
  if(dir_count > SIZE_MAX / dir_count) {
    ghl_error("The ET Legacy reconstruction dimensions are too large.\n");
  }
  const size_t slice_count = dir_count * dir_count;
  if(slice_count > SIZE_MAX / dir_count) {
    ghl_error("The ET Legacy reconstruction dimensions are too large.\n");
  }
  const size_t count = slice_count * dir_count;
  if(count > INT_MAX || count > SIZE_MAX / sizeof(double)) {
    ghl_error(
          "The ET Legacy reconstruction data is too large to index or allocate "
          "safely.\n");
  }

  double *rho = (double *)malloc(sizeof(double) * count);
  double *press = (double *)malloc(sizeof(double) * count);
  double *vx = (double *)malloc(sizeof(double) * count);
  double *vy = (double *)malloc(sizeof(double) * count);
  double *vz = (double *)malloc(sizeof(double) * count);

  double *rhor_trusted = (double *)malloc(sizeof(double) * count);
  double *rhol_trusted = (double *)malloc(sizeof(double) * count);
  double *pressr_trusted = (double *)malloc(sizeof(double) * count);
  double *pressl_trusted = (double *)malloc(sizeof(double) * count);
  double *vxr_trusted = (double *)malloc(sizeof(double) * count);
  double *vxl_trusted = (double *)malloc(sizeof(double) * count);
  double *vyr_trusted = (double *)malloc(sizeof(double) * count);
  double *vyl_trusted = (double *)malloc(sizeof(double) * count);
  double *vzr_trusted = (double *)malloc(sizeof(double) * count);
  double *vzl_trusted = (double *)malloc(sizeof(double) * count);

  double *rhor_pert = (double *)malloc(sizeof(double) * count);
  double *rhol_pert = (double *)malloc(sizeof(double) * count);
  double *pressr_pert = (double *)malloc(sizeof(double) * count);
  double *pressl_pert = (double *)malloc(sizeof(double) * count);
  double *vxr_pert = (double *)malloc(sizeof(double) * count);
  double *vxl_pert = (double *)malloc(sizeof(double) * count);
  double *vyr_pert = (double *)malloc(sizeof(double) * count);
  double *vyl_pert = (double *)malloc(sizeof(double) * count);
  double *vzr_pert = (double *)malloc(sizeof(double) * count);
  double *vzl_pert = (double *)malloc(sizeof(double) * count);

  if(rho == NULL || press == NULL || vx == NULL || vy == NULL || vz == NULL
     || rhor_trusted == NULL || rhol_trusted == NULL || pressr_trusted == NULL
     || pressl_trusted == NULL || vxr_trusted == NULL || vxl_trusted == NULL
     || vyr_trusted == NULL || vyl_trusted == NULL || vzr_trusted == NULL
     || vzl_trusted == NULL || rhor_pert == NULL || rhol_pert == NULL
     || pressr_pert == NULL || pressl_pert == NULL || vxr_pert == NULL
     || vxl_pert == NULL || vyr_pert == NULL || vyl_pert == NULL || vzr_pert == NULL
     || vzl_pert == NULL) {
    ghl_error("Failed to allocate ET Legacy reconstruction test data.\n");
  }

  const bool initial_data_read = fread(rho, sizeof(double), count, infile) == count
                                 && fread(press, sizeof(double), count, infile) == count
                                 && fread(vx, sizeof(double), count, infile) == count
                                 && fread(vy, sizeof(double), count, infile) == count
                                 && fread(vz, sizeof(double), count, infile) == count;

  fclose(infile);

  if(!initial_data_read) {
    ghl_error(
          "An error has occured with reading in initial data. Please check that data\n"
          "is up-to-date with current test version.\n");
  }

  infile = fopen_with_check("ET_Legacy_reconstruction_output.bin","rb");

  FILE* inpert = fopen_with_check("ET_Legacy_reconstruction_output_pert.bin","rb");

  for(int flux_dirn = 0; flux_dirn<3; flux_dirn++) {
    const int num_vars = 3;

    double *vflux[3] = {vx, vy, vz};
    const int xdir = (flux_dirn==0);
    const int ydir = (flux_dirn==1);
    const int zdir = (flux_dirn==2);

    const bool trusted_data_read
          = fread(rhor_trusted, sizeof(double), count, infile) == count
            && fread(rhol_trusted, sizeof(double), count, infile) == count
            && fread(pressr_trusted, sizeof(double), count, infile) == count
            && fread(pressl_trusted, sizeof(double), count, infile) == count
            && fread(vxr_trusted, sizeof(double), count, infile) == count
            && fread(vxl_trusted, sizeof(double), count, infile) == count
            && fread(vyr_trusted, sizeof(double), count, infile) == count
            && fread(vyl_trusted, sizeof(double), count, infile) == count
            && fread(vzr_trusted, sizeof(double), count, infile) == count
            && fread(vzl_trusted, sizeof(double), count, infile) == count;

    if(!trusted_data_read) {
      ghl_error(
            "An error has occured with reading in comparison data. Please check that "
            "data\n"
            "is up-to-date with current test version.\n");
    }

    const bool perturbed_data_read
          = fread(rhor_pert, sizeof(double), count, inpert) == count
            && fread(rhol_pert, sizeof(double), count, inpert) == count
            && fread(pressr_pert, sizeof(double), count, inpert) == count
            && fread(pressl_pert, sizeof(double), count, inpert) == count
            && fread(vxr_pert, sizeof(double), count, inpert) == count
            && fread(vxl_pert, sizeof(double), count, inpert) == count
            && fread(vyr_pert, sizeof(double), count, inpert) == count
            && fread(vyl_pert, sizeof(double), count, inpert) == count
            && fread(vzr_pert, sizeof(double), count, inpert) == count
            && fread(vzl_pert, sizeof(double), count, inpert) == count;

    if(!perturbed_data_read) {
      ghl_error(
            "An error has occured with reading in comparison data. Please check that "
            "data\n"
            "is up-to-date with current test version.\n");
    }

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

          double ftilde[2];
          ghl_compute_ftilde(&params, press_stencil, v_flux_dir, ftilde);

          const double Gamma_eff = eos_Gamma_eff(&eos, rho[index], press[index]);
          ghl_ppm_reconstruction_with_steepening(&params, press_stencil, Gamma_eff, ftilde, rho_stencil, &rhor, &rhol);

          ghl_ppm_reconstruction(ftilde, press_stencil, &pressr, &pressl);

          for(int ind=0; ind<num_vars; ind++)
            ghl_ppm_reconstruction(ftilde, var_data[ind], &var_datar[ind], &var_datal[ind]);

          if( ghl_pert_test_fail(rhor_trusted[index], rhor, rhor_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable rho_r.\n"
                         "  rho_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", rhor_trusted[index], rhor, rhor_pert[index],
                                                     relative_error(rhor_trusted[index], rhor),
                                                     relative_error(rhor_trusted[index], rhor_pert[index]));
          if( ghl_pert_test_fail(pressr_trusted[index], pressr, pressr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable press_r.\n"
                         "  press_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", pressr_trusted[index], pressr, pressr_pert[index],
                                                     relative_error(pressr_trusted[index], pressr),
                                                     relative_error(pressr_trusted[index], pressr_pert[index]));
          if( ghl_pert_test_fail(vxr_trusted[index], var_datar[0], vxr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_r.\n"
                         "  vx_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxr_trusted[index], var_datar[0], vxr_pert[index],
                                                     relative_error(vxr_trusted[index], var_datar[0]),
                                                     relative_error(vxr_trusted[index], vxr_pert[index]));
          if( ghl_pert_test_fail(vyr_trusted[index], var_datar[1], vyr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vy_r.\n"
                         "  vy_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vyr_trusted[index], var_datar[1], vyr_pert[index],
                                                     relative_error(vyr_trusted[index], var_datar[1]),
                                                     relative_error(vyr_trusted[index], vyr_pert[index]));
          if( ghl_pert_test_fail(vzr_trusted[index], var_datar[2], vzr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vz_r.\n"
                         "  vz_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vzr_trusted[index], var_datar[2], vzr_pert[index],
                                                     relative_error(vzr_trusted[index], var_datar[2]),
                                                     relative_error(vzr_trusted[index], vzr_pert[index]));

          if( ghl_pert_test_fail(rhol_trusted[index], rhol, rhol_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable rho_l.\n"
                         "  rho_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", rhol_trusted[index], rhol, rhol_pert[index],
                                                     relative_error(rhol_trusted[index], rhol),
                                                     relative_error(rhol_trusted[index], rhol_pert[index]));
          if( ghl_pert_test_fail(pressl_trusted[index], pressl, pressl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable press_l.\n"
                         "  press_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", pressl_trusted[index], pressl, pressl_pert[index],
                                                     relative_error(pressl_trusted[index], pressl),
                                                     relative_error(pressl_trusted[index], pressl_pert[index]));
          if( ghl_pert_test_fail(vxl_trusted[index], var_datal[0], vxl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_l.\n"
                         "  vx_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxl_trusted[index], var_datal[0], vxl_pert[index],
                                                     relative_error(vxl_trusted[index], var_datal[0]),
                                                     relative_error(vxl_trusted[index], vxl_pert[index]));
          if( ghl_pert_test_fail(vyl_trusted[index], var_datal[1], vyl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vy_l.\n"
                         "  vy_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vyl_trusted[index], var_datal[1], vyl_pert[index],
                                                     relative_error(vyl_trusted[index], var_datal[1]),
                                                     relative_error(vyl_trusted[index], vyl_pert[index]));
          if( ghl_pert_test_fail(vzl_trusted[index], var_datal[2], vzl_pert[index]) )
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
          double ftilde[2];
          ghl_compute_ftilde(&params, press_stencil, v_flux_dir, ftilde);

          for(int ind=0; ind<num_vars2; ind++)
            ghl_ppm_reconstruction(ftilde, var_data[ind], &var_datar[ind], &var_datal[ind]);

          if( ghl_pert_test_fail(vxr_trusted[index], var_datar[0], vxr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_r.\n"
                         "  vx_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxr_trusted[index], var_datar[0], vxr_pert[index],
                                                     relative_error(vxr_trusted[index], var_datar[0]),
                                                     relative_error(vxr_trusted[index], vxr_pert[index]));
          if( ghl_pert_test_fail(vzr_trusted[index], var_datar[1], vzr_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vz_r.\n"
                         "  vz_r trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vzr_trusted[index], var_datar[1], vzr_pert[index],
                                                     relative_error(vzr_trusted[index], var_datar[1]),
                                                     relative_error(vzr_trusted[index], vzr_pert[index]));

          if( ghl_pert_test_fail(vxl_trusted[index], var_datal[0], vxl_pert[index]) )
            ghl_error("Test unit_test_ET_Legacy_reconstruction has failed for variable vx_l.\n"
                         "  vx_l trusted %.14e computed %.14e perturbed %.14e\n"
                         "  rel.err. %.14e %.14e\n", vxl_trusted[index], var_datal[0], vxl_pert[index],
                                                     relative_error(vxl_trusted[index], var_datal[0]),
                                                     relative_error(vxl_trusted[index], vxl_pert[index]));
          if( ghl_pert_test_fail(vzl_trusted[index], var_datal[1], vzl_pert[index]) )
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

double eos_Gamma_eff(const ghl_eos_parameters *restrict eos, const double rho_in, const double press_in) {
  double K, Gamma;
  ghl_hybrid_get_K_and_Gamma(eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return eos->Gamma_th + (Gamma - eos->Gamma_th)*P_cold/press_in;
}
