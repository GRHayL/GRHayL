#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  // Set up test data
  FILE* infile = fopen_with_check("tabulated_flux_input.bin", "rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);

  const double poison = 1e300;

  const char tablepath[] = "LS220_234r_136t_50y_analmu_20091212_SVNr26.h5";
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Ye_min    = 5e-2;
  const double Ye_max    = 5e-1;
  const double T_min     = 1e-2;
  const double T_max     = 1e2;

  ghl_parameters params;
  params.max_Lorentz_factor = W_max;
  params.inv_sq_max_Lorentz_factor = 1.0/SQR(W_max);

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
        tablepath,
        rho_b_min, rho_b_min, rho_b_max,
        Ye_min, Ye_min, Ye_max,
        T_min, T_min, T_max,
        &eos);

  ghl_compute_h_and_cs2 = &ghl_test_compute_h_and_cs2;

  // Allocate memory for metric
  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);
  double *gxx   = (double*) malloc(sizeof(double)*arraylength);
  double *gxy   = (double*) malloc(sizeof(double)*arraylength);
  double *gxz   = (double*) malloc(sizeof(double)*arraylength);
  double *gyy   = (double*) malloc(sizeof(double)*arraylength);
  double *gyz   = (double*) malloc(sizeof(double)*arraylength);
  double *gzz   = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for right face
  double *rho_r = (double*) malloc(sizeof(double)*arraylength);
  double *Ye_r  = (double*) malloc(sizeof(double)*arraylength);
  double *T_r   = (double*) malloc(sizeof(double)*arraylength);
  double *vx_r  = (double*) malloc(sizeof(double)*arraylength);
  double *vy_r  = (double*) malloc(sizeof(double)*arraylength);
  double *vz_r  = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_r  = (double*) malloc(sizeof(double)*arraylength);
  double *By_r  = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_r  = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for left face
  double *rho_l = (double*) malloc(sizeof(double)*arraylength);
  double *Ye_l  = (double*) malloc(sizeof(double)*arraylength);
  double *T_l   = (double*) malloc(sizeof(double)*arraylength);
  double *vx_l  = (double*) malloc(sizeof(double)*arraylength);
  double *vy_l  = (double*) malloc(sizeof(double)*arraylength);
  double *vz_l  = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_l  = (double*) malloc(sizeof(double)*arraylength);
  double *By_l  = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_l  = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for characteristic speeds
  double *cxmin  = (double*) malloc(sizeof(double)*arraylength);
  double *cxmax  = (double*) malloc(sizeof(double)*arraylength);
  double *cymin  = (double*) malloc(sizeof(double)*arraylength);
  double *cymax  = (double*) malloc(sizeof(double)*arraylength);
  double *czmin  = (double*) malloc(sizeof(double)*arraylength);
  double *czmax  = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(lapse, sizeof(double), arraylength, infile);
  key += fread(betax, sizeof(double), arraylength, infile);
  key += fread(betay, sizeof(double), arraylength, infile);
  key += fread(betaz, sizeof(double), arraylength, infile);
  key += fread(gxx,   sizeof(double), arraylength, infile);
  key += fread(gxy,   sizeof(double), arraylength, infile);
  key += fread(gxz,   sizeof(double), arraylength, infile);
  key += fread(gyy,   sizeof(double), arraylength, infile);
  key += fread(gyz,   sizeof(double), arraylength, infile);
  key += fread(gzz,   sizeof(double), arraylength, infile);

  if(key != arraylength*10)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  key  = fread(rho_r, sizeof(double), arraylength, infile);
  key += fread(Ye_r,  sizeof(double), arraylength, infile);
  key += fread(T_r,   sizeof(double), arraylength, infile);
  key += fread(vx_r,  sizeof(double), arraylength, infile);
  key += fread(vy_r,  sizeof(double), arraylength, infile);
  key += fread(vz_r,  sizeof(double), arraylength, infile);
  key += fread(Bx_r,  sizeof(double), arraylength, infile);
  key += fread(By_r,  sizeof(double), arraylength, infile);
  key += fread(Bz_r,  sizeof(double), arraylength, infile);

  key += fread(rho_l, sizeof(double), arraylength, infile);
  key += fread(Ye_l,  sizeof(double), arraylength, infile);
  key += fread(T_l,   sizeof(double), arraylength, infile);
  key += fread(vx_l,  sizeof(double), arraylength, infile);
  key += fread(vy_l,  sizeof(double), arraylength, infile);
  key += fread(vz_l,  sizeof(double), arraylength, infile);
  key += fread(Bx_l,  sizeof(double), arraylength, infile);
  key += fread(By_l,  sizeof(double), arraylength, infile);
  key += fread(Bz_l,  sizeof(double), arraylength, infile);

  key += fread(cxmin, sizeof(double), arraylength, infile);
  key += fread(cxmax, sizeof(double), arraylength, infile);
  key += fread(cymin, sizeof(double), arraylength, infile);
  key += fread(cymax, sizeof(double), arraylength, infile);
  key += fread(czmin, sizeof(double), arraylength, infile);
  key += fread(czmax, sizeof(double), arraylength, infile);

  if(key != arraylength*24)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");
  fclose(infile);

  // Allocate memory for comparison data
  double *trusted_rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_Y_e_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_tau_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_x_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_y_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_z_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_ent_flux = (double*) malloc(sizeof(double)*arraylength);

  double *pert_rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_Y_e_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_tau_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_x_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_y_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_z_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_ent_flux = (double*) malloc(sizeof(double)*arraylength);

  FILE *outfile = fopen_with_check("tabulated_flux_output.bin", "rb");
  FILE *pertfile = fopen_with_check("tabulated_flux_output_pert.bin", "rb");

  double *cmin;
  double *cmax;
  // Function pointer to allow for loop over fluxes
  void (*calculate_HLLE_fluxes)(
        const int direction,
        const ghl_eos_parameters *restrict,
        const ghl_metric_quantities *restrict,
        ghl_primitive_quantities *restrict,
        ghl_primitive_quantities *restrict,
        const double, const double,
        ghl_conservative_quantities *restrict);

  for(int entropy=0; entropy<2; entropy++) {
    calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_tabulated_entropy : &ghl_calculate_HLLE_fluxes_tabulated;
    // Loop over flux directions (x,y,z)
    for(int flux_dirn=0; flux_dirn<3; flux_dirn++) {
      // Set function pointer to specific function for a given direction
      switch(flux_dirn) {
        case 0:
          cmin = cxmin;
          cmax = cxmax;
          break;
        case 1:
          cmin = cymin;
          cmax = cymax;
          break;
        case 2:
          cmin = czmin;
          cmax = czmax;
          break;
      }

      key  = fread(trusted_rho_star_flux, sizeof(double), arraylength, outfile);
      key += fread(trusted_Y_e_flux,      sizeof(double), arraylength, outfile);
      key += fread(trusted_tau_flux,      sizeof(double), arraylength, outfile);
      key += fread(trusted_S_x_flux,      sizeof(double), arraylength, outfile);
      key += fread(trusted_S_y_flux,      sizeof(double), arraylength, outfile);
      key += fread(trusted_S_z_flux,      sizeof(double), arraylength, outfile);
      if(entropy)
        key += fread(trusted_ent_flux,    sizeof(double), arraylength, outfile);

      if(key != arraylength*(6+entropy))
        ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                     "is up-to-date with current test version.\n");

      key  = fread(pert_rho_star_flux, sizeof(double), arraylength, pertfile);
      key += fread(pert_Y_e_flux,      sizeof(double), arraylength, pertfile);
      key += fread(pert_tau_flux,      sizeof(double), arraylength, pertfile);
      key += fread(pert_S_x_flux,      sizeof(double), arraylength, pertfile);
      key += fread(pert_S_y_flux,      sizeof(double), arraylength, pertfile);
      key += fread(pert_S_z_flux,      sizeof(double), arraylength, pertfile);
      if(entropy)
        key += fread(pert_ent_flux,    sizeof(double), arraylength, pertfile);

      if(key != arraylength*(6+entropy))
        ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                     "is up-to-date with current test version.\n");

      for(int index=0; index<arraylength; index++) {

        ghl_metric_quantities metric_face;
        ghl_initialize_metric(
              lapse[index], betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &metric_face);

        ghl_primitive_quantities prims_r, prims_l;
        ghl_initialize_primitives(
              rho_r[index], poison, poison,
              vx_r[index], vy_r[index], vz_r[index],
              Bx_r[index], By_r[index], Bz_r[index],
              poison, Ye_r[index], T_r[index],
              &prims_r);
        ghl_tabulated_compute_P_eps_S_from_T(
              &eos,
              prims_r.rho, prims_r.Y_e, prims_r.temperature,
              &prims_r.press, &prims_r.eps, &prims_r.entropy);

        ghl_initialize_primitives(
              rho_l[index], poison, poison,
              vx_l[index], vy_l[index], vz_l[index],
              Bx_l[index], By_l[index], Bz_l[index],
              poison, Ye_l[index], T_l[index],
              &prims_l);
        ghl_tabulated_compute_P_eps_S_from_T(
              &eos,
              prims_l.rho, prims_l.Y_e, prims_l.temperature,
              &prims_l.press, &prims_l.eps, &prims_l.entropy);

        bool speed_limited;
        ghl_error_codes_t error = ghl_limit_v_and_compute_u0(
              &params, &metric_face, &prims_r, &speed_limited);
        if(error)
          ghl_read_error_codes(error);
        error = ghl_limit_v_and_compute_u0(
              &params, &metric_face, &prims_l, &speed_limited);
        if(error)
          ghl_read_error_codes(error);

        ghl_conservative_quantities cons_fluxes;
        calculate_HLLE_fluxes(
              flux_dirn, &eos, &metric_face,
              &prims_r, &prims_l, cmin[index], cmax[index],
              &cons_fluxes);

        if( ghl_pert_test_fail(trusted_rho_star_flux[index], cons_fluxes.rho, pert_rho_star_flux[index]) )
          ghl_error("Test unit_test_tabulated_flux has failed for variable rho_star_flux.\n"
                       "  rho_star_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_rho_star_flux[index], cons_fluxes.rho, pert_rho_star_flux[index],
                                                   relative_error(trusted_rho_star_flux[index], cons_fluxes.rho),
                                                   relative_error(trusted_rho_star_flux[index], pert_rho_star_flux[index]));

        if( ghl_pert_test_fail(trusted_Y_e_flux[index], cons_fluxes.Y_e, pert_Y_e_flux[index]) )
          ghl_error("Test unit_test_tabulated_flux has failed for variable Y_e_flux.\n"
                       "  Y_e_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_Y_e_flux[index], cons_fluxes.Y_e, pert_Y_e_flux[index],
                                                   relative_error(trusted_Y_e_flux[index], cons_fluxes.Y_e),
                                                   relative_error(trusted_Y_e_flux[index], pert_Y_e_flux[index]));

        if( ghl_pert_test_fail(trusted_tau_flux[index], cons_fluxes.tau, pert_tau_flux[index]) )
          ghl_error("Test unit_test_tabulated_flux has failed for variable tau_flux.\n"
                       "  tau_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_tau_flux[index], cons_fluxes.tau, pert_tau_flux[index],
                                                   relative_error(trusted_tau_flux[index], cons_fluxes.tau),
                                                   relative_error(trusted_tau_flux[index], pert_tau_flux[index]));

        if( ghl_pert_test_fail(trusted_S_x_flux[index], cons_fluxes.SD[0], pert_S_x_flux[index]) )
          ghl_error("Test unit_test_tabulated_flux has failed for variable S_x_flux.\n"
                       "  S_x_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_x_flux[index], cons_fluxes.SD[0], pert_S_x_flux[index],
                                                   relative_error(trusted_S_x_flux[index], cons_fluxes.SD[0]),
                                                   relative_error(trusted_S_x_flux[index], pert_S_x_flux[index]));

        if( ghl_pert_test_fail(trusted_S_y_flux[index], cons_fluxes.SD[1], pert_S_y_flux[index]) )
          ghl_error("Test unit_test_tabulated_flux has failed for variable S_y_flux.\n"
                       "  S_y_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_y_flux[index], cons_fluxes.SD[1], pert_S_y_flux[index],
                                                   relative_error(trusted_S_y_flux[index], cons_fluxes.SD[1]),
                                                   relative_error(trusted_S_y_flux[index], pert_S_y_flux[index]));

        if( ghl_pert_test_fail(trusted_S_z_flux[index], cons_fluxes.SD[2], pert_S_z_flux[index]) )
          ghl_error("Test unit_test_tabulated_flux has failed for variable S_z_flux.\n"
                       "  S_z_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_z_flux[index], cons_fluxes.SD[2], pert_S_z_flux[index],
                                                   relative_error(trusted_S_z_flux[index], cons_fluxes.SD[2]),
                                                   relative_error(trusted_S_z_flux[index], pert_S_z_flux[index]));

        if( entropy && ghl_pert_test_fail(trusted_ent_flux[index], cons_fluxes.entropy, pert_ent_flux[index]) )
          ghl_error("Test unit_test_tabulated_flux has failed for variable ent_flux.\n"
                       "  ent_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_ent_flux[index], cons_fluxes.entropy, pert_ent_flux[index],
                                                   relative_error(trusted_ent_flux[index], cons_fluxes.entropy),
                                                   relative_error(trusted_ent_flux[index], pert_ent_flux[index]));
      }
    } // flux_dir
  } // entropy

  fclose(outfile);
  fclose(pertfile);

  ghl_info("tabulated_flux test has passed!\n");
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(rho_r); free(Ye_r); free(T_r);
  free(vx_r); free(vy_r); free(vz_r);
  free(Bx_r); free(By_r); free(Bz_r);
  free(rho_l); free(Ye_l); free(T_l);
  free(vx_l); free(vy_l); free(vz_l);
  free(Bx_l); free(By_l); free(Bz_l);
  free(cxmin); free(cymin); free(czmin);
  free(cxmax); free(cymax); free(czmax);
  free(trusted_rho_star_flux); free(trusted_tau_flux); free(trusted_Y_e_flux);
  free(trusted_S_x_flux); free(trusted_S_y_flux); free(trusted_S_z_flux);
  free(pert_rho_star_flux); free(pert_tau_flux); free(pert_Y_e_flux);
  free(pert_S_x_flux); free(pert_S_y_flux); free(pert_S_z_flux);
  free(trusted_ent_flux); free(pert_ent_flux);
}
