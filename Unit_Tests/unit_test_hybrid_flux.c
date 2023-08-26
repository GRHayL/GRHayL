#include "unit_tests.h"

int main(int argc, char **argv) {

  // Set up test data
  FILE* infile = fopen_with_check("hybrid_flux_input.bin", "rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);

  const double poison = 1e300;

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        W_max,
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);


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
  double *rho_r   = (double*) malloc(sizeof(double)*arraylength);
  double *press_r = (double*) malloc(sizeof(double)*arraylength);
  double *vx_r    = (double*) malloc(sizeof(double)*arraylength);
  double *vy_r    = (double*) malloc(sizeof(double)*arraylength);
  double *vz_r    = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_r    = (double*) malloc(sizeof(double)*arraylength);
  double *By_r    = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_r    = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for left face
  double *rho_l   = (double*) malloc(sizeof(double)*arraylength);
  double *press_l = (double*) malloc(sizeof(double)*arraylength);
  double *vx_l    = (double*) malloc(sizeof(double)*arraylength);
  double *vy_l    = (double*) malloc(sizeof(double)*arraylength);
  double *vz_l    = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_l    = (double*) malloc(sizeof(double)*arraylength);
  double *By_l    = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_l    = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(lapse,   sizeof(double), arraylength, infile);
  key += fread(betax,   sizeof(double), arraylength, infile);
  key += fread(betay,   sizeof(double), arraylength, infile);
  key += fread(betaz,   sizeof(double), arraylength, infile);
  key += fread(gxx,     sizeof(double), arraylength, infile);
  key += fread(gxy,     sizeof(double), arraylength, infile);
  key += fread(gxz,     sizeof(double), arraylength, infile);
  key += fread(gyy,     sizeof(double), arraylength, infile);
  key += fread(gyz,     sizeof(double), arraylength, infile);
  key += fread(gzz,     sizeof(double), arraylength, infile);

  if(key != arraylength*10)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  key  = fread(rho_r,     sizeof(double), arraylength, infile);
  key += fread(press_r,   sizeof(double), arraylength, infile);
  key += fread(vx_r,      sizeof(double), arraylength, infile);
  key += fread(vy_r,      sizeof(double), arraylength, infile);
  key += fread(vz_r,      sizeof(double), arraylength, infile);
  key += fread(Bx_r,      sizeof(double), arraylength, infile);
  key += fread(By_r,      sizeof(double), arraylength, infile);
  key += fread(Bz_r,      sizeof(double), arraylength, infile);

  key += fread(rho_l,     sizeof(double), arraylength, infile);
  key += fread(press_l,   sizeof(double), arraylength, infile);
  key += fread(vx_l,      sizeof(double), arraylength, infile);
  key += fread(vy_l,      sizeof(double), arraylength, infile);
  key += fread(vz_l,      sizeof(double), arraylength, infile);
  key += fread(Bx_l,      sizeof(double), arraylength, infile);
  key += fread(By_l,      sizeof(double), arraylength, infile);
  key += fread(Bz_l,      sizeof(double), arraylength, infile);

  if(key != arraylength*16)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");
  fclose(infile);

  // Allocate memory for comparison data
  double *trusted_rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_tau_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_x_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_y_flux = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_z_flux = (double*) malloc(sizeof(double)*arraylength);

  double *pert_rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_tau_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_x_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_y_flux = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_z_flux = (double*) malloc(sizeof(double)*arraylength);

  FILE *outfile = fopen_with_check("hybrid_flux_output.bin", "rb");
  FILE *pertfile = fopen_with_check("hybrid_flux_output_pert.bin", "rb");

  // Function pointer to allow for loop over fluxes
  void (*calculate_HLLE_fluxes)(ghl_primitive_quantities *restrict, ghl_primitive_quantities *restrict,
                              const ghl_eos_parameters *restrict, const ghl_metric_quantities *restrict,
                              const double, const double, ghl_conservative_quantities *restrict);

  void (*calculate_characteristic_speed)(ghl_primitive_quantities *restrict, ghl_primitive_quantities *restrict,
                              const ghl_eos_parameters *restrict, const ghl_metric_quantities *restrict, double *restrict, double *restrict);

  for(int entropy=0; entropy<2; entropy++) {
    // Loop over flux directions (x,y,z)
    for(int flux_dirn=0; flux_dirn<3; flux_dirn++) {
      // Set function pointer to specific function for a given direction
      switch(flux_dirn) {
        case 0:
          calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
          calculate_HLLE_fluxes          = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn0_hybrid;
          break;
        case 1:
          calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
          calculate_HLLE_fluxes          = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn1_hybrid;
          break;
        case 2:
          calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
          calculate_HLLE_fluxes          = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn2_hybrid;
          break;
      }

      key  = fread(trusted_rho_star_flux, sizeof(double), arraylength, outfile);
      key += fread(trusted_tau_flux,      sizeof(double), arraylength, outfile);
      key += fread(trusted_S_x_flux,      sizeof(double), arraylength, outfile);
      key += fread(trusted_S_y_flux,      sizeof(double), arraylength, outfile);
      key += fread(trusted_S_z_flux,      sizeof(double), arraylength, outfile);
      
      if(key != arraylength*5)
        ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                     "is up-to-date with current test version.\n");
      
      key  = fread(pert_rho_star_flux, sizeof(double), arraylength, pertfile);
      key += fread(pert_tau_flux,      sizeof(double), arraylength, pertfile);
      key += fread(pert_S_x_flux,      sizeof(double), arraylength, pertfile);
      key += fread(pert_S_y_flux,      sizeof(double), arraylength, pertfile);
      key += fread(pert_S_z_flux,      sizeof(double), arraylength, pertfile);
      
      if(key != arraylength*5)
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
              rho_r[index], press_r[index], poison,
              vx_r[index], vy_r[index], vz_r[index],
              Bx_r[index], By_r[index], Bz_r[index],
              poison, poison, poison, // entropy, Y_e, temp
              &prims_r);

        ghl_initialize_primitives(
              rho_l[index], press_l[index], poison,
              vx_l[index], vy_l[index], vz_l[index],
              Bx_l[index], By_l[index], Bz_l[index],
              poison, poison, poison, // entropy, Y_e, temp
              &prims_l);

      int speed_limit __attribute__((unused)) = ghl_limit_v_and_compute_u0(
            &eos, &metric_face, &prims_r);
      speed_limit = ghl_limit_v_and_compute_u0(
            &eos, &metric_face, &prims_l);

        prims_r.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_r.rho, prims_r.press);
        prims_l.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_l.rho, prims_l.press);

        double cmin, cmax;
        calculate_characteristic_speed(
              &prims_r, &prims_l, &eos,
              &metric_face, &cmin, &cmax);

        ghl_conservative_quantities cons_fluxes;
        calculate_HLLE_fluxes(
              &prims_r, &prims_l, &eos,
              &metric_face, cmin, cmax,
              &cons_fluxes);

        if( ghl_pert_test_fail(trusted_rho_star_flux[index], cons_fluxes.rho, pert_rho_star_flux[index]) )
          ghl_error("Test unit_test_hybrid_flux has failed for variable rho_star_flux.\n"
                       "  rho_star_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_rho_star_flux[index], cons_fluxes.rho, pert_rho_star_flux[index],
                                                   relative_error(trusted_rho_star_flux[index], cons_fluxes.rho),
                                                   relative_error(trusted_rho_star_flux[index], pert_rho_star_flux[index]));
        if( ghl_pert_test_fail(trusted_tau_flux[index], cons_fluxes.tau, pert_tau_flux[index]) )
          ghl_error("Test unit_test_hybrid_flux has failed for variable tau_flux.\n"
                       "  tau_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_tau_flux[index], cons_fluxes.tau, pert_tau_flux[index],
                                                   relative_error(trusted_tau_flux[index], cons_fluxes.tau),
                                                   relative_error(trusted_tau_flux[index], pert_tau_flux[index]));
        if( ghl_pert_test_fail(trusted_S_x_flux[index], cons_fluxes.SD[0], pert_S_x_flux[index]) )
          ghl_error("Test unit_test_hybrid_flux has failed for variable S_x_flux.\n"
                       "  S_x_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_x_flux[index], cons_fluxes.SD[0], pert_S_x_flux[index],
                                                   relative_error(trusted_S_x_flux[index], cons_fluxes.SD[0]),
                                                   relative_error(trusted_S_x_flux[index], pert_S_x_flux[index]));
        if( ghl_pert_test_fail(trusted_S_y_flux[index], cons_fluxes.SD[1], pert_S_y_flux[index]) )
          ghl_error("Test unit_test_hybrid_flux has failed for variable S_y_flux.\n"
                       "  S_y_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_y_flux[index], cons_fluxes.SD[1], pert_S_y_flux[index],
                                                   relative_error(trusted_S_y_flux[index], cons_fluxes.SD[1]),
                                                   relative_error(trusted_S_y_flux[index], pert_S_y_flux[index]));
        if( ghl_pert_test_fail(trusted_S_z_flux[index], cons_fluxes.SD[2], pert_S_z_flux[index]) )
          ghl_error("Test unit_test_hybrid_flux has failed for variable S_z_flux.\n"
                       "  S_z_flux trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_z_flux[index], cons_fluxes.SD[2], pert_S_z_flux[index],
                                                   relative_error(trusted_S_z_flux[index], cons_fluxes.SD[2]),
                                                   relative_error(trusted_S_z_flux[index], pert_S_z_flux[index]));
      }
    } // flux_dir
  } // entropy

  fclose(outfile);
  fclose(pertfile);

  ghl_info("hybrid_flux test has passed!\n");
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(rho_r); free(press_r);
  free(vx_r); free(vy_r); free(vz_r);
  free(Bx_r); free(By_r); free(Bz_r);
  free(rho_l); free(press_l);
  free(vx_l); free(vy_l); free(vz_l);
  free(Bx_l); free(By_l); free(Bz_l);
  free(trusted_rho_star_flux); free(trusted_tau_flux);
  free(trusted_S_x_flux); free(trusted_S_y_flux); free(trusted_S_z_flux);
  free(pert_rho_star_flux); free(pert_tau_flux);
  free(pert_S_x_flux); free(pert_S_y_flux); free(pert_S_z_flux);
}
