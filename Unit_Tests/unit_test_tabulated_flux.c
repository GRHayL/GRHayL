#include "ghl_unit_tests.h"

typedef ghl_error_codes_t (*flux_function)(
      ghl_primitive_quantities *restrict,
      ghl_primitive_quantities *restrict,
      const ghl_eos_parameters *restrict,
      const ghl_metric_quantities *restrict,
      double,
      double,
      ghl_conservative_quantities *restrict);

typedef void (*legacy_flux_function)(
      ghl_primitive_quantities *restrict,
      ghl_primitive_quantities *restrict,
      const ghl_eos_parameters *restrict,
      const ghl_metric_quantities *restrict,
      double,
      double,
      ghl_conservative_quantities *restrict);

static ghl_error_codes_t (*saved_compute_h_and_cs2)(
      const ghl_eos_parameters *restrict,
      ghl_primitive_quantities *restrict,
      double *restrict,
      double *restrict);
static int callback_calls;
static int fail_on_callback;

static ghl_error_codes_t injected_h_and_cs2_failure(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      double *restrict h,
      double *restrict cs2) {
  callback_calls++;
  if(callback_calls == fail_on_callback) {
    return ghl_error_table_max_T;
  }
  return saved_compute_h_and_cs2(eos, prims, h, cs2);
}

static bool conservative_is_poisoned(const ghl_conservative_quantities *restrict cons) {
  return cons->rho == 91.0 && cons->tau == 91.0 && cons->Y_e == 91.0
         && cons->SD[0] == 91.0 && cons->SD[1] == 91.0 && cons->SD[2] == 91.0
         && cons->entropy == 91.0;
}

static bool tabulated_fluxes_differ(
      const ghl_conservative_quantities *restrict a,
      const ghl_conservative_quantities *restrict b,
      const bool entropy) {
  return a->rho != b->rho || a->tau != b->tau || a->Y_e != b->Y_e || a->SD[0] != b->SD[0]
         || a->SD[1] != b->SD[1] || a->SD[2] != b->SD[2]
         || (entropy && a->entropy != b->entropy);
}

static bool flux_value_mismatch(const double expected, const double actual) {
  const double rtol = 8.0e-14;
  const double atol = 1.0e-30;
  return !isfinite(expected) || !isfinite(actual)
         || (fabs(expected - actual) > atol
             && fabs(expected - actual) > rtol * fmax(fabs(expected), fabs(actual)));
}

static bool legacy_floor_mismatch(const double expected, const double actual) {
  // Compatibility envelope for trusted outputs generated before wave-speed flooring.
  const double rtol = 5.0e-7;
  const double atol = 16.0 * DBL_EPSILON;
  return !isfinite(expected) || !isfinite(actual)
         || (fabs(expected - actual) > atol
             && fabs(expected - actual) > rtol * fmax(fabs(expected), fabs(actual)));
}

static bool flux_fixture_mismatch(
      const double expected,
      const double actual,
      const double perturbed,
      const double cmin,
      const double cmax) {
  if(!isfinite(expected) || !isfinite(actual) || !isfinite(perturbed)) {
    return true;
  }
  // The pinned fixtures predate flooring roundoff-negative wave speeds at zero.
  if(cmin < 0.0 || cmax < 0.0) {
    return legacy_floor_mismatch(expected, actual);
  }
  return ghl_pert_test_fail(expected, actual, perturbed);
}

static void check_tabulated_flux_contract(const ghl_eos_parameters *restrict eos) {
  const flux_function functions[2][3]
        = { { ghl_calculate_HLLE_fluxes_dirn0_tabulated_checked,
              ghl_calculate_HLLE_fluxes_dirn1_tabulated_checked,
              ghl_calculate_HLLE_fluxes_dirn2_tabulated_checked },
            { ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy_checked,
              ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy_checked,
              ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy_checked } };
  const legacy_flux_function legacy_functions[2][3]
        = { { ghl_calculate_HLLE_fluxes_dirn0_tabulated,
              ghl_calculate_HLLE_fluxes_dirn1_tabulated,
              ghl_calculate_HLLE_fluxes_dirn2_tabulated },
            { ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy,
              ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy,
              ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy } };
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims_r, prims_l;
  ghl_initialize_primitives(
        2e-4, 0.05, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.3, 1.0, &prims_r);
  prims_r.u0 = 1.0;
  ghl_initialize_primitives(
        1e-4, 0.02, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.2, 1.0, &prims_l);
  prims_l.u0 = 1.0;

  for(int entropy = 0; entropy < 2; entropy++) {
    for(int dir = 0; dir < 3; dir++) {
      ghl_conservative_quantities cons
            = { 91.0, 91.0, 91.0, { 91.0, 91.0, 91.0 }, 91.0 };
      if(functions[entropy][dir](&prims_r, &prims_l, eos, &metric, 0.0, 0.0, &cons)
               != ghl_error_invalid_hlle_wavespeeds
         || !conservative_is_poisoned(&cons)) {
        ghl_error(
              "HLLE zero-speed contract failed for tabulated entropy=%d direction=%d\n",
              entropy, dir);
      }
      if(functions[entropy][dir](&prims_r, &prims_l, eos, &metric, INFINITY, 1.0, &cons)
               != ghl_error_invalid_hlle_wavespeeds
         || !conservative_is_poisoned(&cons)) {
        ghl_error(
              "HLLE nonfinite-speed contract failed for tabulated entropy=%d "
              "direction=%d\n",
              entropy, dir);
      }
      if(functions[entropy][dir](
               &prims_r, &prims_l, eos, &metric, 0.5 / DBL_MAX, 0.0, &cons)
               != ghl_error_invalid_hlle_wavespeeds
         || !conservative_is_poisoned(&cons)) {
        ghl_error(
              "HLLE tiny-speed contract failed for tabulated entropy=%d direction=%d\n",
              entropy, dir);
      }
      if(functions[entropy][dir](&prims_r, &prims_l, eos, &metric, 1e200, 1e200, &cons)
               != ghl_error_invalid_hlle_wavespeeds
         || !conservative_is_poisoned(&cons)) {
        ghl_error(
              "HLLE product-overflow contract failed for tabulated entropy=%d "
              "direction=%d\n",
              entropy, dir);
      }
      if(functions[entropy][dir](&prims_r, &prims_l, eos, &metric, 0.0, 1.0, &cons)
         != ghl_success) {
        ghl_error(
              "HLLE one-sided speed failed for tabulated entropy=%d direction=%d\n",
              entropy, dir);
      }
      legacy_functions[entropy][dir](&prims_r, &prims_l, eos, &metric, 0.0, 1.0, &cons);
      ghl_conservative_quantities zero_floored = { 0 };
      ghl_conservative_quantities residue_floored = { 0 };
      if(functions[entropy][dir](
               &prims_r, &prims_l, eos, &metric, 100.0, 0.0, &zero_floored)
               != ghl_success
         || functions[entropy][dir](
                  &prims_r, &prims_l, eos, &metric, 100.0, -8.0 * DBL_EPSILON,
                  &residue_floored)
                  != ghl_success
         || tabulated_fluxes_differ(&zero_floored, &residue_floored, entropy)) {
        ghl_error(
              "HLLE roundoff-negative floor failed for tabulated entropy=%d "
              "direction=%d\n",
              entropy, dir);
      }
      if(functions[entropy][dir](
               &prims_r, &prims_l, eos, &metric, 0.0, 100.0, &zero_floored)
               != ghl_success
         || functions[entropy][dir](
                  &prims_r, &prims_l, eos, &metric, -8.0 * DBL_EPSILON, 100.0,
                  &residue_floored)
                  != ghl_success
         || tabulated_fluxes_differ(&zero_floored, &residue_floored, entropy)) {
        ghl_error(
              "HLLE roundoff-negative floor failed for tabulated entropy=%d "
              "direction=%d\n",
              entropy, dir);
      }
      if(functions[entropy][dir](&prims_r, &prims_l, eos, &metric, 0.4, 0.7, &cons)
         != ghl_success) {
        ghl_error(
              "HLLE fixed-bound call failed for tabulated entropy=%d direction=%d\n",
              entropy, dir);
      }
      const double invsum = 1.0 / 1.1;
      if(flux_value_mismatch(invsum * (-0.28e-4), cons.rho)
         || flux_value_mismatch(invsum * (-0.28 * (1.2e-4 - 0.3e-4)), cons.tau)
         || flux_value_mismatch(invsum * (-0.28 * (0.6e-4 - 0.2e-4)), cons.Y_e)
         || flux_value_mismatch(invsum * (0.4 * 0.05 + 0.7 * 0.02), cons.SD[dir])
         || flux_value_mismatch(0.0, cons.SD[(dir + 1) % 3])
         || flux_value_mismatch(0.0, cons.SD[(dir + 2) % 3])) {
        ghl_error(
              "Independent asymmetric HLLE check failed for tabulated entropy=%d "
              "direction=%d\n",
              entropy, dir);
      }
      if(entropy && flux_value_mismatch(invsum * (-0.28 * (0.7 - 0.4)), cons.entropy)) {
        ghl_error(
              "Independent entropy HLLE check failed for tabulated direction=%d\n", dir);
      }
      if(functions[entropy][dir](
               &prims_r, &prims_l, eos, &metric, 0.0, 2.0 / DBL_MAX, &cons)
         != ghl_success) {
        ghl_error(
              "HLLE invertible-small speed failed for tabulated entropy=%d "
              "direction=%d\n",
              entropy, dir);
      }

      cons = (ghl_conservative_quantities){
        91.0, 91.0, 91.0, { 91.0, 91.0, 91.0 }, 91.0
      };
      saved_compute_h_and_cs2 = ghl_compute_h_and_cs2;
      ghl_compute_h_and_cs2 = injected_h_and_cs2_failure;
      for(fail_on_callback = 1; fail_on_callback <= 2; fail_on_callback++) {
        callback_calls = 0;
        const ghl_error_codes_t error = functions[entropy][dir](
              &prims_r, &prims_l, eos, &metric, 1.0, 1.0, &cons);
        if(error != ghl_error_table_max_T || !conservative_is_poisoned(&cons)) {
          ghl_error(
                "HLLE callback-error contract failed for tabulated entropy=%d "
                "direction=%d callback=%d\n",
                entropy, dir, fail_on_callback);
        }
      }
      ghl_compute_h_and_cs2 = saved_compute_h_and_cs2;
    }
  }
}

int main(int argc, char **argv) {

  // Set up test data
  FILE* infile = fopen_with_check("tabulated_flux_input.bin", "rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);
  if(key != 1 || arraylength != 10000) {
    ghl_error("Invalid tabulated_flux_input.bin length (expected 10000)\n");
  }

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

  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.clean_sound_speed = true;
  ghl_error_codes_t error = ghl_initialize_tabulated_eos_functions_and_params(
        tablepath,
        rho_b_min, rho_b_min, rho_b_max,
        Ye_min, Ye_min, Ye_max,
        T_min, T_min, T_max,
        &eos);
  ghl_abort_if_error(error);

  ghl_compute_h_and_cs2 = &ghl_test_compute_h_and_cs2;
  check_tabulated_flux_contract(&eos);

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
  ghl_error_codes_t (*calculate_HLLE_fluxes)(
        ghl_primitive_quantities *restrict, ghl_primitive_quantities *restrict,
        const ghl_eos_parameters *restrict, const ghl_metric_quantities *restrict,
        const double, const double, ghl_conservative_quantities *restrict);

  for(int entropy=0; entropy<2; entropy++) {
    // Loop over flux directions (x,y,z)
    for(int flux_dirn=0; flux_dirn<3; flux_dirn++) {
      // Set function pointer to specific function for a given direction
      switch(flux_dirn) {
        case 0:
          cmin = cxmin;
          cmax = cxmax;
          calculate_HLLE_fluxes
                = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy_checked
                            : &ghl_calculate_HLLE_fluxes_dirn0_tabulated_checked;
          break;
        case 1:
          cmin = cymin;
          cmax = cymax;
          calculate_HLLE_fluxes
                = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy_checked
                            : &ghl_calculate_HLLE_fluxes_dirn1_tabulated_checked;
          break;
        case 2:
          cmin = czmin;
          cmax = czmax;
          calculate_HLLE_fluxes
                = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy_checked
                            : &ghl_calculate_HLLE_fluxes_dirn2_tabulated_checked;
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
        error = ghl_tabulated_compute_P_eps_S_from_T(
              &eos, prims_r.rho, prims_r.Y_e, prims_r.temperature, &prims_r.press,
              &prims_r.eps, &prims_r.entropy);
        if(error != ghl_success) {
          ghl_error(
                "Tabulated flux replay failed for right state %d with status %d.\n",
                index, error);
        }

        ghl_initialize_primitives(
              rho_l[index], poison, poison,
              vx_l[index], vy_l[index], vz_l[index],
              Bx_l[index], By_l[index], Bz_l[index],
              poison, Ye_l[index], T_l[index],
              &prims_l);
        error = ghl_tabulated_compute_P_eps_S_from_T(
              &eos, prims_l.rho, prims_l.Y_e, prims_l.temperature, &prims_l.press,
              &prims_l.eps, &prims_l.entropy);
        if(error != ghl_success) {
          ghl_error(
                "Tabulated flux replay failed for left state %d with status %d.\n",
                index, error);
        }

        bool speed_limited = false;
        ghl_error_codes_t error = ghl_limit_v_and_compute_u0(
              &params, &metric_face, &prims_r, &speed_limited);
        ghl_abort_if_error(error);
        error = ghl_limit_v_and_compute_u0(
              &params, &metric_face, &prims_l, &speed_limited);
        ghl_abort_if_error(error);

        ghl_conservative_quantities cons_fluxes;
        error = calculate_HLLE_fluxes(
              &prims_r, &prims_l, &eos, &metric_face, cmin[index], cmax[index],
              &cons_fluxes);
        ghl_abort_if_error(error);

        if(flux_fixture_mismatch(
                 trusted_rho_star_flux[index], cons_fluxes.rho,
                 pert_rho_star_flux[index], cmin[index], cmax[index])) {
          ghl_error(
                "Test unit_test_tabulated_flux has failed for variable rho_star_flux.\n"
                "  rho_star_flux trusted %.14e computed %.14e perturbed %.14e\n"
                "  rel.err. %.14e %.14e\n",
                trusted_rho_star_flux[index], cons_fluxes.rho, pert_rho_star_flux[index],
                relative_error(trusted_rho_star_flux[index], cons_fluxes.rho),
                relative_error(trusted_rho_star_flux[index], pert_rho_star_flux[index]));
        }

        if(flux_fixture_mismatch(
                 trusted_Y_e_flux[index], cons_fluxes.Y_e, pert_Y_e_flux[index],
                 cmin[index], cmax[index])) {
          ghl_error(
                "Test unit_test_tabulated_flux has failed for variable Y_e_flux.\n"
                "  Y_e_flux trusted %.14e computed %.14e perturbed %.14e\n"
                "  rel.err. %.14e %.14e\n",
                trusted_Y_e_flux[index], cons_fluxes.Y_e, pert_Y_e_flux[index],
                relative_error(trusted_Y_e_flux[index], cons_fluxes.Y_e),
                relative_error(trusted_Y_e_flux[index], pert_Y_e_flux[index]));
        }

        if(flux_fixture_mismatch(
                 trusted_tau_flux[index], cons_fluxes.tau, pert_tau_flux[index],
                 cmin[index], cmax[index])) {
          ghl_error(
                "Test unit_test_tabulated_flux has failed for variable tau_flux.\n"
                "  tau_flux trusted %.14e computed %.14e perturbed %.14e\n"
                "  rel.err. %.14e %.14e\n",
                trusted_tau_flux[index], cons_fluxes.tau, pert_tau_flux[index],
                relative_error(trusted_tau_flux[index], cons_fluxes.tau),
                relative_error(trusted_tau_flux[index], pert_tau_flux[index]));
        }

        if(flux_fixture_mismatch(
                 trusted_S_x_flux[index], cons_fluxes.SD[0], pert_S_x_flux[index],
                 cmin[index], cmax[index])) {
          ghl_error(
                "Test unit_test_tabulated_flux has failed for variable S_x_flux.\n"
                "  S_x_flux trusted %.14e computed %.14e perturbed %.14e\n"
                "  rel.err. %.14e %.14e\n",
                trusted_S_x_flux[index], cons_fluxes.SD[0], pert_S_x_flux[index],
                relative_error(trusted_S_x_flux[index], cons_fluxes.SD[0]),
                relative_error(trusted_S_x_flux[index], pert_S_x_flux[index]));
        }

        if(flux_fixture_mismatch(
                 trusted_S_y_flux[index], cons_fluxes.SD[1], pert_S_y_flux[index],
                 cmin[index], cmax[index])) {
          ghl_error(
                "Test unit_test_tabulated_flux has failed for variable S_y_flux.\n"
                "  S_y_flux trusted %.14e computed %.14e perturbed %.14e\n"
                "  rel.err. %.14e %.14e\n",
                trusted_S_y_flux[index], cons_fluxes.SD[1], pert_S_y_flux[index],
                relative_error(trusted_S_y_flux[index], cons_fluxes.SD[1]),
                relative_error(trusted_S_y_flux[index], pert_S_y_flux[index]));
        }

        if(flux_fixture_mismatch(
                 trusted_S_z_flux[index], cons_fluxes.SD[2], pert_S_z_flux[index],
                 cmin[index], cmax[index])) {
          ghl_error(
                "Test unit_test_tabulated_flux has failed for variable S_z_flux.\n"
                "  S_z_flux trusted %.14e computed %.14e perturbed %.14e\n"
                "  rel.err. %.14e %.14e\n",
                trusted_S_z_flux[index], cons_fluxes.SD[2], pert_S_z_flux[index],
                relative_error(trusted_S_z_flux[index], cons_fluxes.SD[2]),
                relative_error(trusted_S_z_flux[index], pert_S_z_flux[index]));
        }

        if(entropy
           && flux_fixture_mismatch(
                 trusted_ent_flux[index], cons_fluxes.entropy, pert_ent_flux[index],
                 cmin[index], cmax[index])) {
          ghl_error(
                "Test unit_test_tabulated_flux has failed for variable ent_flux.\n"
                "  ent_flux trusted %.14e computed %.14e perturbed %.14e\n"
                "  rel.err. %.14e %.14e\n",
                trusted_ent_flux[index], cons_fluxes.entropy, pert_ent_flux[index],
                relative_error(trusted_ent_flux[index], cons_fluxes.entropy),
                relative_error(trusted_ent_flux[index], pert_ent_flux[index]));
        }
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
