#include "ghl_nrpyeos_tabulated.h"
#include "ghl_unit_tests.h"

#include <float.h>
#include <hdf5.h>
#include <stdint.h>

typedef enum {
  compose_direct,
  compose_log_energy,
  compose_log_pressure,
  compose_cs2,
  compose_dedt,
  compose_dpderho,
  compose_dpdrhoe,
} compose_conversion;

typedef struct {
  const char *name;
  NRPyEOS_keys key;
  compose_conversion conversion;
} compose_field;

static const compose_field fields[] = {
  { "Abar", NRPyEOS_Abar_key, compose_direct },
  { "Xa", NRPyEOS_X_a_key, compose_direct },
  { "Xh", NRPyEOS_X_h_key, compose_direct },
  { "Xn", NRPyEOS_X_n_key, compose_direct },
  { "Xp", NRPyEOS_X_p_key, compose_direct },
  { "Zbar", NRPyEOS_Zbar_key, compose_direct },
  { "cs2", NRPyEOS_cs2_key, compose_cs2 },
  { "dedt", NRPyEOS_depsdT_key, compose_dedt },
  { "dpderho", NRPyEOS_dPdeps_key, compose_dpderho },
  { "dpdrhoe", NRPyEOS_dPdrho_key, compose_dpdrhoe },
  { "entropy", NRPyEOS_entropy_key, compose_direct },
  { "gamma", NRPyEOS_Gamma_key, compose_direct },
  { "logenergy", NRPyEOS_eps_key, compose_log_energy },
  { "logpress", NRPyEOS_press_key, compose_log_pressure },
  { "mu_e", NRPyEOS_mu_e_key, compose_direct },
  { "mu_n", NRPyEOS_mu_n_key, compose_direct },
  { "mu_p", NRPyEOS_mu_p_key, compose_direct },
  { "muhat", NRPyEOS_muhat_key, compose_direct },
  { "munu", NRPyEOS_munu_key, compose_direct },
};

static void check_hdf5(const herr_t status, const char *operation) {
  if(status < 0) {
    ghl_error("HDF5 operation failed: %s\n", operation);
  }
}

static void check_close(
      const char *quantity,
      const double expected,
      const double actual,
      const double rtol,
      const double atol) {
  if(!isfinite(expected) || !isfinite(actual)
     || (fabs(expected - actual) > atol
         && fabs(expected - actual) > rtol * fmax(fabs(expected), fabs(actual)))) {
    ghl_error("%s mismatch: expected %.17e, got %.17e\n", quantity, expected, actual);
  }
}

static void read_int_scalar(hid_t file, const char *name, int32_t *value) {
  hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if(dataset < 0) {
    ghl_error("Could not open dataset '%s'\n", name);
  }
  check_hdf5(
        H5Dread(dataset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, value), name);
  check_hdf5(H5Dclose(dataset), name);
}

static void read_double_dataset(
      hid_t file,
      const char *name,
      const int expected_rank,
      const hsize_t *expected_dims,
      double *values) {
  hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if(dataset < 0) {
    ghl_error("Could not open dataset '%s'\n", name);
  }
  hid_t space = H5Dget_space(dataset);
  if(space < 0) {
    ghl_error("Could not inspect dataset '%s'\n", name);
  }
  const int rank = H5Sget_simple_extent_ndims(space);
  if(rank != expected_rank) {
    ghl_error(
          "Dataset '%s' rank mismatch: expected %d, got %d\n", name, expected_rank,
          rank);
  }
  hsize_t dims[3] = { 0, 0, 0 };
  check_hdf5(H5Sget_simple_extent_dims(space, dims, NULL), name);
  for(int d = 0; d < expected_rank; d++) {
    if(dims[d] != expected_dims[d]) {
      ghl_error(
            "Dataset '%s' dimension %d mismatch: expected %llu, got %llu\n", name, d,
            (unsigned long long)expected_dims[d], (unsigned long long)dims[d]);
    }
  }
  check_hdf5(
        H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values),
        name);
  check_hdf5(H5Sclose(space), name);
  check_hdf5(H5Dclose(dataset), name);
}

static double convert_value(const compose_conversion conversion, const double value) {
  switch(conversion) {
    case compose_direct:
      return value;
    case compose_log_energy:
      return value * log(10.0) + log(CGS_TO_CODE_ENERGY);
    case compose_log_pressure:
      return value * log(10.0) + log(CGS_TO_CODE_PRESSURE);
    case compose_cs2:
      return value * CGS_TO_CODE_LENGTH * CGS_TO_CODE_LENGTH
             / (CGS_TO_CODE_TIME * CGS_TO_CODE_TIME);
    case compose_dedt:
      return value * CGS_TO_CODE_ENERGY;
    case compose_dpderho:
      return value * CGS_TO_CODE_PRESSURE / CGS_TO_CODE_ENERGY;
    case compose_dpdrhoe:
      return value * CGS_TO_CODE_PRESSURE / CGS_TO_CODE_DENSITY;
  }
  ghl_error("Unknown CompOSE conversion %d\n", conversion);
  return NAN;
}

static double midpoint_value(
      const ghl_eos_parameters *restrict eos,
      const NRPyEOS_keys key,
      const int ir,
      const int it,
      const int iy) {
  double result = 0.0;
  for(int dy = 0; dy < 2; dy++) {
    for(int dt = 0; dt < 2; dt++) {
      for(int dr = 0; dr < 2; dr++) {
        result += eos->table_all[NRPYEOS_IDX3D(eos, ir + dr, it + dt, iy + dy, key)];
      }
    }
  }
  return result / 8.0;
}

static void check_six_value_interpolator(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      const int ir,
      const int it,
      const int iy) {
  double muhat = NAN, mu_e = NAN, mu_p = NAN, mu_n = NAN, X_n = NAN, X_p = NAN;
  ghl_error_codes_t err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, rho, Y_e, T, &muhat, &mu_e, &mu_p, &mu_n, &X_n, &X_p);
  ghl_abort_if_error(err);

  const NRPyEOS_keys keys[6] = {
    NRPyEOS_muhat_key, NRPyEOS_mu_e_key, NRPyEOS_mu_p_key,
    NRPyEOS_mu_n_key,  NRPyEOS_X_n_key,  NRPyEOS_X_p_key,
  };
  const double actual[6] = { muhat, mu_e, mu_p, mu_n, X_n, X_p };
  const char *names[6] = { "muhat", "mu_e", "mu_p", "mu_n", "Xn", "Xp" };
  for(int i = 0; i < 6; i++) {
    check_close(
          names[i], midpoint_value(eos, keys[i], ir, it, iy), actual[i], 4.0e-13,
          2.0e-13);
  }

  err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, nextafter(eos->table_rho_min, 0.0), Y_e, T, &muhat, &mu_e, &mu_p, &mu_n,
        &X_n, &X_p);
  if(err != ghl_error_table_min_rho) {
    ghl_error("Six-value interpolator accepted density below table\n");
  }
  err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, nextafter(eos->table_rho_max, INFINITY), Y_e, T, &muhat, &mu_e, &mu_p,
        &mu_n, &X_n, &X_p);
  if(err != ghl_error_table_max_rho) {
    ghl_error("Six-value interpolator accepted density above table\n");
  }
  err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, rho, nextafter(eos->table_Y_e_min, -INFINITY), T, &muhat, &mu_e, &mu_p,
        &mu_n, &X_n, &X_p);
  if(err != ghl_error_table_min_ye) {
    ghl_error("Six-value interpolator accepted electron fraction below table\n");
  }
  err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, rho, nextafter(eos->table_Y_e_max, INFINITY), T, &muhat, &mu_e, &mu_p,
        &mu_n, &X_n, &X_p);
  if(err != ghl_error_table_max_ye) {
    ghl_error("Six-value interpolator accepted electron fraction above table\n");
  }
  err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, rho, Y_e, nextafter(eos->table_T_min, 0.0), &muhat, &mu_e, &mu_p, &mu_n,
        &X_n, &X_p);
  if(err != ghl_error_table_min_T) {
    ghl_error("Six-value interpolator accepted temperature below table\n");
  }
  err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, rho, Y_e, nextafter(eos->table_T_max, INFINITY), &muhat, &mu_e, &mu_p,
        &mu_n, &X_n, &X_p);
  if(err != ghl_error_table_max_T) {
    ghl_error("Six-value interpolator accepted temperature above table\n");
  }
}

static void check_inverses(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      const double stored_h) {
  double P = NAN, eps = NAN, S = NAN, cs2 = NAN;
  ghl_error_codes_t err
        = ghl_tabulated_compute_P_eps_S_cs2_from_T(eos, rho, Y_e, T, &P, &eps, &S, &cs2);
  ghl_abort_if_error(err);
  const double physical_h = 1.0 + eps + P / rho;
  if(!(P > 0.0) || !(eps + eos->energy_shift > 0.0) || !(physical_h > 0.0)
     || !isfinite(physical_h) || !(stored_h > 0.0) || !(cs2 > 0.0 && cs2 < 1.0)) {
    ghl_error("Forward EOS returned inadmissible state\n");
  }

  const double alternate_T
        = T == eos->table_T_min ? eos->table_T_max : eos->table_T_min;
  const double initial_T
        = exp(0.75 * log(T) + 0.25 * log(alternate_T));
  if(initial_T == T || initial_T < eos->table_T_min
        || initial_T > eos->table_T_max) {
    ghl_error("Inverse-search initial temperature is not distinct and valid\n");
  }
  const int inverse_keys[4] = {
    NRPyEOS_eps_key, NRPyEOS_press_key, NRPyEOS_entropy_key,
    NRPyEOS_enthalpy_key,
  };
  double initial_auxiliary[4] = {NAN, NAN, NAN, NAN};
  err = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(
        eos, 4, rho, Y_e, initial_T, inverse_keys, initial_auxiliary);
  ghl_abort_if_error(err);
  const double target_search_values[4] = {
    log(eps + eos->energy_shift), log(P), S, log(stored_h),
  };
  const double initial_search_values[4] = {
    log(initial_auxiliary[0] + eos->energy_shift),
    log(initial_auxiliary[1]),
    initial_auxiliary[2],
    log(initial_auxiliary[3]),
  };
  for(int i = 0; i < 4; i++) {
    if(fabs(initial_search_values[i] - target_search_values[i])
          <= eos->root_finding_precision * fabs(target_search_values[i])) {
      ghl_error("Inverse search %d would take the initial-match shortcut\n", i);
    }
  }

  double recovered_T = initial_T;
  double recovered_P = NAN;
  err = ghl_tabulated_compute_P_T_from_eps(
        eos, rho, Y_e, eps, &recovered_P, &recovered_T);
  ghl_abort_if_error(err);
  check_close("eps inverse temperature", T, recovered_T, 2.0e-9, 2.0e-12);
  check_close("eps inverse pressure", P, recovered_P, 2.0e-9, 2.0e-13);

  double recovered_eps = NAN, recovered_S = NAN;
  recovered_T = initial_T;
  err = ghl_tabulated_compute_eps_S_T_from_P(
        eos, rho, Y_e, P, &recovered_eps, &recovered_S, &recovered_T);
  ghl_abort_if_error(err);
  check_close("pressure inverse temperature", T, recovered_T, 2.0e-9, 2.0e-12);
  check_close("pressure inverse energy", eps, recovered_eps, 2.0e-9, 2.0e-13);
  check_close("pressure inverse entropy", S, recovered_S, 2.0e-9, 2.0e-13);

  recovered_T = initial_T;
  recovered_P = recovered_eps = NAN;
  err = ghl_tabulated_compute_P_eps_T_from_S(
        eos, rho, Y_e, S, &recovered_P, &recovered_eps, &recovered_T);
  ghl_abort_if_error(err);
  check_close("entropy inverse temperature", T, recovered_T, 2.0e-9, 2.0e-12);
  check_close("entropy inverse pressure", P, recovered_P, 2.0e-9, 2.0e-13);
  check_close("entropy inverse energy", eps, recovered_eps, 2.0e-9, 2.0e-13);

  recovered_T = initial_T;
  recovered_P = recovered_eps = recovered_S = NAN;
  err = ghl_tabulated_compute_P_eps_S_T_from_h(
        eos, rho, Y_e, stored_h, &recovered_P, &recovered_eps, &recovered_S,
        &recovered_T);
  ghl_abort_if_error(err);
  check_close("enthalpy inverse temperature", T, recovered_T, 2.0e-9, 2.0e-12);
  check_close("enthalpy inverse pressure", P, recovered_P, 2.0e-9, 2.0e-13);
  check_close("enthalpy inverse energy", eps, recovered_eps, 2.0e-9, 2.0e-13);
  check_close("enthalpy inverse entropy", S, recovered_S, 2.0e-9, 2.0e-13);
}

static void
check_synthetic_pressure_pava_inverse(const ghl_eos_parameters *restrict eos) {
  if(eos->N_rho != 5 || eos->N_T != 4 || eos->N_Ye != 3) {
    return;
  }

  const int ir = 2;
  const int target_it = 2;
  const int initial_it = 1;
  const int iy = 1;
  const double rho = exp(eos->table_logrho[ir]);
  const double Y_e = eos->table_Y_e[iy];
  const double target_T = exp(eos->table_logT[target_it]);
  const double initial_T = exp(eos->table_logT[initial_it]);

  double target_P = NAN, target_eps = NAN, target_S = NAN, target_cs2 = NAN;
  ghl_error_codes_t err = ghl_tabulated_compute_P_eps_S_cs2_from_T(
        eos, rho, Y_e, target_T, &target_P, &target_eps, &target_S, &target_cs2);
  ghl_abort_if_error(err);

  double recovered_T = initial_T;
  double recovered_eps = NAN, recovered_S = NAN;
  err = ghl_tabulated_compute_eps_S_T_from_P(
        eos, rho, Y_e, target_P, &recovered_eps, &recovered_S, &recovered_T);
  ghl_abort_if_error(err);
  check_close(
        "pressure PAVA inverse temperature", target_T, recovered_T, 2.0e-9, 2.0e-12);
  check_close(
        "pressure PAVA inverse energy", target_eps, recovered_eps, 2.0e-9, 2.0e-13);
  check_close("pressure PAVA inverse entropy", target_S, recovered_S, 2.0e-9, 2.0e-13);
}

static void check_downstream_consumers(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      const double P,
      const double eps,
      const double S,
      const double cs2) {
  const double velocity = 0.125;
  const double W = 1.0 / sqrt(1.0 - velocity * velocity);
  const double h = 1.0 + eps + P / rho;
  const double D = rho * W;
  const double momentum = rho * h * W * W * velocity;

  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&metric, &metric_aux);

  ghl_primitive_quantities original;
  ghl_initialize_primitives(
        rho, P, eps, velocity, 0.0, 0.0, 0.0, 0.0, 0.0, S, Y_e, T, &original);
  original.u0 = W;

  const ghl_con2prim_id_t no_backups[3] = {
    ghl_con2prim_id_None,
    ghl_con2prim_id_None,
    ghl_con2prim_id_None,
  };
  ghl_parameters params;
  ghl_initialize_params(
        ghl_con2prim_id_Palenzuela1D, no_backups, false, true, true, 1.0e100, 10.0, 0.0,
        &params);
  ghl_conservative_quantities conservative;
  ghl_initialize_conservatives(
        D, rho * h * W * W - P - D, momentum, 0.0, 0.0, S * W, D * Y_e, &conservative);
  ghl_primitive_quantities recovered = { 0 };
  recovered.u0 = W;
  ghl_con2prim_diagnostics diagnostics;
  ghl_initialize_diagnostics(&diagnostics);
  ghl_error_codes_t err = ghl_con2prim_tabulated_multi_method(
        &params, eos, &metric, &metric_aux, &conservative, &recovered, &diagnostics);
  ghl_abort_if_error(err);
  if(diagnostics.which_routine != ghl_con2prim_id_Palenzuela1D || diagnostics.backup[0]
     || diagnostics.backup[1] || diagnostics.backup[2] || diagnostics.nn_guess_used
     || diagnostics.speed_limited || diagnostics.tau_fix || diagnostics.Stilde_fix) {
    ghl_error("Regularized-table Con2Prim used a fallback or repair path\n");
  }
  check_close("Con2Prim density", rho, recovered.rho, 2.0e-8, 2.0e-13);
  check_close("Con2Prim pressure", P, recovered.press, 2.0e-8, 2.0e-13);
  check_close("Con2Prim energy", eps, recovered.eps, 2.0e-8, 2.0e-13);
  check_close("Con2Prim temperature", T, recovered.temperature, 2.0e-8, 2.0e-12);
  check_close("Con2Prim entropy", S, recovered.entropy, 2.0e-8, 2.0e-13);
  check_close("Con2Prim electron fraction", Y_e, recovered.Y_e, 2.0e-10, 2.0e-13);
  check_close("Con2Prim x velocity", velocity, recovered.vU[0], 2.0e-8, 2.0e-13);
  check_close("Con2Prim y velocity", 0.0, recovered.vU[1], 0.0, 2.0e-13);
  check_close("Con2Prim z velocity", 0.0, recovered.vU[2], 0.0, 2.0e-13);

  ghl_primitive_quantities left = original;
  ghl_primitive_quantities right = original;
  double cmin = NAN, cmax = NAN;
  ghl_calculate_characteristic_speed_dirn0(&right, &left, eos, &metric, &cmin, &cmax);
  const double sound_speed = sqrt(cs2);
  const double lambda_plus = (velocity + sound_speed) / (1.0 + velocity * sound_speed);
  const double lambda_minus = (velocity - sound_speed) / (1.0 - velocity * sound_speed);
  check_close(
        "positive characteristic speed", fmax(lambda_plus, 0.0), cmax, 2.0e-12, 2.0e-13);
  check_close(
        "negative characteristic speed", fmax(-lambda_minus, 0.0), cmin, 2.0e-12,
        2.0e-13);

  ghl_conservative_quantities flux;
  ghl_calculate_HLLE_fluxes_dirn0_tabulated(
        &right, &left, eos, &metric, cmin, cmax, &flux);
  check_close("mass flux", D * velocity, flux.rho, 2.0e-12, 2.0e-13);
  check_close("electron fraction flux", D * Y_e * velocity, flux.Y_e, 2.0e-12, 2.0e-13);
  check_close(
        "energy flux", (rho * h * W * W - D) * velocity, flux.tau, 2.0e-12, 2.0e-13);
  check_close("x momentum flux", momentum * velocity + P, flux.SD[0], 2.0e-12, 2.0e-13);
  check_close("y momentum flux", 0.0, flux.SD[1], 0.0, 2.0e-13);
  check_close("z momentum flux", 0.0, flux.SD[2], 0.0, 2.0e-13);
  ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy(
        &right, &left, eos, &metric, cmin, cmax, &flux);
  check_close("entropy flux", S * W * velocity, flux.entropy, 2.0e-12, 2.0e-13);

  ghl_metric_quantities derivative_x = { 0 };
  ghl_metric_quantities derivative_y = { 0 };
  ghl_metric_quantities derivative_z = { 0 };
  derivative_x.lapse = 0.01;
  derivative_y.lapse = -0.02;
  derivative_z.lapse = 0.03;
  ghl_extrinsic_curvature curvature;
  ghl_initialize_extrinsic_curvature(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &curvature);
  ghl_conservative_quantities source;
  ghl_calculate_source_terms(
        eos, &original, &metric, &derivative_x, &derivative_y, &derivative_z, &curvature,
        &source);
  const double energy_density = rho * h * W * W - P;
  check_close(
        "x momentum source", -0.01 * energy_density, source.SD[0], 2.0e-12, 2.0e-13);
  check_close(
        "y momentum source", 0.02 * energy_density, source.SD[1], 2.0e-12, 2.0e-13);
  check_close(
        "z momentum source", -0.03 * energy_density, source.SD[2], 2.0e-12, 2.0e-13);
  check_close("energy source", -0.01 * momentum, source.tau, 2.0e-12, 2.0e-13);

  const ghl_neutrino_optical_depths optical_depths = {
    .nue = { 0.01, 0.02 },
    .anue = { 0.03, 0.04 },
    .nux = { 0.05, 0.06 },
  };
  ghl_neutrino_opacities opacities;
  double R_source = NAN, Q_source = NAN;
  err = NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
        eos, rho, Y_e, T, &optical_depths, &opacities, &R_source, &Q_source);
  ghl_abort_if_error(err);
  const double leakage_values[8] = {
    opacities.nue[0], opacities.nue[1], opacities.anue[0], opacities.anue[1],
    opacities.nux[0], opacities.nux[1], R_source,          Q_source,
  };
  // Fixed regression goldens for the two qualified table dimensions.
  // The tolerance covers the production TEH fits.
  static const double analytic_expected[8] = {
    1.29265593927696647e-03,  2.06728788589579618e-03,  2.10755462947867636e-02,
    3.43634428155397380e-02,  1.27577464301964489e-03,  2.05212276868065074e-03,
    -1.28383988566382909e-17, -5.54094954548167941e-21,
  };
  static const double sro141_expected[8] = {
    2.58057523839381463e-05,  4.16887059307055122e-05,  6.22605413998306228e-05,
    9.86193518381971975e-05,  2.58186010234441992e-05,  4.15300141804711669e-05,
    -8.13615027117866091e-19, -6.05367209040713332e-22,
  };
  const double *expected = NULL;
  if(eos->N_rho == 5 && eos->N_T == 4 && eos->N_Ye == 3) {
    expected = analytic_expected;
  }
  else if(eos->N_rho == 391 && eos->N_T == 163 && eos->N_Ye == 66) {
    expected = sro141_expected;
  }
  else {
    ghl_error("No combined-leakage golden exists for table dimensions\n");
  }
  static const char *names[8] = {
    "nue opacity 0", "nue opacity 1", "anue opacity 0", "anue opacity 1",
    "nux opacity 0", "nux opacity 1", "R source",       "Q source",
  };
  for(int i = 0; i < 8; i++) {
    check_close(names[i], expected[i], leakage_values[i], 6.0e-4, DBL_MIN);
  }
}

int main(int argc, char **argv) {
  if(argc != 2) {
    ghl_error("Usage: %s REGULARIZED_TABLE.h5\n", argv[0]);
  }

  hid_t file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file < 0) {
    ghl_error("Could not open table '%s'\n", argv[1]);
  }
  int32_t n_rho = 0, n_T = 0, n_Ye = 0, have_rel_cs2 = 0;
  read_int_scalar(file, "pointsrho", &n_rho);
  read_int_scalar(file, "pointstemp", &n_T);
  read_int_scalar(file, "pointsye", &n_Ye);
  read_int_scalar(file, "have_rel_cs2", &have_rel_cs2);
  if(n_rho < 2 || n_T < 2 || n_Ye < 2 || have_rel_cs2 != 1) {
    ghl_error("Invalid regularized table dimensions or sound-speed convention\n");
  }

  double *logrho = malloc((size_t)n_rho * sizeof(*logrho));
  double *logT = malloc((size_t)n_T * sizeof(*logT));
  double *Y_e = malloc((size_t)n_Ye * sizeof(*Y_e));
  if(logrho == NULL || logT == NULL || Y_e == NULL) {
    ghl_error("Could not allocate axis buffers\n");
  }
  const hsize_t rho_dims[1] = { (hsize_t)n_rho };
  const hsize_t T_dims[1] = { (hsize_t)n_T };
  const hsize_t Ye_dims[1] = { (hsize_t)n_Ye };
  read_double_dataset(file, "logrho", 1, rho_dims, logrho);
  read_double_dataset(file, "logtemp", 1, T_dims, logT);
  read_double_dataset(file, "ye", 1, Ye_dims, Y_e);

  const double rho_atm = pow(10.0, logrho[n_rho / 2]) * CGS_TO_CODE_DENSITY;
  const double T_atm = pow(10.0, logT[n_T / 2]);
  const double Ye_atm = Y_e[n_Ye / 2];
  ghl_eos_parameters eos = { 0 };
  ghl_error_codes_t err = ghl_initialize_tabulated_eos_functions_and_params(
        argv[1], rho_atm, -1.0, -1.0, Ye_atm, -1.0, -1.0, T_atm, -1.0, -1.0, &eos);
  ghl_abort_if_error(err);
  if(eos.N_rho != n_rho || eos.N_T != n_T || eos.N_Ye != n_Ye || eos.clean_sound_speed
     || eos.c2p_nn != NULL) {
    ghl_error("Initializer changed regularized-table contract\n");
  }
  const ghl_eos_parameters *const eos_ptr = &eos;

  for(int i = 0; i < n_rho; i++) {
    check_close(
          "density axis", logrho[i] * log(10.0) + log(CGS_TO_CODE_DENSITY),
          eos.table_logrho[i], 4.0e-15, 4.0e-15);
  }
  for(int i = 0; i < n_T; i++) {
    check_close(
          "temperature axis", logT[i] * log(10.0), eos.table_logT[i], 4.0e-15, 4.0e-15);
  }
  for(int i = 0; i < n_Ye; i++) {
    check_close("electron-fraction axis", Y_e[i], eos.table_Y_e[i], 0.0, 0.0);
  }

  double raw_energy_shift = NAN;
  const hsize_t scalar_dims[1] = { 1 };
  read_double_dataset(file, "energy_shift", 1, scalar_dims, &raw_energy_shift);
  check_close(
        "energy shift", raw_energy_shift * CGS_TO_CODE_ENERGY, eos.energy_shift, 4.0e-15,
        4.0e-15);

  const size_t npoints = (size_t)n_rho * (size_t)n_T * (size_t)n_Ye;
  double *raw = malloc(npoints * sizeof(*raw));
  if(raw == NULL) {
    ghl_error("Could not allocate field buffer\n");
  }
  const hsize_t field_dims[3] = {
    (hsize_t)n_Ye,
    (hsize_t)n_T,
    (hsize_t)n_rho,
  };
  for(size_t field = 0; field < sizeof(fields) / sizeof(fields[0]); field++) {
    read_double_dataset(file, fields[field].name, 3, field_dims, raw);
    for(size_t point = 0; point < npoints; point++) {
      const double expected = convert_value(fields[field].conversion, raw[point]);
      const double actual
            = eos.table_all[fields[field].key + NRPyEOS_ntablekeys * point];
      check_close(fields[field].name, expected, actual, 8.0e-14, 2.0e-14);
    }
  }

  for(int iy = 0; iy < n_Ye; iy++) {
    for(int it = 0; it < n_T; it++) {
      for(int ir = 0; ir < n_rho; ir++) {
        const int point = NRPYEOS_IDX1D(eos_ptr, ir, it, iy);
        const double rho = exp(eos.table_logrho[ir]);
        const double P = exp(
              eos.table_all[NRPYEOS_IDX3D(eos_ptr, ir, it, iy, NRPyEOS_press_key)]);
        const double eps
              = exp(eos.table_all[NRPYEOS_IDX3D(eos_ptr, ir, it, iy, NRPyEOS_eps_key)])
                - eos.energy_shift;
        const double expected_logh = log(1.0 + eps + P / rho);
        check_close(
              "derived enthalpy", expected_logh,
              eos.table_all[NRPYEOS_IDX3D(eos_ptr, ir, it, iy, NRPyEOS_enthalpy_key)],
              8.0e-14, 2.0e-14);
        check_close(
              "enthalpy inversion storage", expected_logh, eos.table_logh[point],
              8.0e-14, 2.0e-14);
      }
    }
  }

  const int ir = (n_rho - 1) / 2;
  const int it = (n_T - 1) / 2;
  const int iy = (n_Ye - 1) / 2;
  const double rho = exp(0.5 * (eos.table_logrho[ir] + eos.table_logrho[ir + 1]));
  const double T = exp(0.5 * (eos.table_logT[it] + eos.table_logT[it + 1]));
  const double ye = 0.5 * (eos.table_Y_e[iy] + eos.table_Y_e[iy + 1]);

  double P = NAN, eps = NAN, S = NAN, cs2 = NAN;
  err = ghl_tabulated_compute_P_eps_S_cs2_from_T(&eos, rho, ye, T, &P, &eps, &S, &cs2);
  ghl_abort_if_error(err);
  check_close(
        "off-grid pressure", exp(midpoint_value(&eos, NRPyEOS_press_key, ir, it, iy)), P,
        4.0e-13, 2.0e-13);
  check_close(
        "off-grid energy",
        exp(midpoint_value(&eos, NRPyEOS_eps_key, ir, it, iy)) - eos.energy_shift, eps,
        4.0e-13, 2.0e-13);
  check_close(
        "off-grid entropy", midpoint_value(&eos, NRPyEOS_entropy_key, ir, it, iy), S,
        4.0e-13, 2.0e-13);
  check_close(
        "off-grid sound speed", midpoint_value(&eos, NRPyEOS_cs2_key, ir, it, iy), cs2,
        4.0e-13, 2.0e-13);
  check_six_value_interpolator(&eos, rho, ye, T, ir, it, iy);
  check_downstream_consumers(&eos, rho, ye, T, P, eps, S, cs2);

  double logh_min = 0.0, logh_max = 0.0;
  for(int dy = 0; dy < 2; dy++) {
    for(int dr = 0; dr < 2; dr++) {
      logh_min += eos.table_all[NRPYEOS_IDX3D(
            eos_ptr, ir + dr, 0, iy + dy, NRPyEOS_enthalpy_key)];
      logh_max += eos.table_all[NRPYEOS_IDX3D(
            eos_ptr, ir + dr, n_T - 1, iy + dy, NRPyEOS_enthalpy_key)];
    }
  }
  check_inverses(&eos, rho, ye, eos.table_T_min, exp(logh_min / 4.0));
  check_inverses(
        &eos, rho, ye, T, exp(midpoint_value(&eos, NRPyEOS_enthalpy_key, ir, it, iy)));
  check_inverses(&eos, rho, ye, eos.table_T_max, exp(logh_max / 4.0));
  check_synthetic_pressure_pava_inverse(&eos);

  free(raw);
  free(logrho);
  free(logT);
  free(Y_e);
  check_hdf5(H5Fclose(file), "close table");
  ghl_tabulated_free_memory(&eos);
  ghl_info("CompOSE regularized-table integration test passed\n");
  return 0;
}
