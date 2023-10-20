#include "unit_tests.h"
#include <math.h>

int main(int argc, char **argv) {

  if(argc != 2) {
    ghl_info("Usage: %s <EOS table path>\n", argv[0]);
    exit(1);
  }

  const char *tablepath = argv[1];

  srand(0);
  const bool random_data = true;

  ghl_parameters params;
  params.max_lorenz_factor = 10.0;
  params.inv_sq_max_lorenz_factor = 1.0 / SQR(params.max_lorenz_factor);

  // ghl_eos_t eos_type = ghl_eos_hybrid;
  ghl_eos_t eos_type = ghl_eos_tabulated;

  ghl_eos_parameters eos;
  const double rho_atm = 1e-12;

  if(eos_type == ghl_eos_hybrid) {
    const double Gamma = 2.0;
    const double K = 1.0;
    ghl_initialize_hybrid_eos_functions_and_params(
        rho_atm, 1e-12, 1e-1, 1, &rho_atm, &Gamma, K, Gamma, &eos);
  }
  else if(eos_type == ghl_eos_tabulated) {
    const double Y_e_atm = 0.5;
    const double T_atm = 1e-2;
    ghl_initialize_tabulated_eos_functions_and_params(
        tablepath, rho_atm, -1, -1, Y_e_atm, -1, -1, T_atm, -1, -1, &eos);
    eos.root_finding_precision = 1e-10;
  }

  double alp = 1, betaU[3] = { 0, 0, 0 }, gDD[6] = { 1, 0, 0, 1, 0, 1 };

  if(random_data) {
    alp = randf(0.1, 1);
    betaU[0] = randf(-0.1, 0.1);
    betaU[1] = randf(-0.1, 0.1);
    betaU[2] = randf(-0.1, 0.1);
    gDD[0] += randf(0, 0.3);
    gDD[1] += randf(-0.1, 0.1);
    gDD[2] += randf(-0.1, 0.1);
    gDD[3] += randf(0, 0.3);
    gDD[4] += randf(-0.1, 0.1);
    gDD[5] += randf(0, 0.3);
  }

  ghl_metric_quantities metric;
  ghl_initialize_metric(
      alp, betaU[0], betaU[1], betaU[2], gDD[0], gDD[1], gDD[2], gDD[3], gDD[4], gDD[5], &metric);

  ghl_primitive_quantities prims_r, prims_l;
  if(eos_type == ghl_eos_tabulated) {
    prims_r.rho = exp(randf(log(eos.table_rho_min), log(eos.table_rho_max)));
    prims_r.temperature = exp(randf(log(eos.table_T_min), log(eos.table_T_max)));
    prims_r.Y_e = randf(eos.table_Y_e_min, eos.table_Y_e_max);
    prims_r.vU[0] = 0; // randf(-1,1);
    prims_r.vU[1] = 0; // randf(-1,1);
    prims_r.vU[2] = 0; // randf(-1,1);
    prims_l.rho = exp(randf(log(eos.table_rho_min), log(eos.table_rho_max)));
    prims_l.temperature = exp(randf(log(eos.table_T_min), log(eos.table_T_max)));
    prims_l.Y_e = randf(eos.table_Y_e_min, eos.table_Y_e_max);
    prims_l.vU[0] = 0; // randf(-1,1);
    prims_l.vU[1] = 0; // randf(-1,1);
    prims_l.vU[2] = 0; // randf(-1,1);
    ghl_tabulated_compute_P_eps_S_from_T(
        &eos, prims_r.rho, prims_r.Y_e, prims_r.temperature, &prims_r.press, &prims_r.eps,
        &prims_r.entropy);
    ghl_tabulated_compute_P_eps_S_from_T(
        &eos, prims_l.rho, prims_l.Y_e, prims_l.temperature, &prims_l.press, &prims_l.eps,
        &prims_l.entropy);
  }
  else if(eos_type == ghl_eos_hybrid) {
    prims_r.rho = 1.2345e-4;
    prims_l.rho = 5.4321e-3;
    prims_r.press = prims_r.rho * prims_r.rho;
    prims_l.press = prims_l.rho * prims_l.rho;
    prims_r.vU[0] = prims_l.vU[0] = 0;
    prims_r.vU[1] = prims_l.vU[1] = 0;
    prims_r.vU[2] = prims_l.vU[2] = 0;
  }

  prims_r.BU[0] = 0; // prims_r.press*pow(10, randf(-5, -3));
  prims_r.BU[1] = 0; // prims_r.press*pow(10, randf(-5, -3));
  prims_r.BU[2] = 0; // prims_r.press*pow(10, randf(-5, -3));
  prims_l.BU[0] = 0; // prims_l.press*pow(10, randf(-5, -3));
  prims_l.BU[1] = 0; // prims_l.press*pow(10, randf(-5, -3));
  prims_l.BU[2] = 0; // prims_l.press*pow(10, randf(-5, -3));
  ghl_limit_v_and_compute_u0(&params, &metric, &prims_r);
  ghl_limit_v_and_compute_u0(&params, &metric, &prims_l);

  printf(
      "Metric: %e, %e %e %e, %e %e %e %e %e %e\n", metric.lapse, metric.betaU[0], metric.betaU[1],
      metric.betaU[2], metric.gammaDD[0][0], metric.gammaDD[0][1], metric.gammaDD[0][2],
      metric.gammaDD[1][1], metric.gammaDD[1][2], metric.gammaDD[2][2]);

  printf(
      "R: %e %e %e -> %e %e %e, %e %e %e, %e %e %e\n", prims_r.rho, prims_r.Y_e,
      prims_r.temperature, prims_r.press, prims_r.eps, prims_r.entropy, prims_r.vU[0],
      prims_r.vU[1], prims_r.vU[2], prims_r.BU[0], prims_r.BU[1], prims_r.BU[2]);

  printf(
      "L: %e %e %e -> %e %e %e, %e %e %e, %e %e %e\n", prims_l.rho, prims_l.Y_e,
      prims_l.temperature, prims_l.press, prims_l.eps, prims_l.entropy, prims_l.vU[0],
      prims_l.vU[1], prims_l.vU[2], prims_l.BU[0], prims_l.BU[1], prims_l.BU[2]);

  double cmin, cmax;
  ghl_conservative_quantities flux;

  ghl_calculate_characteristic_speed_dirn0(&prims_r, &prims_l, &eos, &metric, &cmin, &cmax);
  if(eos_type == ghl_eos_tabulated) {
    ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy(
        &prims_r, &prims_l, &eos, &metric, cmin, cmax, &flux);
    fprintf(
        stderr, "0 %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
        cmin, cmax, flux.rho, flux.tau, flux.SD[0], flux.SD[1], flux.SD[2], flux.Y_e,
        flux.entropy);
  }
  else {
    ghl_calculate_HLLE_fluxes_dirn0_hybrid(&prims_r, &prims_l, &eos, &metric, cmin, cmax, &flux);
    fprintf(stderr, "0 %22.15e %22.15e\n", cmin, cmax);
  }

  ghl_calculate_characteristic_speed_dirn1(&prims_r, &prims_l, &eos, &metric, &cmin, &cmax);
  if(eos_type == ghl_eos_tabulated) {
    ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy(
        &prims_r, &prims_l, &eos, &metric, cmin, cmax, &flux);
    fprintf(
        stderr, "1 %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
        cmin, cmax, flux.rho, flux.tau, flux.SD[0], flux.SD[1], flux.SD[2], flux.Y_e,
        flux.entropy);
  }
  else {
    ghl_calculate_HLLE_fluxes_dirn1_hybrid(&prims_r, &prims_l, &eos, &metric, cmin, cmax, &flux);
    fprintf(stderr, "1 %22.15e %22.15e\n", cmin, cmax);
  }

  ghl_calculate_characteristic_speed_dirn2(&prims_r, &prims_l, &eos, &metric, &cmin, &cmax);
  if(eos_type == ghl_eos_tabulated) {
    ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy(
        &prims_r, &prims_l, &eos, &metric, cmin, cmax, &flux);
    fprintf(
        stderr, "2 %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
        cmin, cmax, flux.rho, flux.tau, flux.SD[0], flux.SD[1], flux.SD[2], flux.Y_e,
        flux.entropy);
  }
  else {
    ghl_calculate_HLLE_fluxes_dirn2_hybrid(&prims_r, &prims_l, &eos, &metric, cmin, cmax, &flux);
    fprintf(stderr, "2 %22.15e %22.15e\n", cmin, cmax);
  }

  if(eos_type == ghl_eos_tabulated) {
    ghl_tabulated_free_memory(&eos);
  }

  return 0;
}
