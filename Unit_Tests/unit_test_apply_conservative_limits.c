#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  FILE* infile = fopen_with_check("metric_Bfield_initial_data.bin","rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);
  if( key != 1 || arraylength < 1 )
    ghl_error("An error has occured with reading the grid size. "
                 "Please check that metric_initial_data.bin"
                 "is up-to-date with current test version.\n");

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const ghl_con2prim_id_t None = ghl_con2prim_id_None;
  const ghl_con2prim_id_t backup_routine[3] = {None, None, None};
  const bool evolve_entropy = false;
  const bool evolve_temperature = false;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const double W_max = 10.0;

  const int neos = 1;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, W_max, 0.0, &params);

  ghl_eos_parameters eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);


  // Allocate memory for the metric data
  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(lapse, sizeof(double), arraylength, infile);
  key += fread(betax, sizeof(double), arraylength, infile);
  key += fread(betay, sizeof(double), arraylength, infile);
  key += fread(betaz, sizeof(double), arraylength, infile);

  key += fread(gxx, sizeof(double), arraylength, infile);
  key += fread(gxy, sizeof(double), arraylength, infile);
  key += fread(gxz, sizeof(double), arraylength, infile);
  key += fread(gyy, sizeof(double), arraylength, infile);
  key += fread(gyz, sizeof(double), arraylength, infile);
  key += fread(gzz, sizeof(double), arraylength, infile);

  key += fread(Bx, sizeof(double), arraylength, infile);
  key += fread(By, sizeof(double), arraylength, infile);
  key += fread(Bz, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*13)
    ghl_error("An error has occured with reading in metric data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the initial conservative data
  double *rho_star = (double*) malloc(sizeof(double)*arraylength);
  double *tau = (double*) malloc(sizeof(double)*arraylength);
  double *S_x = (double*) malloc(sizeof(double)*arraylength);
  double *S_y = (double*) malloc(sizeof(double)*arraylength);
  double *S_z = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("apply_conservative_limits_input.bin","rb");
  key  = fread(rho_star, sizeof(double), arraylength, infile);
  key += fread(tau, sizeof(double), arraylength, infile);
  key += fread(S_x, sizeof(double), arraylength, infile);
  key += fread(S_y, sizeof(double), arraylength, infile);
  key += fread(S_z, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*5)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the trusted conservative data
  double *rho_star_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *tau_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_trusted = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("apply_conservative_limits_output.bin","rb");
  key  = fread(rho_star_trusted, sizeof(double), arraylength, infile);
  key += fread(tau_trusted, sizeof(double), arraylength, infile);
  key += fread(S_x_trusted, sizeof(double), arraylength, infile);
  key += fread(S_y_trusted, sizeof(double), arraylength, infile);
  key += fread(S_z_trusted, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*5)
    ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the perturbed conservative data
  double *rho_star_pert = (double*) malloc(sizeof(double)*arraylength);
  double *tau_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_pert = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("apply_conservative_limits_output_pert.bin","rb");
  key  = fread(rho_star_pert, sizeof(double), arraylength, infile);
  key += fread(tau_pert, sizeof(double), arraylength, infile);
  key += fread(S_x_pert, sizeof(double), arraylength, infile);
  key += fread(S_y_pert, sizeof(double), arraylength, infile);
  key += fread(S_z_pert, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*5)
    ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  const double poison = 0.0/0.0;

  for(int i=0;i<arraylength;i++) {
    // Define the various GRHayL structs for the unit tests
    ghl_con2prim_diagnostics diagnostics;
    diagnostics.tau_fix = true;
    diagnostics.Stilde_fix = true;
    diagnostics.speed_limited = true;
    diagnostics.backup[0] = true;
    diagnostics.backup[1] = true;
    diagnostics.backup[2] = true;
    diagnostics.nn_guess_used = true;
    diagnostics.which_routine = ghl_con2prim_id_Font1D;
    diagnostics.n_iter = -1;
    ghl_initialize_diagnostics(&diagnostics);
    if(diagnostics.tau_fix || diagnostics.Stilde_fix || diagnostics.speed_limited
          || diagnostics.backup[0] || diagnostics.backup[1] || diagnostics.backup[2]
          || diagnostics.nn_guess_used || diagnostics.n_iter != 0
          || diagnostics.which_routine != ghl_con2prim_id_None) {
      ghl_error("ghl_initialize_diagnostics did not clear every field\n");
    }
    ghl_metric_quantities metric_adm;
    ghl_primitive_quantities prims;
    ghl_conservative_quantities cons;

    // Read initial data accompanying trusted output
    ghl_initialize_metric(
          lapse[i],
          betax[i], betay[i], betaz[i],
          gxx[i], gxy[i], gxz[i],
          gyy[i], gyz[i], gzz[i],
          &metric_adm);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&metric_adm, &metric_aux);

    ghl_initialize_primitives(
          poison, poison, poison,
          poison, poison, poison,
          Bx[i], By[i], Bz[i],
          poison, poison, poison,
          &prims);

    ghl_initialize_conservatives(
          rho_star[i], tau[i],
          S_x[i], S_y[i], S_z[i],
          poison, poison, &cons);
    const ghl_conservative_quantities cons_before = cons;

    //This applies inequality fixes on the conservatives
    if(i == arraylength-1 || i == arraylength-2)
      params.psi6threshold = 0.0;
    ghl_apply_conservative_limits(&params, &eos, &metric_adm, &prims, &cons, &diagnostics);
    if(i == arraylength-1 || i == arraylength-2)
      params.psi6threshold = Psi6threshold;

    const bool tau_changed = cons.tau != cons_before.tau;
    const bool momentum_changed = cons.SD[0] != cons_before.SD[0]
                               || cons.SD[1] != cons_before.SD[1]
                               || cons.SD[2] != cons_before.SD[2];
    if(diagnostics.tau_fix != tau_changed
          || diagnostics.Stilde_fix != momentum_changed) {
      ghl_error("conservative-limit diagnostics disagree with actual changes at point %d\n", i);
    }

    ghl_conservative_quantities cons_trusted, cons_pert;
    ghl_initialize_conservatives(
          rho_star_trusted[i], tau_trusted[i],
          S_x_trusted[i], S_y_trusted[i], S_z_trusted[i],
          poison, poison, &cons_trusted);

    ghl_initialize_conservatives(
          rho_star_pert[i], tau_pert[i],
          S_x_pert[i], S_y_pert[i], S_z_pert[i],
          poison, poison, &cons_pert);


    ghl_pert_test_fail_conservatives(params.evolve_entropy, &cons_trusted, &cons, &cons_pert);
  }

  ghl_metric_quantities sticky_metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &sticky_metric);
  ghl_primitive_quantities sticky_prims;
  ghl_initialize_primitives(
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &sticky_prims);
  ghl_conservative_quantities sticky_cons;
  ghl_initialize_conservatives(
        1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &sticky_cons);
  ghl_con2prim_diagnostics sticky_diagnostics;
  ghl_initialize_diagnostics(&sticky_diagnostics);
  sticky_diagnostics.tau_fix = true;
  sticky_diagnostics.Stilde_fix = true;
  ghl_apply_conservative_limits(
        &params, &eos, &sticky_metric, &sticky_prims, &sticky_cons,
        &sticky_diagnostics);
  if(!sticky_diagnostics.tau_fix || !sticky_diagnostics.Stilde_fix) {
    ghl_error("conservative-limit diagnostics did not preserve incoming true values\n");
  }

  const double tau_inputs[3] = {
    0.5*eos.tau_atm, eos.tau_atm, 10.0*eos.tau_atm
  };
  for(int i = 0; i < 3; ++i) {
    ghl_conservative_quantities floor_cons;
    ghl_initialize_conservatives(
          1.0, tau_inputs[i], 0.0, 0.0, 0.0, 0.0, 0.0, &floor_cons);
    ghl_con2prim_diagnostics floor_diagnostics;
    ghl_initialize_diagnostics(&floor_diagnostics);
    ghl_apply_conservative_limits(
          &params, &eos, &sticky_metric, &sticky_prims, &floor_cons,
          &floor_diagnostics);
    const bool expect_floor = i == 0;
    if(floor_diagnostics.tau_fix != expect_floor
          || floor_diagnostics.Stilde_fix
          || floor_cons.tau != (expect_floor ? eos.tau_atm : tau_inputs[i])) {
      ghl_error("tau atmosphere floor contract failed for synthetic case %d\n", i);
    }
  }

  ghl_conservative_quantities momentum_cons;
  ghl_initialize_conservatives(
        1.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, &momentum_cons);
  ghl_con2prim_diagnostics momentum_diagnostics;
  ghl_initialize_diagnostics(&momentum_diagnostics);
  ghl_apply_conservative_limits(
        &params, &eos, &sticky_metric, &sticky_prims, &momentum_cons,
        &momentum_diagnostics);
  if(momentum_diagnostics.tau_fix || !momentum_diagnostics.Stilde_fix
        || !(momentum_cons.SD[0] < 10.0)) {
    ghl_error("isolated low-field momentum correction was not diagnosed\n");
  }

  ghl_primitive_quantities high_psi_prims = sticky_prims;
  high_psi_prims.BU[0] = 1.0;
  ghl_conservative_quantities high_psi_cons;
  ghl_initialize_conservatives(
        1.0, 10.0, 100.0, 0.0, 0.0, 0.0, 0.0, &high_psi_cons);
  ghl_con2prim_diagnostics high_psi_diagnostics;
  ghl_initialize_diagnostics(&high_psi_diagnostics);
  const double saved_psi6threshold = params.psi6threshold;
  params.psi6threshold = 0.0;
  ghl_apply_conservative_limits(
        &params, &eos, &sticky_metric, &high_psi_prims, &high_psi_cons,
        &high_psi_diagnostics);
  params.psi6threshold = saved_psi6threshold;
  if(high_psi_diagnostics.tau_fix || !high_psi_diagnostics.Stilde_fix
        || !(high_psi_cons.SD[0] < 100.0)) {
    ghl_error("isolated high-psi6 momentum correction was not diagnosed\n");
  }

  ghl_primitive_quantities magnetic_prims = sticky_prims;
  magnetic_prims.BU[0] = 2.0;
  ghl_conservative_quantities magnetic_cons;
  ghl_initialize_conservatives(
        1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &magnetic_cons);
  ghl_con2prim_diagnostics magnetic_diagnostics;
  ghl_initialize_diagnostics(&magnetic_diagnostics);
  ghl_apply_conservative_limits(
        &params, &eos, &sticky_metric, &magnetic_prims, &magnetic_cons,
        &magnetic_diagnostics);
  const double magnetic_tau = eos.tau_atm + 2.0;
  if(!magnetic_diagnostics.tau_fix || magnetic_diagnostics.Stilde_fix
        || magnetic_cons.tau != magnetic_tau) {
    ghl_error("isolated magnetic-energy tau correction was not diagnosed\n");
  }

  ghl_primitive_quantities psi6_prims = sticky_prims;
  psi6_prims.BU[0] = 1.0;
  ghl_eos_parameters psi6_eos = eos;
  psi6_eos.tau_atm = 1.0;
  ghl_conservative_quantities psi6_cons;
  const double psi6_tau = 0.5 + 1.0005*psi6_eos.tau_atm;
  ghl_initialize_conservatives(
        1.0, psi6_tau, 0.0, 0.0, 0.0, 0.0, 0.0, &psi6_cons);
  ghl_con2prim_diagnostics psi6_diagnostics;
  ghl_initialize_diagnostics(&psi6_diagnostics);
  const double saved_psi6threshold_for_tau = params.psi6threshold;
  params.psi6threshold = 0.0;
  ghl_apply_conservative_limits(
        &params, &psi6_eos, &sticky_metric, &psi6_prims, &psi6_cons,
        &psi6_diagnostics);
  params.psi6threshold = saved_psi6threshold_for_tau;
  const double corrected_psi6_tau = 0.5 + 1.001*psi6_eos.tau_atm;
  if(!psi6_diagnostics.tau_fix || psi6_diagnostics.Stilde_fix
        || psi6_cons.tau != corrected_psi6_tau) {
    ghl_error("isolated high-psi6 tau correction was not diagnosed\n");
  }
  ghl_info("ghl_apply_conservative_limits function test has passed!\n");
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(Bx); free(By); free(Bz);
  free(rho_star); free(tau);
  free(S_x); free(S_y); free(S_z);
  free(rho_star_trusted); free(tau_trusted);
  free(S_x_trusted); free(S_y_trusted); free(S_z_trusted);
  free(rho_star_pert); free(tau_pert);
  free(S_x_pert); free(S_y_pert); free(S_z_pert);
}
