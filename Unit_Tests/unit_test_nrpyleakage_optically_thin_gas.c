#include "ghl_unit_tests.h"

//#define GENERATE_ASCII_DATA
#define Y_E 0
#define EPS 1

static int rhs_call_count;

static inline
void
ghl_pert_test_fail_computed_values(
    FILE *fp_unperturbed,
    FILE *fp_perturbed,
    const double t,
    const double Y_e,
    const double eps,
    const double T ) {

  double t_unperturbed, Y_e_unperturbed, eps_unperturbed, T_unperturbed;
  double t_perturbed, Y_e_perturbed, eps_perturbed, T_perturbed;

  // Read time, Y_e, eps, temperature from first file
  int err = 0;
  err += fread(&t_unperturbed  , sizeof(double), 1, fp_unperturbed);
  err += fread(&Y_e_unperturbed, sizeof(double), 1, fp_unperturbed);
  err += fread(&eps_unperturbed, sizeof(double), 1, fp_unperturbed);
  err += fread(&T_unperturbed  , sizeof(double), 1, fp_unperturbed);
  if( err != 4 )
    ghl_error("Failed to read unperturbed data from file\n");

  // Read time, Y_e, eps, temperature from second file
  err = 0;
  err += fread(&t_perturbed  , sizeof(double), 1, fp_perturbed);
  err += fread(&Y_e_perturbed, sizeof(double), 1, fp_perturbed);
  err += fread(&eps_perturbed, sizeof(double), 1, fp_perturbed);
  err += fread(&T_perturbed  , sizeof(double), 1, fp_perturbed);
  if( err != 4 )
    ghl_error("Failed to read perturbed data from file\n");

  // Perform validation
  if(ghl_pert_test_fail(t_unperturbed, t, t_perturbed)) {
    ghl_error(
          "Validation failed for t at t = %.17e: trusted %.17e, computed %.17e, "
          "perturbed %.17e\n",
          t, t_unperturbed, t, t_perturbed);
  }
  if(ghl_pert_test_fail(Y_e_unperturbed, Y_e, Y_e_perturbed)) {
    ghl_error(
          "Validation failed for Y_e at t = %.17e: trusted %.17e, computed %.17e, "
          "perturbed %.17e\n",
          t, Y_e_unperturbed, Y_e, Y_e_perturbed);
  }
  if(ghl_pert_test_fail(eps_unperturbed, eps, eps_perturbed)) {
    ghl_error(
          "Validation failed for eps at t = %.17e: trusted %.17e, computed %.17e, "
          "perturbed %.17e\n",
          t, eps_unperturbed, eps, eps_perturbed);
  }
  if(ghl_pert_test_fail(T_unperturbed, T, T_perturbed)) {
    ghl_error(
          "Validation failed for T at t = %.17e: trusted %.17e, computed %.17e, "
          "perturbed %.17e\n",
          t, T_unperturbed, T, T_perturbed);
  }
}

static inline ghl_error_codes_t
rhs(const ghl_eos_parameters *restrict eos,
    const double rho,
    const double Y_e,
    const double eps,
    const double T,
    double *restrict rhs_gfs) {
  ghl_neutrino_optical_depths tau = {{0,0},{0,0},{0,0}};
  ghl_neutrino_opacities kappa;
  double R_source, Q_source;
  rhs_call_count++;
  const ghl_error_codes_t error =
    NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(eos, rho, Y_e, T,
                                                                  &tau, &kappa, &R_source, &Q_source);
  if(error != ghl_success) {
    return error;
  }

  rhs_gfs[Y_E] = R_source/rho;
  rhs_gfs[EPS] = Q_source/rho;
  return ghl_success;
}

static inline ghl_error_codes_t rk4_step_ode(
      const ghl_eos_parameters *restrict eos,
      const double dt,
      const double rho,
      double *restrict gfs,
      double *restrict T) {

  // RK4 (no explicit time dependence on rhs):
  //
  // k1 = dt rhs(y(t))
  // k2 = dt rhs(y(t)+k1/2)
  // k3 = dt rhs(y(t)+k2/2)
  // k4 = dt rhs(y(t)+k3)
  //
  // y(t+dt) = y(t) + (1/6)( k1 + 2(k2 + k3) + k4 )
  double Y_e, eps, k1[2]={0,0}, k2[2]={0,0}, k3[2]={0,0}, k4[2]={0,0};

  // RK4 - substep 1
  *T = eos->T_max;
  Y_e = gfs[Y_E];
  eps = gfs[EPS];
  ghl_error_codes_t error = ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  if(error != ghl_success) {
    return error;
  }
  error = rhs(eos, rho, Y_e, eps, *T, k1);
  if(error != ghl_success) {
    return error;
  }

  // RK4 - substep 2;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + 0.5*dt*k1[Y_E];
  eps = gfs[EPS] + 0.5*dt*k1[EPS];
  error = ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  if(error != ghl_success) {
    return error;
  }
  error = rhs(eos, rho, Y_e, eps, *T, k2);
  if(error != ghl_success) {
    return error;
  }

  // RK4 - substep 3;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + 0.5*dt*k2[Y_E];
  eps = gfs[EPS] + 0.5*dt*k2[EPS];
  error = ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  if(error != ghl_success) {
    return error;
  }
  error = rhs(eos, rho, Y_e, eps, *T, k3);
  if(error != ghl_success) {
    return error;
  }

  // RK4 - substep 4;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + dt*k3[Y_E];
  eps = gfs[EPS] + dt*k3[EPS];
  error = ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  if(error != ghl_success) {
    return error;
  }
  error = rhs(eos, rho, Y_e, eps, *T, k4);
  if(error != ghl_success) {
    return error;
  }

  // RK4 - update step
  for(int i=0;i<2;i++)
    gfs[i] += (dt/6.0)*( k1[i] + 2.0*( k2[i] + k3[i] ) + k4[i] );
  return ghl_success;
}

static ghl_error_codes_t generate_one_fixture(
      const ghl_eos_parameters *restrict eos,
      const int perturb,
      const char *restrict filename,
      const double t_final,
      const double dt,
      const int n_steps) {
  double initial_rho = 1e-12;
  double initial_Y_e = 0.5;
  double initial_T = 1.0;
  if(perturb) {
    initial_rho *= (1 + randf(-1, 1) * 1e-14);
    initial_Y_e *= (1 + randf(-1, 1) * 1e-14);
    initial_T *= (1 + randf(-1, 1) * 1e-14);
  }

  double eps;
  ghl_error_codes_t error = ghl_tabulated_compute_eps_from_T(
        eos, initial_rho, initial_Y_e, initial_T, &eps);
  if(error != ghl_success) {
    return error;
  }

  FILE *fp = fopen_with_check(filename, "wb");
  double gfs[2] = { initial_Y_e, eps };
  fwrite(&n_steps, sizeof(int), 1, fp);
  if(!perturb) {
    fwrite(&dt, sizeof(double), 1, fp);
    fwrite(&t_final, sizeof(double), 1, fp);
    fwrite(&initial_rho, sizeof(double), 1, fp);
    fwrite(&initial_T, sizeof(double), 1, fp);
    fwrite(&gfs[Y_E], sizeof(double), 1, fp);
    fwrite(&gfs[EPS], sizeof(double), 1, fp);
  }
  double t = 0.0;
  for(int n = 0; n < n_steps; n++) {
    double T;
    error = rk4_step_ode(eos, dt, initial_rho, gfs, &T);
    if(error != ghl_success) {
      fclose(fp);
      return error;
    }
    t += dt;
    fwrite(&t, sizeof(double), 1, fp);
    fwrite(&gfs[Y_E], sizeof(double), 1, fp);
    fwrite(&gfs[EPS], sizeof(double), 1, fp);
    fwrite(&T, sizeof(double), 1, fp);
  }
  fclose(fp);
  return ghl_success;
}

void
generate_test_data(const ghl_eos_parameters *restrict eos) {

  const double t_final = 0.5*NRPyLeakage_units_cgs_to_geom_T;
  const double dt      = 0.001*NRPyLeakage_units_cgs_to_geom_T;
  const int n_steps    = (int)(t_final/dt+0.5);

  for(int perturb=0;perturb<=1;perturb++) {
    char filename[64];
    if( perturb )
      sprintf(filename, "nrpyleakage_optically_thin_gas_perturbed.bin");
    else
      sprintf(filename, "nrpyleakage_optically_thin_gas_unperturbed.bin");

    ghl_abort_if_error(
          generate_one_fixture(eos, perturb, filename, t_final, dt, n_steps));
    ghl_info("Finished %s evolution\n", perturb ? "perturbed" : "unperturbed");
  }
}

static int injected_T_failure_call;
static int injected_T_call_count;
static ghl_error_codes_t (*saved_compute_T_from_eps)(
      const ghl_eos_parameters *restrict,
      double,
      double,
      double,
      double *restrict);

static ghl_error_codes_t injected_compute_T_from_eps(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict T) {
  injected_T_call_count++;
  if(injected_T_call_count == injected_T_failure_call) {
    return ghl_error_table_max_T;
  }
  return saved_compute_T_from_eps(eos, rho, Y_e, eps, T);
}

static ghl_error_codes_t injected_compute_eps_from_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict eps) {
  (void)eos;
  (void)rho;
  (void)Y_e;
  (void)T;
  (void)eps;
  return ghl_error_table_max_T;
}

static void check_lookup_failures(const ghl_eos_parameters *restrict eos) {
  double initial_eps;
  ghl_abort_if_error(
        ghl_tabulated_compute_eps_from_T(eos, 1e-12, 0.5, 1.0, &initial_eps));
  saved_compute_T_from_eps = ghl_tabulated_compute_T_from_eps;
  ghl_tabulated_compute_T_from_eps = injected_compute_T_from_eps;
  for(int failure_call = 1; failure_call <= 4; failure_call++) {
    double gfs[2] = { 0.5, initial_eps };
    double T;
    injected_T_failure_call = failure_call;
    injected_T_call_count = 0;
    rhs_call_count = 0;
    const ghl_error_codes_t error = rk4_step_ode(eos, 1.0e-8, 1e-12, gfs, &T);
    if(error != ghl_error_table_max_T || rhs_call_count != failure_call - 1) {
      ghl_tabulated_compute_T_from_eps = saved_compute_T_from_eps;
      ghl_error(
            "RK lookup failure %d returned %d after %d RHS calls\n", failure_call, error,
            rhs_call_count);
    }
  }
  ghl_tabulated_compute_T_from_eps = saved_compute_T_from_eps;

  const char marker_path[] = "nrpyleakage_initial_lookup_failure_marker.bin";
  const unsigned char marker[] = { 0x47, 0x48, 0x4c };
  FILE *marker_file = fopen_with_check(marker_path, "wb");
  fwrite(marker, sizeof(marker), 1, marker_file);
  fclose(marker_file);
  ghl_error_codes_t (*saved_compute_eps_from_T)(
        const ghl_eos_parameters *restrict, double, double, double, double *restrict)
        = ghl_tabulated_compute_eps_from_T;
  ghl_tabulated_compute_eps_from_T = injected_compute_eps_from_T;
  const ghl_error_codes_t error = generate_one_fixture(eos, 0, marker_path, 0.0, 1.0, 0);
  ghl_tabulated_compute_eps_from_T = saved_compute_eps_from_T;
  unsigned char observed[sizeof(marker)] = { 0 };
  marker_file = fopen_with_check(marker_path, "rb");
  const size_t items = fread(observed, sizeof(observed), 1, marker_file);
  fclose(marker_file);
  remove(marker_path);
  if(error != ghl_error_table_max_T || items != 1
     || memcmp(marker, observed, sizeof(marker)) != 0) {
    ghl_error("Initial energy lookup failure modified its output file\n");
  }
}

void
run_unit_test(const ghl_eos_parameters *restrict eos) {
  check_lookup_failures(eos);
  int n1, n2;

  FILE *fp_unpert = fopen_with_check("nrpyleakage_optically_thin_gas_unperturbed.bin", "rb");
  FILE *fp_pert   = fopen_with_check("nrpyleakage_optically_thin_gas_perturbed.bin", "rb");

  int err = 0;
  err += fread(&n1, sizeof(int), 1, fp_unpert);
  err += fread(&n2, sizeof(int), 1, fp_pert  );
  if( err != 2 || n1 != 500 || n2 != 500 ) {
    fclose(fp_unpert); fclose(fp_pert);
    ghl_error("Invalid nrpyleakage_optically_thin_gas_{unperturbed,perturbed}.bin length "
              "(err: %d, n1: %d, n2: %d; expected 500)\n",
                 err, n1, n2);
  }

  const int n_steps = n1;
  double dt, t_final, initial_rho, initial_T, gfs[2];
  err  = 0;
  err += fread(&dt         , sizeof(double), 1, fp_unpert);
  err += fread(&t_final    , sizeof(double), 1, fp_unpert);
  err += fread(&initial_rho, sizeof(double), 1, fp_unpert);
  err += fread(&initial_T  , sizeof(double), 1, fp_unpert);
  err += fread(&gfs[Y_E]   , sizeof(double), 1, fp_unpert);
  err += fread(&gfs[EPS]   , sizeof(double), 1, fp_unpert);
  if( err != 6 ) {
    fclose(fp_unpert); fclose(fp_pert);
    ghl_error("Failed to read initial data from unperturbed data file\n");
  }
  double t = 0.0;
  for(int n=0;n<n_steps;n++) {
    double T;
    ghl_abort_if_error(rk4_step_ode(eos, dt, initial_rho, gfs, &T));
    t += dt;
    ghl_pert_test_fail_computed_values(fp_unpert, fp_pert, t, gfs[Y_E], gfs[EPS], T);
  }
}

#include "nrpyleakage_main.h"
