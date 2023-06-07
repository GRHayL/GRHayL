#include "neutrinos.h"
#include "unit_tests.h"

//#define GENERATE_ASCII_DATA
#define Y_E 0
#define EPS 1

static inline
void
validate_computed_values(
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
  validate(t_unperturbed  , t  , t_perturbed  );
  validate(Y_e_unperturbed, Y_e, Y_e_perturbed);
  validate(eps_unperturbed, eps, eps_perturbed);
  validate(T_unperturbed  , T  , T_perturbed  );
}

static inline
void
rhs(const eos_parameters *restrict eos,
    const double rho,
    const double Y_e,
    const double eps,
    const double T,
    double *restrict rhs_gfs ) {
  neutrino_optical_depths tau = {{0,0},{0,0},{0,0}};
  neutrino_opacities kappa;
  double R_source, Q_source;
  NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(eos, rho, Y_e, T,
                                                                &tau, &kappa, &R_source, &Q_source);

  rhs_gfs[Y_E] = R_source/rho;
  rhs_gfs[EPS] = Q_source/rho;
}

static inline
void
rk4_step_ode(
    const eos_parameters *restrict eos,
    const double dt,
    const double rho,
    double *restrict gfs,
    double *restrict T ) {

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
  ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  rhs(eos, rho, Y_e, eps, *T, k1);

  // RK4 - substep 2;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + 0.5*dt*k1[Y_E];
  eps = gfs[EPS] + 0.5*dt*k1[EPS];
  ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  rhs(eos, rho, Y_e, eps, *T, k2);

  // RK4 - substep 3;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + 0.5*dt*k2[Y_E];
  eps = gfs[EPS] + 0.5*dt*k2[EPS];
  ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  rhs(eos, rho, Y_e, eps, *T, k3);

  // RK4 - substep 4;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + dt*k3[Y_E];
  eps = gfs[EPS] + dt*k3[EPS];
  ghl_tabulated_compute_T_from_eps(eos, rho, Y_e, eps, T);
  rhs(eos, rho, Y_e, eps, *T, k4);

  // RK4 - update step
  for(int i=0;i<2;i++)
    gfs[i] += (dt/6.0)*( k1[i] + 2.0*( k2[i] + k3[i] ) + k4[i] );
}

void
generate_test_data(const eos_parameters *restrict eos) {

  const double t_final = 0.5*NRPyLeakage_units_cgs_to_geom_T;
  const double dt      = 0.001*NRPyLeakage_units_cgs_to_geom_T;
  const int n_steps    = (int)(t_final/dt+0.5);

  for(int perturb=0;perturb<=1;perturb++) {
    char filename[64];
    if( perturb )
      sprintf(filename, "nrpyleakage_optically_thin_gas_perturbed.bin");
    else
      sprintf(filename, "nrpyleakage_optically_thin_gas_unperturbed.bin");

    FILE *fp = fopen_with_check(filename, "wb");

    double initial_rho = 1e-12;
    double initial_Y_e = 0.5;
    double initial_T   = 1.0;
    if( perturb ) {
      initial_rho *= (1+randf(-1,1)*1e-14);
      initial_Y_e *= (1+randf(-1,1)*1e-14);
      initial_T   *= (1+randf(-1,1)*1e-14);
    }

    double eps;
    ghl_tabulated_compute_eps_from_T(eos, initial_rho, initial_Y_e, initial_T, &eps);

    double gfs[2] = {initial_Y_e, eps};
    fwrite(&n_steps  , sizeof(int)   , 1, fp);
    if( !perturb ) {
      fwrite(&dt         , sizeof(double), 1, fp);
      fwrite(&t_final    , sizeof(double), 1, fp);
      fwrite(&initial_rho, sizeof(double), 1, fp);
      fwrite(&initial_T  , sizeof(double), 1, fp);
      fwrite(&gfs[Y_E]   , sizeof(double), 1, fp);
      fwrite(&gfs[EPS]   , sizeof(double), 1, fp);
    }
    double t = 0.0;
    for(int n=0;n<n_steps;n++) {
      double T;
      rk4_step_ode(eos, dt, initial_rho, gfs, &T);
      t += dt;
      fwrite(&t        , sizeof(double), 1, fp);
      fwrite(&gfs[Y_E] , sizeof(double), 1, fp);
      fwrite(&gfs[EPS] , sizeof(double), 1, fp);
      fwrite(&T        , sizeof(double), 1, fp);
    }
    fclose(fp);
    ghl_info("Finished %s evolution\n", perturb ? "perturbed" : "unperturbed");
  }
}

void
run_unit_test(const eos_parameters *restrict eos) {
  int n1, n2;

  FILE *fp_unpert = fopen_with_check("nrpyleakage_optically_thin_gas_unperturbed.bin", "rb");
  FILE *fp_pert   = fopen_with_check("nrpyleakage_optically_thin_gas_perturbed.bin", "rb");

  int err = 0;
  err += fread(&n1, sizeof(int), 1, fp_unpert);
  err += fread(&n2, sizeof(int), 1, fp_pert  );
  if( err != 2 || n1 != n2 ) {
    fclose(fp_unpert); fclose(fp_pert);
    ghl_error("Problem reading number of steps from file (err: %d, n1: %d, n2: %d)\n",
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
    rk4_step_ode(eos, dt, initial_rho, gfs, &T);
    t += dt;
    validate_computed_values(fp_unpert, fp_pert, t, gfs[Y_E], gfs[EPS], T);
  }
}

#include "nrpyleakage_main.h"
