#include "GRHayL.h"
#include "Neutrinos.h"
#include "NRPyLeakage.h"
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
    grhayl_error("Failed to read unperturbed data from file\n");

  // Read time, Y_e, eps, temperature from second file
  err = 0;
  err += fread(&t_perturbed  , sizeof(double), 1, fp_perturbed);
  err += fread(&Y_e_perturbed, sizeof(double), 1, fp_perturbed);
  err += fread(&eps_perturbed, sizeof(double), 1, fp_perturbed);
  err += fread(&T_perturbed  , sizeof(double), 1, fp_perturbed);
  if( err != 4 )
    grhayl_error("Failed to read perturbed data from file\n");

  // Perform validation
  validate(t_unperturbed  , t  , t_perturbed  );
  validate(Y_e_unperturbed, Y_e, Y_e_perturbed);
  validate(eps_unperturbed, eps, eps_perturbed);
  validate(T_unperturbed  , T  , T_perturbed  );
}

static inline
void
rhs(const eos_parameters *restrict eos,
    const double rho_b,
    const double Y_e,
    const double eps,
    const double T,
    double *restrict rhs_gfs ) {
  neutrino_optical_depths tau = {{0,0},{0,0},{0,0}};
  neutrino_opacities kappa;
  double R_source, Q_source;
  NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(eos, rho_b, Y_e, T,
                                                                &tau, &kappa, &R_source, &Q_source);

  rhs_gfs[Y_E] = R_source/rho_b;
  rhs_gfs[EPS] = Q_source/rho_b;
}

static inline
void
rk4_step_ode(
    const eos_parameters *restrict eos,
    const double dt,
    const double rho_b,
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
  eos->tabulated_compute_T_from_eps(eos, rho_b, Y_e, eps, T);
  rhs(eos, rho_b, Y_e, eps, *T, k1);

  // RK4 - substep 2;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + 0.5*dt*k1[Y_E];
  eps = gfs[EPS] + 0.5*dt*k1[EPS];
  eos->tabulated_compute_T_from_eps(eos, rho_b, Y_e, eps, T);
  rhs(eos, rho_b, Y_e, eps, *T, k2);

  // RK4 - substep 3;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + 0.5*dt*k2[Y_E];
  eps = gfs[EPS] + 0.5*dt*k2[EPS];
  eos->tabulated_compute_T_from_eps(eos, rho_b, Y_e, eps, T);
  rhs(eos, rho_b, Y_e, eps, *T, k3);

  // RK4 - substep 4;
  *T = eos->T_max;
  Y_e = gfs[Y_E] + dt*k3[Y_E];
  eps = gfs[EPS] + dt*k3[EPS];
  eos->tabulated_compute_T_from_eps(eos, rho_b, Y_e, eps, T);
  rhs(eos, rho_b, Y_e, eps, *T, k4);

  // RK4 - update step
  for(int i=0;i<2;i++)
    gfs[i] += (dt/6.0)*( k1[i] + 2.0*( k2[i] + k3[i] ) + k4[i] );
}

int main( int argc, char **argv ) {

  if( argc != 3 ) {
    grhayl_info("Correct usage is: %s <eos table> <test key>\n", *argv);
    grhayl_info("Available test keys:\n");
    grhayl_info("  0 : generate unperturbed\n");
    grhayl_info("  1 : generate perturbed\n");
    grhayl_info("  2 : perform validation test\n");
    exit(1);
  }

  const char *tablepath = argv[1];
  int test_key = atoi(argv[2]);
  double initial_rho_b = 1e-12;
  double initial_Y_e   = 0.5;
  double initial_T     = 1.0;
  if( test_key == 1 ) {
    initial_rho_b *= (1+randf(-1,1)*1e-14);
    initial_Y_e   *= (1+randf(-1,1)*1e-14);
    initial_T     *= (1+randf(-1,1)*1e-14);
  }
  else if( test_key != 0 && test_key != 2 )
    grhayl_error("Unsupported test key %d\n", test_key);

  const double W_max     = 10.0;
  const double rho_b_atm = 1e-12;
  const double rho_b_min = -1;
  const double rho_b_max = -1;
  const double Y_e_atm   = 0.5;
  const double Y_e_min   = -1;
  const double Y_e_max   = -1;
  const double T_atm     = 1e-2;
  const double T_min     = -1;
  const double T_max     = -1;

  eos_parameters eos;
  initialize_tabulated_eos_functions_and_params(tablepath, W_max,
                                                rho_b_atm, rho_b_min, rho_b_max,
                                                Y_e_atm, Y_e_min, Y_e_max,
                                                T_atm, T_min, T_max, &eos);
  eos.root_finding_precision=1e-10;

  const double t_final = 0.5*NRPyLeakage_units_cgs_to_geom_T;
  const double dt      = 0.001*NRPyLeakage_units_cgs_to_geom_T;
  const int n_steps    = (int)(t_final/dt+0.5);
  double t = 0.0;

  double eps;
  eos.tabulated_compute_eps_from_T(&eos, initial_rho_b, initial_Y_e, initial_T, &eps);

  double gfs[2] = {initial_Y_e, eps};

#ifdef GENERATE_ASCII_DATA
  FILE *fp = fopen("opticallythingas_semi_analytic_grhayl.txt","w");
  fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T, gfs[Y_E], gfs[EPS], initial_T);
#else
  FILE *fp1=NULL, *fp2=NULL;
  if( test_key != 2 ) {
    if( test_key )
      fp1 = fopen("nrpyleakage_perturbed.bin","wb");
    else
      fp1 = fopen("nrpyleakage_unperturbed.bin","wb");

    fwrite(&t        , sizeof(double), 1, fp1);
    fwrite(&gfs[Y_E] , sizeof(double), 1, fp1);
    fwrite(&gfs[EPS] , sizeof(double), 1, fp1);
    fwrite(&initial_T, sizeof(double), 1, fp1);
  }
  else {
    fp1 = fopen("nrpyleakage_unperturbed.bin", "rb");
    fp2 = fopen("nrpyleakage_perturbed.bin"  , "rb");
  }
#endif
  for(int n=0;n<n_steps;n++) {
    double T;
    rk4_step_ode(&eos, dt, initial_rho_b, gfs, &T);
    t += dt;
#ifdef GENERATE_ASCII_DATA
    fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T, gfs[Y_E], gfs[EPS], T);
#else
    if( test_key == 2 )
      validate_computed_values(fp1, fp2, t, gfs[Y_E], gfs[EPS], T);
    else {
      fwrite(&t        , sizeof(double), 1, fp1);
      fwrite(&gfs[Y_E] , sizeof(double), 1, fp1);
      fwrite(&gfs[EPS] , sizeof(double), 1, fp1);
      fwrite(&T        , sizeof(double), 1, fp1);
    }
#endif
  }
  fclose(fp1);
  if( test_key == 2 ) fclose(fp2);
  eos.tabulated_free_memory(&eos);
}
