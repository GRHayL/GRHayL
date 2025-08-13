#include "ghl_unit_tests.h"
#include "ghl_radiation.h"
#include "ghl_m1.h"


void test_source_update();

//copied from unit_test_nrpyleakage_optically_thin_gas.c
int main (int argc, char **argv) {
  
  //initialize the closure of choice
  ghl_error_codes_t initialize_check = ghl_initialize_m1_closure(Minerbo);
  ghl_read_error_codes(initialize_check);

  printf("Run Source update test\n");
  test_source_update();

  return 0;
}

void test_source_update() {
  
  // Set metric quantities to Minkowski
  ghl_metric_quantities metric;
  ghl_initialize_metric(1, 0, 0, 0, //lapse, shift
                           1, 0, 0, // gxx gxy gxz
                              1, 0, //     gyy gyz
                                 1, //         gzz
                           &metric);

  ghl_ADM_aux_quantities adm_aux;
  ghl_compute_ADM_auxiliaries(&metric, &adm_aux);
  
  // This section sets up the initial prims, does not matter in this test case
  const int backup_routine[3] = {None,None,None};
  const bool evolve_entropy = false;
  const bool evolve_temperature = false;
  const bool calc_prims_guess = false;
  const double Psi6threshold = 1e100;
  const double W_max = 10.0;
  const double Lorenz_damping_factor = 0.0;

  const int neos = 1;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  double vx = 0.87;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, W_max, Lorenz_damping_factor, &params);

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);

  // Set prim
  // Ideal fluid P = (Gamma - 1) rho eps
  ghl_primitive_quantities prims;
  ghl_initialize_primitives(0.0, 0.0, 0.0, // rho p eps
                            vx , 0.0, 0.0, // vx vy vz
                            0.0, 0.0, 0.0, // Bx By Bz
                            0.0, 0.0, 0.0,   // ent Ye T
                        &prims);

  //This applies limits on the primitives
  bool speed_limited = false;
  ghl_enforce_primitive_limits_and_compute_u0(&params, &eos, &metric, &prims, &speed_limited);

  double E_thin = 3.0;
  ghl_radiation_flux_vector F4_thin;
  F4_thin.D[0] = 0.0;
  F4_thin.D[1] = E_thin;
  F4_thin.D[2] = 0.0;
  F4_thin.D[3] = 0.0;
  

  // For verifying test cases, these are not needed for source_update
  double u4U[4] = { prims.u0, prims.u0 * prims.vU[0], prims.u0 * prims.vU[1],
                          prims.u0 * prims.vU[2] };
  double n4D[4] = { -metric.lapse, 0, 0, 0 };
  double W = 0;
  ///////////////////////////////////////
  double cdt = 0.1;
  ghl_m1_closure_t closure = Minerbo;
  gsl_root_fsolver * gsl_solver_1d =
            gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_multiroot_fdfsolver * gsl_solver_nd =
            gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridsj, 4);
  ghl_m1_thc_params thc_params = {0};
  thc_params.source_thick_limit = 20.0;
  thc_params.source_scat_limit = -1.0;
  thc_params.source_maxiter = 64;
  thc_params.source_epsabs = 1.0e-15;
  thc_params.source_epsrel = 1.0e-5;

  
  // Initialize E and F4 that will be updated by ghl_source_update
  double E_new;
  double E_star_0 = E_thin;
  ghl_radiation_flux_vector F4_new;
  ghl_radiation_flux_vector F4_star_0 = {0};
  F4_star_0.D[0] = F4_thin.D[0];
  F4_star_0.D[1] = F4_thin.D[1];
  F4_star_0.D[2] = F4_thin.D[2];
  F4_star_0.D[3] = F4_thin.D[3];
  ghl_radiation_con_source_vector rF_source_0 = {0};

  // Initialize param struct to 0
  ghl_m1_powell_params p_TestZero     = {0}; 
  ghl_m1_powell_params p_TestExplicit = {0};
  ghl_m1_powell_params p_TestImplicit = {0};

  ///////////////////////////////// Test Zero case ////////////////////////////////
  // All source is zero
  double chi_TestZero = 0.0;
  double eta_TestZero = 0.0;
  double kabs_TestZero = 0.0;
  double kscat_TestZero = 0.0;

  // Initialize the parameters
  E_new = E_thin;
  F4_new.D[0] = F4_thin.D[0];
  F4_new.D[1] = F4_thin.D[1];
  F4_new.D[2] = F4_thin.D[2];
  F4_new.D[3] = F4_thin.D[3];

  printf("\nSource test: TestZero\n");
  // int ierr_TestZero = ghl_source_update(&thc_params, &p_TestZero, &E_new, &F4_new);
  int ierr_TestZero = ghl_source_update(
      &thc_params, chi_TestZero, eta_TestZero, kabs_TestZero, kscat_TestZero, cdt,
      &metric, &adm_aux, &prims,
      E_star_0, &F4_star_0,
      closure, gsl_solver_1d, gsl_solver_nd,
      &E_new, &F4_new);
  printf("\nResults:\n");
  printf("E_new: expected=%f  solution=%f  \n",E_thin, E_new);
  printf("F4_new: expected=(%f, %f, %f, %f)  solution=(%f, %f, %f, %f)\n",
         F4_thin.D[0], F4_thin.D[1], F4_thin.D[2], F4_thin.D[3],
         F4_new.D[0], F4_new.D[1], F4_new.D[2], F4_new.D[3]);
  /////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////// Test Explicit case ////////////////////////////////
  // cdt * kabs < 1 && cdt * kscat < 1
  // All source is zero
  double chi_TestExplicit = 0.0;
  double eta_TestExplicit = 0.5;
  double kabs_TestExplicit = 0.0;
  double kscat_TestExplicit = 0.0;

  // Initialize the parameters
  E_new = E_thin;
  F4_new.D[0] = F4_thin.D[0];
  F4_new.D[1] = F4_thin.D[1];
  F4_new.D[2] = F4_thin.D[2];
  F4_new.D[3] = F4_thin.D[3];

  printf("\nSource test: TestExplicit\n");
  int ierr_TestExplicit = ghl_source_update(
      &thc_params, chi_TestExplicit, eta_TestExplicit, kabs_TestExplicit, kscat_TestExplicit, cdt,
      &metric, &adm_aux, &prims,
      E_star_0, &F4_star_0,
      closure, gsl_solver_1d, gsl_solver_nd,
      &E_new, &F4_new);

  printf("\nResults:\n");
  double rE_source_0 = eta_TestExplicit * prims.u0;
  double E_TestExplicit = E_thin + cdt * rE_source_0;
  double Fx_TestExplicit = F4_thin.D[1] + cdt * rE_source_0 * vx;
  printf("E_new: expected=%f  solution=%f  \n",E_TestExplicit, E_new);
  printf("F4_new: expected=(%f, %f, %f, %f)  solution=(%f, %f, %f, %f)\n",
         0.0, Fx_TestExplicit, 0.0, 0.0,
         F4_new.D[0], F4_new.D[1], F4_new.D[2], F4_new.D[3]);
  /////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////// Test Implicit case ////////////////////////////////
  // not (cdt * kabs < 1 && cdt * kscat < 1)
  double chi_TestImplicit    = 19.0;
  double eta_TestImplicit    = 0.5;
  double kabs_TestImplicit   = 11;
  double kscat_TestImplicit  = 7.0;

  E_new = E_thin;
  F4_new.D[0] = F4_thin.D[0];
  F4_new.D[1] = F4_thin.D[1];
  F4_new.D[2] = F4_thin.D[2];
  F4_new.D[3] = F4_thin.D[3];

  printf("\nSource test: TestImplicit\n");
  int ierr_TestImplicit = ghl_source_update(
      &thc_params, chi_TestImplicit, eta_TestImplicit, kabs_TestImplicit, kscat_TestImplicit, cdt,
      &metric, &adm_aux, &prims,
      E_star_0, &F4_star_0,
      closure, gsl_solver_1d, gsl_solver_nd,
      &E_new, &F4_new);
  //TODO: This is only explicit limit, need to find a test case with exact implicit result
  double E_TestImplicit = E_new + cdt * eta_TestImplicit * prims.u0;
  double Fx_TestImplicit = F4_thin.D[1] + cdt * rE_source_0 * vx;
  printf("\nResults:\n");
  printf("E_new: expected=%f  solution=%f  \n",E_TestImplicit, E_new);
  printf("F4_new: expected=(%f, %f, %f, %f)  solution=(%f, %f, %f, %f)\n",
         0.0, Fx_TestImplicit, 0.0, 0.0,
         F4_new.D[0], F4_new.D[1], F4_new.D[2], F4_new.D[3]);
  /////////////////////////////////////////////////////////////////////////////////////
  }
  
