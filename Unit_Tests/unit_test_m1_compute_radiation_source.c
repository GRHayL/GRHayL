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

  // Initialize the pressure tensor to 0, will be updated at prepare_closure
  ghl_radiation_pressure_tensor P4_0;
  for (int a = 0; a < 4; a++){
    for (int b = 0; b < 4; b++){
      P4_0.DD[a][b] = 0.0;
    }
  }
  // Initialize E and F4 that will be updated by ghl_source_update
  double E_new;
  ghl_radiation_flux_vector F4_new;

  ghl_radiation_flux_vector F4_star_0;
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
  p_TestZero.chi    = 0.0; // initial chi value
  p_TestZero.eta    = 0.0;
  p_TestZero.kabs   = 0.0;
  p_TestZero.kscat  = 0.0;
  
  p_TestZero.E              = E_thin;
  p_TestZero.F4             = &F4_thin;
  p_TestZero.P4             = &P4_0;
  p_TestZero.metric         = &metric;
  p_TestZero.adm_aux        = &adm_aux;
  p_TestZero.prims          = &prims;
  p_TestZero.rF_source      = &rF_source_0;
  // p_TestZero.E_new          = E_thin;
  // p_TestZero.F4_new         = &F4_new_0;
  p_TestZero.E_star         = E_thin;
  p_TestZero.F4_star        = &F4_star_0;
  p_TestZero.cdt            = cdt;
  p_TestZero.closure        = closure;
  p_TestZero.gsl_solver_1d  = gsl_solver_1d;
  p_TestZero.gsl_solver_nd  = gsl_solver_nd;

  E_new = E_thin;
  F4_new.D[0] = F4_thin.D[0];
  F4_new.D[1] = F4_thin.D[1];
  F4_new.D[2] = F4_thin.D[2];
  F4_new.D[3] = F4_thin.D[3];

  printf("\nSource test: TestZero\n");
  int ierr_TestZero = ghl_source_update(&thc_params, &p_TestZero, &E_new, &F4_new);
  printf("E_new: expected=%f  solution=%f  \n",E_thin, E_new);
  printf("F4_new: expected=(%f, %f, %f, %f)  solution=(%f, %f, %f, %f)\n",
         F4_thin.D[0], F4_thin.D[1], F4_thin.D[2], F4_thin.D[3],
         F4_new.D[0], F4_new.D[1], F4_new.D[2], F4_new.D[3]);
  /////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////// Test Explicit case ////////////////////////////////
  // cdt * kabs < 1 && cdt * kscat < 1
  p_TestExplicit.chi    = 0.0;
  p_TestExplicit.eta    = 0.5;
  p_TestExplicit.kabs   = 0.0;
  p_TestExplicit.kscat  = 0.0;
  p_TestExplicit.E              = E_thin;
  p_TestExplicit.F4             = &F4_thin;
  p_TestExplicit.P4             = &P4_0;
  p_TestExplicit.metric         = &metric;
  p_TestExplicit.adm_aux        = &adm_aux;
  p_TestExplicit.prims          = &prims;
  p_TestExplicit.E_star         = E_thin;
  p_TestExplicit.F4_star        = &F4_star_0;
  p_TestExplicit.rF_source      = &rF_source_0;
  p_TestExplicit.cdt            = cdt;
  p_TestExplicit.closure        = closure;
  p_TestExplicit.gsl_solver_1d  = gsl_solver_1d;
  p_TestExplicit.gsl_solver_nd  = gsl_solver_nd;

  E_new = E_thin;
  F4_new.D[0] = F4_thin.D[0];
  F4_new.D[1] = F4_thin.D[1];
  F4_new.D[2] = F4_thin.D[2];
  F4_new.D[3] = F4_thin.D[3];

  printf("\nSource test: TestImplicit\n");
  int ierr_TestExplicit = ghl_source_update(&thc_params, &p_TestExplicit, &E_new, &F4_new);
  // Eq (29)-0 rE_source = - alp * volfrom * (S^mu n_mu) = - eta * u^0
  double rE_source_0 = + p_TestExplicit.eta * p_TestExplicit.prims->u0;
  double E_TestExplicit = E_thin + p_TestExplicit.cdt * rE_source_0;
  printf("E_new: expected=%f  solution=%f  \n",E_TestExplicit, E_new);
  printf("F4_new: expected=(%f, %f, %f, %f)  solution=(%f, %f, %f, %f)\n",
         F4_thin.D[0], F4_thin.D[1], F4_thin.D[2], F4_thin.D[3],
         F4_new.D[0], F4_new.D[1], F4_new.D[2], F4_new.D[3]);
  /////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////// Test Implicit case ////////////////////////////////
  // not (cdt * kabs < 1 && cdt * kscat < 1)
  p_TestImplicit.chi    = 0.0;
  p_TestImplicit.eta    = 120.0;
  p_TestImplicit.kabs   = 17.0;
  p_TestImplicit.kscat  = 19.0;
  p_TestImplicit.E              = E_thin;
  p_TestImplicit.F4             = &F4_thin;
  p_TestImplicit.P4             = &P4_0;
  p_TestImplicit.metric         = &metric;
  p_TestImplicit.adm_aux        = &adm_aux;
  p_TestImplicit.prims          = &prims;
  p_TestImplicit.E_star         = E_thin;
  p_TestImplicit.F4_star        = &F4_star_0;
  p_TestImplicit.rF_source      = &rF_source_0;
  p_TestImplicit.cdt            = cdt;
  p_TestImplicit.closure        = closure;
  p_TestImplicit.gsl_solver_1d  = gsl_solver_1d;
  p_TestImplicit.gsl_solver_nd  = gsl_solver_nd;

  E_new = E_thin;
  F4_new.D[0] = F4_thin.D[0];
  F4_new.D[1] = F4_thin.D[1];
  F4_new.D[2] = F4_thin.D[2];
  F4_new.D[3] = F4_thin.D[3];

  printf("\nSource test: TestExplicit\n");
  int ierr_TestImplicit = ghl_source_update(&thc_params, &p_TestImplicit, &E_new, &F4_new);
  //TODO: This is only explicit limit, need to find a test case with exact implicit result
  double E_TestImplicit = E_new + p_TestImplicit.cdt * p_TestImplicit.eta * p_TestImplicit.prims->u0;
  printf("E_new: expected=%f  solution=%f  \n",E_TestImplicit, E_new);
  printf("F4_new: expected=(%f, %f, %f, %f)  solution=(%f, %f, %f, %f)\n",
         F4_thin.D[0], F4_thin.D[1], F4_thin.D[2], F4_thin.D[3],
         F4_new.D[0], F4_new.D[1], F4_new.D[2], F4_new.D[3]);
  /////////////////////////////////////////////////////////////////////////////////////


  }
  
