#include "ghl_unit_tests.h"
#include "ghl.h"
#include "ghl_m1.h"
#include "ghl_radiation.h"
#include <assert.h>

//copied from unit_test_nrpyleakage_optically_thin_gas.c
int main (int argc, char **argv) {
  
  //initialize the closure of choice
  ghl_error_codes_t initialize_check = ghl_initialize_m1_closure(Minerbo);
  ghl_read_error_codes(initialize_check);

  printf("\nRun Calc update test\n");
  test_calc_update();
  printf("'unit_test_m1_calc_update' assert passed\n");

  return 0;
}

int test_calc_update(){
  // Set metric quantities to Minkowski
  ghl_metric_quantities metric;
  ghl_initialize_metric(1, 0, 0, 0, //lapse, shift
                           1, 0, 0, // gxx gxy gxz
                              1, 0, //     gyy gyz
                                 1, //         gzz
                           &metric);

  ghl_ADM_aux_quantities adm_aux;
  ghl_compute_ADM_auxiliaries(&metric, &adm_aux);

  ghl_metric_quantities metric_derivs_x = {0};
  ghl_metric_quantities metric_derivs_y = {0};
  ghl_metric_quantities metric_derivs_z = {0};
  metric_derivs_x.lapse = 0.5;
  metric_derivs_y.lapse = 0.7;
  metric_derivs_z.lapse = 0.9;
  ghl_extrinsic_curvature curv = {0};
  curv.K[0][0] = 1.3;

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
  F4_thin.D[1] = E_thin; //When set exactly to E, the root is found but not swapped
  F4_thin.D[2] = 0.0;
  F4_thin.D[3] = 0.0;
  
  
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
  ghl_radiation_flux_vector F4_new;
  F4_new.D[0] = F4_thin.D[0];
  F4_new.D[1] = F4_thin.D[1];
  F4_new.D[2] = F4_thin.D[2];
  F4_new.D[3] = F4_thin.D[3];

  #define NGROUPS 1
  #define NSPECIES 1
  ghl_m1_conservative_quantities rU = {0};
  ghl_m1_conservative_quantities rU_p = {0};
  int ig = 0;
  rU.rE[ig] = E_thin;
  rU.rF4[ig].D[0] = F4_thin.D[0];
  rU.rF4[ig].D[1] = F4_thin.D[1];
  rU.rF4[ig].D[2] = F4_thin.D[2];
  rU.rF4[ig].D[3] = F4_thin.D[3];
  rU.rN[ig] = 0.0;
  rU_p.rE[ig] = E_thin;
  rU_p.rF4[ig].D[0] = F4_thin.D[0];
  rU_p.rF4[ig].D[1] = F4_thin.D[1];
  rU_p.rF4[ig].D[2] = F4_thin.D[2];
  rU_p.rF4[ig].D[3] = F4_thin.D[3];
  rU_p.rN[ig] = 0.0;
  
  double abs_0[NGROUPS*NSPECIES]    = {0};
  double abs_1[NGROUPS*NSPECIES]    = {0};
  double eta_0[NGROUPS*NSPECIES]    = {0};
  double eta_1[NGROUPS*NSPECIES]    = {0};
  double scat_0[NGROUPS*NSPECIES]   = {0};
  double scat_1[NGROUPS*NSPECIES]   = {0};

  // Assign test case here
  // non-zero test values for absorption (abs), emission (eta) and scattering (scat)
  abs_0[0] = 0.0;
  abs_1[0] = 0.0;
  eta_0[0] = 0.0;
  eta_1[0] = 0.1;
  scat_0[0] = 0.0;
  scat_1[0] = 0.0;
  
  ghl_m1_conservative_quantities rhs_flux = {0};
  ghl_m1_conservative_quantities rhs_gr = {0};
  
  ghl_m1_update(
    cdt,
    &metric,
    &curv,
    &adm_aux,
    &prims,
    &thc_params,
    &abs_0,
    &abs_1,
    &eta_0,
    &eta_1,
    &scat_0,
    &scat_1,
    &rU,
    &rU_p,
    &rhs_flux,
    &rhs_gr
  );

  /* 
   * The expected update from the source term is given by the following equations:
   * E_new   = E_old   + cdt * S^0 = E_old   + cdt * eta * gamma
   * F^x_new = F^x_old + cdt * S^1 = F^x_old + cdt * eta * gamma * vx
   */
  const double gamma = 1.0 / sqrt(1.0 - vx * vx);
  const double expected_E_from_source  = rU_p.rE[ig]       + cdt * eta_1[0] * gamma;
  const double expected_Fx_from_source = rU_p.rF4[ig].D[1] + cdt * eta_1[0] * gamma * vx;
  const double tol_E = 1.0e-6;
  const double tol_Fx = 1.0e-6;
  assert(fabs(rU.rE[ig] - expected_E_from_source) < tol_E);
  assert(fabs(rU.rF4[ig].D[1] - expected_Fx_from_source) < tol_Fx);
} // test_calc_update
