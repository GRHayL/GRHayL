#include "radiation.h"
#include "unit_tests.h"
#include "ghl_m1.h"

//copied from unit_test_nrpyleakage_optically_thin_gas.c
int main (int argc, char **argv) {
  printf("Run Optically thin test\n");
  test_optically_thin_advection();

  printf("Run Optically thick test\n");
  test_optically_thick_advection();
  return 0;
}

// Optically thick advection //TODO: this is not finished
void test_optically_thick_advection(){

  m1_root_params test_params;
  
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
                            0.87, 0.0, 0.0, // vx vy vz
                            0.0, 0.0, 0.0, // Bx By Bz
                            0.0, 0.0, 0.0,   // ent Ye T
                        &prims);

  //This applies limits on the primitives
  ghl_enforce_primitive_limits_and_compute_u0(&params, &eos, &metric, &prims);


  double E = 3.0;

  ghl_radiation_flux_vector F4;
  F4.D[0] = 0.5*E;
  F4.D[1] = -0.5*E;
  F4.D[2] = 0.0;
  F4.D[3] = 0.0;

  ghl_radiation_pressure_tensor P4;
  for (int a = 0; a < 4; a++){
    for (int b = 0; b < 4; b++){
      P4.DD[a][b] = 0.0;
    }
  }
  test_params.metric   = &metric;
  test_params.adm_aux  = &adm_aux;
  test_params.prims    = &prims;
  test_params.E        = E;
  test_params.F4       = &F4;
  test_params.P4       = &P4;
  ghl_radiation_rootSolve_closure(&test_params);
}

// Optically thin advection in Sec. 4.1
void test_optically_thin_advection(){

  m1_root_params test_params;
  
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
                            0.87, 0.0, 0.0, // vx vy vz
                            0.0, 0.0, 0.0, // Bx By Bz
                            0.0, 0.0, 0.0,   // ent Ye T
                        &prims);

  //This applies limits on the primitives
  ghl_enforce_primitive_limits_and_compute_u0(&params, &eos, &metric, &prims);


  double E = 3.0;

  ghl_radiation_flux_vector F4;
  F4.D[0] = 0.0;
  F4.D[1] = 0.99999999*E; //When set exactly to E, the root is found but not swapped
  F4.D[2] = 0.0;
  F4.D[3] = 0.0;

  ghl_radiation_pressure_tensor P4;
  for (int a = 0; a < 4; a++){
    for (int b = 0; b < 4; b++){
      P4.DD[a][b] = 0.0;
    }
  }
  test_params.metric   = &metric;
  test_params.adm_aux  = &adm_aux;
  test_params.prims    = &prims;
  test_params.E        = E;
  test_params.F4       = &F4;
  test_params.P4       = &P4;
  ghl_radiation_rootSolve_closure(&test_params);
}
