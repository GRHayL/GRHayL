#include "ghl_unit_tests.h"
#include "ghl_radiation.h"
#include "ghl_m1.h"
#include <assert.h>
void test_flux_pointwise();

//copied from unit_test_nrpyleakage_optically_thin_gas.c
int main (int argc, char **argv) {
  
  //initialize the closure of choice
  ghl_error_codes_t initialize_check = ghl_initialize_m1_closure(Minerbo);
  ghl_read_error_codes(initialize_check);

  printf("\nRun M1 Flux test\n");
  test_flux_pointwise();
  printf("'unit_test_m1_compute_radiation_flux' assert passed\n");

  return 0;
}

void test_flux_pointwise() {
  int grid_Nx = 10, grid_Ny = 10, grid_Nz = 10;
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

  

  // For verifying test cases, these are not needed for source_update
  // double u4U[4] = { prims.u0, prims.u0 * prims.vU[0], prims.u0 * prims.vU[1],
  //                         prims.u0 * prims.vU[2] };
  // double n4D[4] = { -metric.lapse, 0, 0, 0 };
  // double W = 0;

  ///////////////////////////////////////
  double cdt = 0.1;
  ghl_m1_thc_params thc_params;
  thc_params.source_thick_limit = 20.0;
  thc_params.source_scat_limit = -1.0;
  thc_params.source_maxiter = 64;
  thc_params.source_epsabs = 1.0e-15;
  thc_params.source_epsrel = 1.0e-5;

  #define NGROUPS 1
  #define NSPECIES 1

  #define NSTENCIL_CONS 5
  #define NSTENCIL_FLUX 3

  #define STENCIL_LENGTH (NGROUPS * NSPECIES * NSTENCIL_CONS)
  int array_length = NGROUPS * NSPECIES * grid_Nx * grid_Ny * grid_Nz;
  double stencil_rE[STENCIL_LENGTH]                         = {0};
  ghl_radiation_flux_vector stencil_rF4[STENCIL_LENGTH]     = {0};
  ghl_radiation_pressure_tensor stencil_rP4[STENCIL_LENGTH] = {0};
  double stencil_rN[STENCIL_LENGTH]                         = {0};
  ghl_metric_quantities stencil_metric[NSTENCIL_FLUX]      = {0};
  ghl_ADM_aux_quantities stencil_adm_aux[NSTENCIL_FLUX]    = {0};
  ghl_primitive_quantities stencil_prims[NSTENCIL_FLUX]    = {0};


  

  // Add IDX macros after testing
  #define IDX_GRID(ig, ix, iy, iz) ((ig)*NGROUPS*NSPECIES*grid_Nz*grid_Ny*grid_Nx \
                                  + (ix)*NGROUPS*NSPECIES*grid_Nz*grid_Ny \
                                  + (iy)*NGROUPS*NSPECIES*grid_Nz \
                                  + (iz))

  ///////////////////////////////// Test Zero case ////////////////////////////////
  // All source is zero
  double kabs_TestZero = 5.0;
  double kscat_TestZero = 7.0;

  double E_thin = 3.0;
  ghl_radiation_flux_vector F4_thin;
  F4_thin.D[0] = 0.0;
  F4_thin.D[1] = E_thin;
  F4_thin.D[2] = E_thin;
  F4_thin.D[3] = 0.0;
  ghl_radiation_pressure_tensor P4_0;
  for (int a = 0; a < 4; a++){
    for (int b = 0; b < 4; b++){
      P4_0.DD[a][b] = 0.0;
    }
  }

  // variable length arrays cannot be initialized
  printf(" - Pointwise flux test: Constant gradient\n");
  double E_grid[array_length];
  ghl_radiation_flux_vector F4_grid[array_length];
  ghl_radiation_pressure_tensor P4_grid[array_length];
  double N_grid[array_length];
  double F4xdx = 0.3, F4ydx = 0.5, F4zdx = 0.7;
  double P4xxdx = 0.9, P4xydx = 0.11, P4xzdx = 0.13;

  for(int ig=0; ig<NGROUPS*NSPECIES; ig++) {
    for(int iz=2; iz<grid_Nz-2; iz++) {
      for(int iy=2; iy<grid_Ny-2; iy++) {
        for(int ix=2; ix<grid_Nx-2; ix++) {
          E_grid[IDX_GRID(ig, ix, iy, iz)] = E_thin;
          F4_grid[IDX_GRID(ig, ix, iy, iz)] = F4_thin;
          F4_grid[IDX_GRID(ig, ix, iy, iz)].D[1] += ix * F4xdx; // F4x varies with x
          F4_grid[IDX_GRID(ig, ix, iy, iz)].D[2] += ix * F4ydx; // F4y varies with x
          F4_grid[IDX_GRID(ig, ix, iy, iz)].D[3] += ix * F4zdx; // F4z varies with x
          P4_grid[IDX_GRID(ig, ix, iy, iz)] = P4_0;
          P4_grid[IDX_GRID(ig, ix, iy, iz)].DD[1][1] += P4xxdx * ix; // P4xx varies with x
          P4_grid[IDX_GRID(ig, ix, iy, iz)].DD[1][2] += P4xydx * ix; // P4xy varies with x
          P4_grid[IDX_GRID(ig, ix, iy, iz)].DD[1][3] += P4xzdx * ix; // P4xz varies with x
          N_grid[IDX_GRID(ig, ix, iy, iz)] = 0.0;
        }
      }
    }
  }

  // initialize stencil arrays
  int inspect_ix = 5;
  int inspect_iy = 5;
  int inspect_iz = 5;
  for(int ig=0; ig<NGROUPS*NSPECIES; ig++){
    for(int k=0; k<NSTENCIL_CONS; k++){
      int idx = ig*NSTENCIL_CONS + k;
      // Very Important: flux uses P4.UD and F4.U, we need to initialize them
      stencil_rE[idx]  =  E_grid[IDX_GRID(ig, inspect_ix + k - 1, inspect_iy, inspect_iz)];
      stencil_rF4[idx] = F4_grid[IDX_GRID(ig, inspect_ix + k - 1, inspect_iy, inspect_iz)];
      stencil_rP4[idx] = P4_grid[IDX_GRID(ig, inspect_ix + k - 1, inspect_iy, inspect_iz)];
      stencil_rN[idx]  =  N_grid[IDX_GRID(ig, inspect_ix + k - 1, inspect_iy, inspect_iz)];
      get_FU(&adm_aux, &stencil_rF4[idx]);
      get_PUD_and_PUU(&adm_aux,&stencil_rP4[idx]);
    }
  }
  for (int k=0; k<NSTENCIL_FLUX; k++){
    stencil_metric[k]   = metric;
    stencil_adm_aux[k]  = adm_aux;
    stencil_prims[k]    = prims;
  }
  double h = 0.1;
  double dx = h, dy = h, dz = h;
  double delta[3] = {dx,dy,dz};
  int dir = 0; // x-direction
  double abs_1[NGROUPS*NSPECIES] = {0.0};
  double scat_1[NGROUPS*NSPECIES] = {0.0};
  for (int ig=0; ig<NGROUPS*NSPECIES; ig++){
    abs_1[ig] = kabs_TestZero;
    scat_1[ig] = kscat_TestZero;
  }

  // double* rhs = (double*) malloc(sizeof(double)*NGROUPS*NSPECIES*5);
  ghl_m1_conservative_quantities rhs_flux = {0};
  ghl_m1_conservative_quantities rhs_flux_expected = {0};
  double rhs_expected[NGROUPS*NSPECIES*5] = {0};
  for (int ig=0; ig<NGROUPS*NSPECIES; ig++){
    rhs_flux_expected.rN[ig] = 0.0;
    rhs_flux_expected.rF4[ig].D[1] = -P4xxdx / dx;
    rhs_flux_expected.rF4[ig].D[2] = -P4xydx / dx;
    rhs_flux_expected.rF4[ig].D[3] = -P4xzdx / dx;
    rhs_flux_expected.rE[ig] = -F4xdx  / dx;
  }
  
  int ierr_TestZero = ghl_m1_calc_flux(
    cdt,
    stencil_metric,
    stencil_adm_aux,
    stencil_prims,
    &thc_params,
    stencil_rE,
    stencil_rF4,
    stencil_rP4,
    stencil_rN,
    delta,
    abs_1,
    scat_1,
    dir,
    &rhs_flux
  );

  const double rhs_tol = 1e-13;
  for(int ig=0; ig<NGROUPS*NSPECIES; ig++){
    assert(fabs(rhs_flux.rN[ig]       - rhs_flux_expected.rN[ig]) < rhs_tol);
    assert(fabs(rhs_flux.rF4[ig].D[1] - rhs_flux_expected.rF4[ig].D[1]) < rhs_tol);
    assert(fabs(rhs_flux.rF4[ig].D[2] - rhs_flux_expected.rF4[ig].D[2]) < rhs_tol);
    assert(fabs(rhs_flux.rF4[ig].D[3] - rhs_flux_expected.rF4[ig].D[3]) < rhs_tol);
    assert(fabs(rhs_flux.rE[ig]       - rhs_flux_expected.rE[ig]) < rhs_tol);
  }
} // test_flux_pointwise

  
