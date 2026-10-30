#include "ghl_unit_tests.h"
#include "ghl_radiation.h"
#include "ghl_m1.h"
#include <assert.h>
void test_advection_grid();

//copied from unit_test_nrpyleakage_optically_thin_gas.c
int main (int argc, char **argv) {
  
  // thin-limit example: a small x-directed perturbation on a uniform background
  ghl_error_codes_t initialize_check = ghl_initialize_m1_closure(Thin);
  ghl_read_error_codes(initialize_check);

  printf("Run thin-limit M1 advection grid test\n");
  test_advection_grid();

  return 0;
}


void test_advection_grid() {
  // 1D advection test using ghl_M1_update_grid (based on test 4.1)
  // Keep the transverse directions uniform so this is a single x-directed advection mode.
  
  
  // Set metric quantities to Minkowski
    ghl_metric_quantities metric;
    ghl_initialize_metric(1, 0, 0, 0, //lapse, shift
                    1, 0, 0, // gxx gxy gxz
                      1, 0, //     gyy gyz
                        1, //         gzz
                    &metric);

    ghl_ADM_aux_quantities adm_aux;
    ghl_extrinsic_curvature curv;
    ghl_compute_ADM_auxiliaries(&metric, &adm_aux);
    printf("Minkowski metric initialized: lapse=%g, betaU[0]=%g, betaU[1]=%g, betaU[2]=%g betaU[3]=%g, gammaDD[0][0]=%g, gammaDD[0][1]=%g, gammaDD[0][2]=%g, gammaDD[1][1]=%g, gammaDD[1][2]=%g, gammaDD[2][2]=%g\n", metric.lapse, metric.betaU[0], metric.betaU[1], metric.betaU[2], metric.betaU[3], metric.gammaDD[0][0], metric.gammaDD[0][1], metric.gammaDD[0][2], metric.gammaDD[1][1], metric.gammaDD[1][2], metric.gammaDD[2][2]);
    ghl_initialize_extrinsic_curvature(
        0, 0, 0,
        0, 0, 0,
        &curv);

  const int grid_Nx = 3, grid_Ny = 3, grid_Nz = 200;
  const double h = 0.01;
  double cdt = 0.3 * h; // timestep
  printf("gridsize h = %e, timestep = %e\n", h, cdt);
  
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

  ///////////////////////////////////////
  
  ghl_m1_thc_params thc_params;
  thc_params.source_thick_limit = 20.0;
  thc_params.source_scat_limit = -1.0;
  thc_params.source_maxiter = 64;
  thc_params.source_epsabs = 1.0e-15;
  thc_params.source_epsrel = 1.0e-5;
  thc_params.minmod_theta = 1.0;
  thc_params.mindiss = 1.0;

  #define NGROUPS 1
  #define NSPECIES 1

  const int array_length = NGROUPS * NSPECIES * grid_Nx * grid_Ny * grid_Nz;

  double *E_grid = (double*) calloc(array_length, sizeof(double));
  ghl_radiation_flux_vector *F4_grid = (ghl_radiation_flux_vector*) calloc(array_length, sizeof(ghl_radiation_flux_vector));
  double *N_grid = (double*) calloc(array_length, sizeof(double));


  

  // Add IDX macros after testing
  #define IDX_GRID(ig, ix, iy, iz) ((ig)*NGROUPS*NSPECIES*grid_Nx*grid_Ny*grid_Nz \
                                  + (ix)*NGROUPS*NSPECIES*grid_Ny*grid_Nz \
                                  + (iy)*NGROUPS*NSPECIES*grid_Nz \
                                  + (iz))

  // initialize a smooth x-advecting pulse on top of a nonzero background
  // thin-limit behavior is reached by setting the closure above and keeping
  // the flow nearly aligned with the radiation flux direction.

  const double dx = h, dy = h, dz = h;
  const double bg_E = 3.0;
  const double pert_amp = 0.01;
  const double x_center = -0.1;
  const double y_center = 0.0;
  const double z_center = 0.0;
  const double xmin = -0.5 * (grid_Nx - 1) * dx;
  const double xmax =  0.5 * (grid_Nx - 1) * dx;
  const double ymin = -0.5 * (grid_Ny - 1) * dy;
  const double ymax =  0.5 * (grid_Ny - 1) * dy;
  const double zmin = -0.5 * (grid_Nz - 1) * dz;
  const double zmax =  0.5 * (grid_Nz - 1) * dz;
  const double sigma_x = 0.1;

  for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
    for (int ix = 0; ix < grid_Nx; ++ix) {
      for (int iy = 0; iy < grid_Ny; ++iy) {
        for (int iz = 0; iz < grid_Nz; ++iz) {
          int idx = IDX_GRID(ig, ix, iy, iz);
          double x = (ix - 0.5 * (grid_Nx - 1)) * dx;
          double y = (iy - 0.5 * (grid_Ny - 1)) * dy;
          double z = (iz - 0.5 * (grid_Nz - 1)) * dz;
          double E0 = (z >= -0.5) ? 1.0 : 0.0;

          E_grid[idx] = E0;
          F4_grid[idx].D[1] = 0.0;  // Fx
          F4_grid[idx].D[2] = 0.0;  // Fy
          F4_grid[idx].D[3] = E0;   // Fz = E

          N_grid[idx] = 1.0;         // A simple choice; N is not specified by this test
          get_FU(&adm_aux, &F4_grid[idx]);
        }
      }
    }
  }

  // build metric/prims grids
  ghl_metric_quantities *metric_grid = (ghl_metric_quantities*) calloc(grid_Nx*grid_Ny*grid_Nz, sizeof(ghl_metric_quantities));
  ghl_ADM_aux_quantities *adm_aux_grid = (ghl_ADM_aux_quantities*) calloc(grid_Nx*grid_Ny*grid_Nz, sizeof(ghl_ADM_aux_quantities));
  ghl_extrinsic_curvature *curv_grid = (ghl_extrinsic_curvature*) calloc(grid_Nx*grid_Ny*grid_Nz, sizeof(ghl_extrinsic_curvature));
  ghl_primitive_quantities *prims_grid = (ghl_primitive_quantities*) calloc(grid_Nx*grid_Ny*grid_Nz, sizeof(ghl_primitive_quantities));

  for (int ix=0;ix<grid_Nx;ix++){
    for (int iy=0;iy<grid_Ny;iy++){
      for (int iz=0;iz<grid_Nz;iz++){
        int idx = ix*grid_Ny*grid_Nz + iy*grid_Nz + iz;
        double z = (iz - 0.5 * (grid_Nz - 1)) * dz;
        double vz = (z < 0.0) ? 0.87 : -0.87;
        ghl_primitive_quantities prims;
        ghl_initialize_primitives(0.0, 0.0, 0.0,
                                  0.0, 0.0, vz,   // vx, vy, vz
                                  0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0,
                                  &prims);

         //This applies limits on the primitives
        bool speed_limited = false;
        ghl_enforce_primitive_limits_and_compute_u0(&params, &eos, &metric, &prims, &speed_limited);
        metric_grid[idx] = metric;
        adm_aux_grid[idx] = adm_aux;
        curv_grid[idx] = curv;
        prims_grid[idx] = prims;
      }
    }
  }

  // Prepare opacity grids (zero)
  double *abs_0_grid  = (double*) calloc(NGROUPS*NSPECIES*grid_Nx*grid_Ny*grid_Nz, sizeof(double));
  double *abs_1_grid  = (double*) calloc(NGROUPS*NSPECIES*grid_Nx*grid_Ny*grid_Nz, sizeof(double));
  double *eta_0_grid  = (double*) calloc(NGROUPS*NSPECIES*grid_Nx*grid_Ny*grid_Nz, sizeof(double));
  double *eta_1_grid  = (double*) calloc(NGROUPS*NSPECIES*grid_Nx*grid_Ny*grid_Nz, sizeof(double));
  double *scat_0_grid = (double*) calloc(NGROUPS*NSPECIES*grid_Nx*grid_Ny*grid_Nz, sizeof(double));
  double *scat_1_grid = (double*) calloc(NGROUPS*NSPECIES*grid_Nx*grid_Ny*grid_Nz, sizeof(double));
  for (int ig=0; ig<NGROUPS*NSPECIES; ++ig){
    for (int ix=0; ix<grid_Nx; ++ix){
      for (int iy=0; iy<grid_Ny; ++iy){
        for (int iz=0; iz<grid_Nz; ++iz){
          int idx = IDX_GRID(ig, ix, iy, iz);
          abs_0_grid[idx] = 0.0;
          abs_1_grid[idx] = 0.0;
          eta_0_grid[idx] = 0.0;
          eta_1_grid[idx] = 0.0;
          scat_0_grid[idx] = 0.0;
          scat_1_grid[idx] = 0.0;
        }
      }
    }
  }



  // Prepare conservative grids
  ghl_m1_conservative_quantities *rU_grid = (ghl_m1_conservative_quantities*) calloc(array_length, sizeof(ghl_m1_conservative_quantities));
  ghl_m1_conservative_quantities *rU_p_grid = (ghl_m1_conservative_quantities*) calloc(array_length, sizeof(ghl_m1_conservative_quantities));
  ghl_m1_conservative_quantities *rU_rhs_grid = (ghl_m1_conservative_quantities*) calloc(array_length, sizeof(ghl_m1_conservative_quantities));

  // Prepare primitive/conservative grids per cell with the same x-direction flow in every cell.
  for (int ig=0; ig<NGROUPS*NSPECIES; ++ig){
    for (int ix=0; ix<grid_Nx; ++ix){
      for (int iy=0; iy<grid_Ny; ++iy){
        for (int iz=0; iz<grid_Nz; ++iz){
          int idx = IDX_GRID(ig, ix, iy, iz);
          rU_p_grid[idx].rE[ig] = E_grid[idx];
          rU_p_grid[idx].rF4[ig].D[1] = F4_grid[idx].D[1];
          rU_p_grid[idx].rF4[ig].D[2] = F4_grid[idx].D[2];
          rU_p_grid[idx].rF4[ig].D[3] = F4_grid[idx].D[3];
          rU_p_grid[idx].rN[ig] = N_grid[idx];
          // initialize current with same values
          rU_grid[idx] = rU_p_grid[idx];
        }
      }
    }
  }

  // Run several update steps and write the requested 1D slice files.
  // x-slice keeps y = z = 0, y-slice keeps x = z = 0, z-slice keeps x = y = 0.
  double *rhs_grid = (double*) calloc(grid_Nx*grid_Ny*grid_Nz*5, sizeof(double));
  FILE *fp_x = fopen("/media/miguel/localdata/bt22824/simulations/m1_advection_thin_x.asc", "w");
  FILE *fp_y = fopen("/media/miguel/localdata/bt22824/simulations/m1_advection_thin_y.asc", "w");
  FILE *fp_z = fopen("/media/miguel/localdata/bt22824/simulations/m1_advection_thin_z.asc", "w");
  if(fp_x == NULL || fp_y == NULL || fp_z == NULL) {
    fprintf(stderr, "Failed to open one or more slice output files for writing\n");
    exit(1);
  }
  const double final_time = 1.0;
  const int nsteps = (int) (final_time / cdt);
  // const int nsteps = 4;
  printf("Total time steps = %d\n", nsteps);
  for(int istep = 0; istep < nsteps; ++istep) {
    printf("\rProgress: %3f%%", istep / (double)nsteps * 100);
    fflush(stdout);
    
    const double time = istep * cdt;
    const int ix_center = (grid_Nx - 1) / 2;
    const int iy_center = (grid_Ny - 1) / 2;
    const int iz_center = (grid_Nz - 1) / 2;

    for(int ix = 0; ix < grid_Nx; ++ix) {
      const int idx = IDX_GRID(0, ix, iy_center, iz_center);
      const double x = (ix - 0.5 * (grid_Nx - 1)) * dx;
      const double E = rU_grid[idx].rE[0];
      const double Fx = rU_grid[idx].rF4[0].D[1];
      const double Fy = rU_grid[idx].rF4[0].D[2];
      const double Fz = rU_grid[idx].rF4[0].D[3];
      const double N = rU_grid[idx].rN[0];
      fprintf(fp_x, "%g %g %g %g %g %g %g\n",
              time, x, E, Fx, Fy, Fz, N);
    }
    fprintf(fp_x, "\n");

    for(int iy = 0; iy < grid_Ny; ++iy) {
      const int idx = IDX_GRID(0, ix_center, iy, iz_center);
      const double y = (iy - 0.5 * (grid_Ny - 1)) * dy;
      const double E = rU_grid[idx].rE[0];
      const double Fx = rU_grid[idx].rF4[0].D[1];
      const double Fy = rU_grid[idx].rF4[0].D[2];
      const double Fz = rU_grid[idx].rF4[0].D[3];
      const double N = rU_grid[idx].rN[0];
      fprintf(fp_y, "%g %g %g %g %g %g %g\n",
              time, y, E, Fx, Fy, Fz, N);
    }
    fprintf(fp_y, "\n");

    for(int iz = 0; iz < grid_Nz; ++iz) {
      const int idx = IDX_GRID(0, ix_center, iy_center, iz);
      const double z = (iz - 0.5 * (grid_Nz - 1)) * dz;
      const double E = rU_grid[idx].rE[0];
      const double Fx = rU_grid[idx].rF4[0].D[1];
      const double Fy = rU_grid[idx].rF4[0].D[2];
      const double Fz = rU_grid[idx].rF4[0].D[3];
      const double N = rU_grid[idx].rN[0];
      fprintf(fp_z, "%g %g %g %g %g %g %g\n",
              time, z, E, Fx, Fy, Fz, N);
    }
    fprintf(fp_z, "\n");

    if(istep + 1 < nsteps) {
      for(int iy = 0; iy < grid_Ny; ++iy) {
        for(int iz = 0; iz < grid_Nz; ++iz) {
          for(int ix = 0; ix < grid_Nx; ++ix) {
            const int idx = IDX_GRID(0, ix, iy, iz);
            E_grid[idx] = rU_grid[idx].rE[0];
            F4_grid[idx].D[1] = rU_grid[idx].rF4[0].D[1];
            F4_grid[idx].D[2] = rU_grid[idx].rF4[0].D[2];
            F4_grid[idx].D[3] = rU_grid[idx].rF4[0].D[3];
            N_grid[idx] = rU_grid[idx].rN[0];
            rU_p_grid[idx] = rU_grid[idx];
          }
        }
      }

      int ierr = ghl_m1_update_grid(
            grid_Nx, grid_Ny, grid_Nz,
            0, grid_Nx-1, 0, grid_Ny-1, 0, grid_Nz-1,
            cdt,
            &thc_params,
            (double[]){dx, dy, dz},
            metric_grid,
            curv_grid,
            adm_aux_grid,
            prims_grid,
            // E_grid,
            // F4_grid,
            // N_grid,
            abs_0_grid,
            abs_1_grid,
            eta_0_grid,
            eta_1_grid,
            scat_0_grid,
            scat_1_grid,
            rU_grid,
            rU_p_grid);
    } // if(istep + 1 < nsteps)
  }// for istep
  printf("\n");
  fclose(fp_x);
  fclose(fp_y);
  fclose(fp_z);

  // free allocated memory
  free(E_grid);
  free(F4_grid);
  free(N_grid);
  free(metric_grid);
  free(adm_aux_grid);
  free(curv_grid);
  free(prims_grid);
  free(abs_0_grid);
  free(abs_1_grid);
  free(eta_0_grid);
  free(eta_1_grid);
  free(scat_0_grid);
  free(scat_1_grid);
  free(rU_grid);
  free(rU_p_grid);
  free(rU_rhs_grid);
  free(rhs_grid);
}