#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
#include "GRHayL.h"
#include "con2prim.h"
/*
 * // main() function:
 * // Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
 * // Step 1: Write test data to gridfunctions
 * // Step 2: Overwrite all data in ghost zones with NaNs
 * // Step 3: Apply boundary conditions
 * // Step 4: Print gridfunction data after curvilinear boundary conditions have been applied
 * // Step 5: Free all allocated memory
 */
int main(int argc, const char *argv[]) {

    
    
     /* ################################################
                     Initialize GRHayL
     ################################################ */

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  const int update_Tmunu = 1; //IGM default

  const int neos = 1;
  const double W_max = 10.0; //IGM default
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300; //IGM default
  const double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(Noble2D, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess, Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);


  /* ################################################
                     Initialize GRHayL
     ################################################ */
    
    
  griddata_struct griddata;
  set_Cparameters_to_default(&griddata.params);
  griddata.params.outer_bc_type = RADIATION_OUTER_BCS;

  // Step 0a: Read command-line input, error out if nonconformant
  if(argc != 5 || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < NGHOSTS) {
    printf("Error: Expected one command-line argument: ./ScalarWaveCurvilinear_Playground Nx0 Nx1 Nx2,\n");
    printf("where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
    printf("Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
    exit(1);
  }
  // Step 0b: Set up numerical grid structure, first in space...
  const int Nxx[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };
  if(Nxx[0]%2 != 0 || Nxx[1]%2 != 0 || Nxx[2]%2 != 0) {
    printf("Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
    printf("       For example, in case of angular directions, proper symmetry zones will not exist.\n");
    exit(1);
  }

  // Step 0c: Set free parameters, overwriting Cparameters defaults
  //          by hand or with command-line input, as desired.
#include "free_parameters.h"

  // Step 0d: Uniform coordinate grids are stored to *xx[3]
  // Step 0d.i: Set bcstruct
  {
    int EigenCoord;
    EigenCoord = 1;
    // Step 0d.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //             chosen Eigen-CoordSystem.
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
    // Step 0e: Find ghostzone mappings; set up bcstruct
    bcstruct_set_up(&griddata.params, griddata.xx, &griddata.bcstruct);
    // Step 0e.i: Free allocated space for xx[][] array
    for(int i=0;i<3;i++) free(griddata.xx[i]);

    // Step 0f: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //          params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //          chosen (non-Eigen) CoordSystem.
    EigenCoord = 0;
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
  }

  // Step 0g: Time coordinate parameters
  const REAL t_final =  0.7*domain_size; /* Final time is set so that at t=t_final,
                                          * data at the origin have not been corrupted
                                          * by the approximate outer boundary condition */

  // Step 0h: Set timestep based on smallest proper distance between gridpoints and CFL factor
  REAL dt = find_timestep(&griddata.params, griddata.xx, CFL_FACTOR);
  //printf("# Timestep set to = %e\n",(double)dt);
  int N_final = (int)(t_final / dt + 0.5); // The number of points in time.
                                           // Add 0.5 to account for C rounding down
                                           // typecasts to integers.
  int output_every_N = (int)((REAL)N_final/800.0);
  if(output_every_N == 0) output_every_N = 1;

  // Step 0i: Error out if the number of auxiliary gridfunctions outnumber evolved gridfunctions.
  //              This is a limitation of the RK method. You are always welcome to declare & allocate
  //              additional gridfunctions by hand.
  if(NUM_AUX_GFS > NUM_EVOL_GFS) {
    printf("Error: NUM_AUX_GFS > NUM_EVOL_GFS. Either reduce the number of auxiliary gridfunctions,\n");
    printf("       or allocate (malloc) by hand storage for *diagnostic_output_gfs. \n");
    exit(1);
  }

  // Step 0j: Declare struct for gridfunctions and allocate memory for y_n_gfs gridfunctions
  MoL_malloc_y_n_gfs(&griddata.params, &griddata.gridfuncs);

  // Step 0.l: Set up initial data to be exact solution at time=0:
  //griddata.params.time = 0.0; //exact_solution_all_points(&griddata.params, griddata.xx, griddata.gridfuncs.y_n_gfs);

  // Step 0.m: Allocate memory for non y_n_gfs. We do this here to free up
  //         memory for setting up initial data (for cases in which initial
  //         data setup is memory intensive.)
  MoL_malloc_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  
  const int Nxx_plus_2NGHOSTS0 = griddata.params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata.params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata.params.Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  for(int which_gf=0; which_gf<NUM_AUXEVOL_GFS; which_gf++) for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
        griddata.gridfuncs.auxevol_gfs[IDX4ptS(which_gf, ijk)] = 0.0/0.0;
      }
        
  for(int which_gf=0; which_gf<NUM_EVOL_GFS; which_gf++) for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
      griddata.gridfuncs.y_n_gfs[IDX4ptS(which_gf, ijk)] = 0.0/0.0;
      griddata.gridfuncs.k_odd_gfs[IDX4ptS(which_gf, ijk)] = 0.0/0.0;
      griddata.gridfuncs.k_even_gfs[IDX4ptS(which_gf, ijk)] = 0.0/0.0;
    }

// Step 1: Set up initial data to be exact solution at time=0:
  REAL time = 0.0;
  calculate_metric_gfs(&griddata.params, griddata.xx, griddata.gridfuncs.auxevol_gfs);
  const char *initial_data_option = argv[4];
  initial_data(initial_data_option, &griddata.params, griddata.xx, griddata.gridfuncs.y_n_gfs, griddata.gridfuncs.auxevol_gfs);
  compute_B_and_Bstagger_from_A(&griddata.params, griddata.gridfuncs.y_n_gfs, griddata.gridfuncs.auxevol_gfs);

  primitive_quantities prims;
  con2prim_diagnostics diagnostics;
  metric_quantities metric;

  for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {

    initialize_primitives(
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(RHOBGF, ijk)], 
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(PGF, ijk)], 
             0.0/0.0, /*griddata.gridfuncs.auxevol_gfs[IDX4ptS(EPSILONGF, ijk)],*/
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(VU0GF, ijk)], 
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(VU1GF, ijk)], 
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(VU2GF, ijk)],
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(BU0GF, ijk)], 
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(BU1GF, ijk)], 
             griddata.gridfuncs.auxevol_gfs[IDX4ptS(BU2GF, ijk)],
             0.0/0.0,
             0.0/0.0,
             0.0/0.0,
             &prims);

    REAL cs2;
    eos.compute_h_and_cs2(&eos, &prims, 
                          &griddata.gridfuncs.auxevol_gfs[IDX4ptS(HGF, ijk)],
                          &cs2);

    if(prims.vx==0 && prims.vy==0 && prims.vz==0) {
        griddata.gridfuncs.auxevol_gfs[IDX4ptS(U4UTGF, ijk)] = 1.0;
    } else {
        limit_v_and_compute_u0(&eos,
                               &metric,
                               &prims,
                               &diagnostics);
        printf("check = %d\n", diagnostics.failure_checker);
   }
  }

  prims_to_cons(&griddata.params, 
                 griddata.xx, 
                 griddata.gridfuncs.auxevol_gfs, 
                 griddata.gridfuncs.y_n_gfs);  
  
  // Step 2c: Output relative error between exact & numerical at center of grid.
  //const int i0mid=Nxx_plus_2NGHOSTS0/2;
  const int i1mid=Nxx_plus_2NGHOSTS1/2;
  const int i2mid=Nxx_plus_2NGHOSTS2/2;
  char filename[100];
  sprintf(filename,"output/out%d__%s.txt", Nxx_plus_2NGHOSTS0, initial_data_option);
  FILE *out2D = fopen(filename, "w");
  for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
  const int idx = IDX3S(i0,i1mid,i2mid);
  fprintf(out2D,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
          griddata.xx[0][i0],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(STILDED0GF,idx)],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(STILDED1GF,idx)],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(STILDED2GF,idx)],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(RHO_STARGF,idx)],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(TAU_TILDEGF,idx)],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(AD0GF,idx)],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(AD1GF,idx)],
          griddata.gridfuncs.y_n_gfs[IDX4ptS(AD2GF,idx)],                   
          griddata.gridfuncs.auxevol_gfs[IDX4ptS(VU0GF,idx)],
          griddata.gridfuncs.auxevol_gfs[IDX4ptS(VU1GF,idx)],
          griddata.gridfuncs.auxevol_gfs[IDX4ptS(VU2GF,idx)],                                  
          griddata.gridfuncs.auxevol_gfs[IDX4ptS(RHOBGF,idx)],
          griddata.gridfuncs.auxevol_gfs[IDX4ptS(PGF,idx)]);
    }
    fclose(out2D);

/*
  for(int n=0;n<=N_final;n++)
    { // Main loop to progress forward in time.

    // Step 1: Set current time to correct value & compute exact solution
    //griddata.params.time = ((REAL)n)*dt;

    // Step 3: Evolve initial data forward in time using Method of Lines with RK4 algorithm,
    //         applying outer boundary conditions.
    // Step 3.b: Step forward one timestep (t -> t+dt) in time using
    //           chosen RK-like MoL timestepping algorithm
    MoL_step_forward_in_time(&griddata, dt);

  } // End main loop to progress forward in time.
  
  */

  // Step 4: Free all allocated memory

  free(griddata.bcstruct.inner_bc_array);
  for(int ng=0;ng<NGHOSTS*3;ng++) free(griddata.bcstruct.pure_outer_bc_array[ng]);
  MoL_free_memory_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  MoL_free_memory_non_y_n_gfs(&griddata.params, &griddata.gridfuncs);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
  return 0;
}
