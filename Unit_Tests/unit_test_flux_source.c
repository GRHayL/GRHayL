#include "unit_tests.h"
#include "flux_source_unit_test.h"

static inline void compute_h_and_cs2(struct eos_parameters const *restrict eos,
                                     primitive_quantities const *restrict prims,
                                     double *restrict h,
                                     double *restrict cs2) {


// CCTK_REAL h_ = U[PRESSURE]*U[VX]/U[VZ];
*h = prims->press*prims->vx / prims->vz;

// CCTK_REAL c_s_squared  = U[PRESSURE]*U[VZ]/(h);
*cs2 = prims->rho*prims->vz/(*h);
// printf("works!!\n");
}

#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

static inline void calculate_metric_face_values(const int flux_dirn, 
                                  const int metric_faces_gfs[],
                                  const int metric_gfs[],
                                  double *restrict auxevol_gfs) {

  LOOP_REGION(2, Nxx_plus_2NGHOSTS0 - 1,
              2, Nxx_plus_2NGHOSTS1 - 1,
              2, Nxx_plus_2NGHOSTS2 - 1) {

    int idxm2 = IDX3S(i0 - 2*kronecker_delta[flux_dirn+1][0], 
                      i1 - 2*kronecker_delta[flux_dirn+1][1], 
                      i2 - 2*kronecker_delta[flux_dirn+1][2]);
    
    int idxm1 = IDX3S(i0 - 1*kronecker_delta[flux_dirn+1][0], 
                      i1 - 1*kronecker_delta[flux_dirn+1][1], 
                      i2 - 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    for(int which_gf=0; which_gf<10; which_gf++){
      
      int gf_face = metric_faces_gfs[which_gf];
      int gf = metric_gfs[which_gf];
      auxevol_gfs[IDX4ptS(gf_face, idx)] = COMPUTE_FCVAL(auxevol_gfs[IDX4ptS(gf, idxm2)],
                                                         auxevol_gfs[IDX4ptS(gf, idxm1)],
                                                         auxevol_gfs[IDX4ptS(gf, idx  )],
                                                         auxevol_gfs[IDX4ptS(gf, idxp1)]);
    }
  }
}

/*
 * // main() function:
 * // Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
 * // Step 1: Write test data to gridfunctions
 * // Step 2: Overwrite all data in ghost zones with NaNs
 * // Step 3: Apply curvilinear boundary conditions
 * // Step 4: Print gridfunction data after curvilinear boundary conditions have been applied
 * // Step 5: Free all allocated memory
 */
int main(int argc, char **argv) {

  eos_parameters eos;
  eos.compute_h_and_cs2 = &compute_h_and_cs2;  
  
  // Step 0.m: Allocate memory for non y_n_gfs. We do this here to free up
  //         memory for setting up initial data (for cases in which initial
  //         data setup is memory intensive.)
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  double *restrict     evol_gfs    = (double *restrict)malloc(sizeof(double) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  double *restrict etk_evol_gfs    = (double *restrict)malloc(sizeof(double) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  double *restrict auxevol_gfs = (double *restrict)malloc(sizeof(double) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);

  // Now we need to properly populate our structs for the unit test.

  primitive_quantities prims, prims_r, prims_l;

  conservative_quantities cons_fluxes, cons_sources;
  
  metric_quantities metric, metric_face;
  
  extrinsic_curvature curv;
  
  metric_derivatives metric_derivs;

  // We begin by generating random inital data.

  for(int which_gf=0; which_gf<NUM_AUXEVOL_GFS; which_gf++) for(int idx=0; idx<Nxx_plus_2NGHOSTS_tot; idx++) {
    auxevol_gfs[IDX4ptS(which_gf, idx)] = 0.0 / 0.0; // ((double)rand() - (double)rand())/(double)(RAND_MAX);
        }

  for(int which_gf=0; which_gf<NUM_EVOL_GFS; which_gf++) for(int idx=0; idx<Nxx_plus_2NGHOSTS_tot; idx++) {
    evol_gfs[IDX4ptS(which_gf, idx)] = 0.0/0.0;
        }

  /*
  After reading in the same random data used in the original IGM code, now we
  calculate face cell-face values and derivatives at cell centers
  */

  int metric_faces_gfs[10] = {ALPHA_FACEGF, BETA_FACEU0GF, BETA_FACEU1GF, BETA_FACEU2GF,
                          GAMMA_FACEDD00GF, GAMMA_FACEDD01GF, GAMMA_FACEDD02GF,
                          GAMMA_FACEDD11GF, GAMMA_FACEDD12GF, GAMMA_FACEDD22GF};

  int metric_gfs[10] = {ALPHAGF, BETAU0GF, BETAU1GF, BETAU2GF,
                          GAMMADD00GF, GAMMADD01GF, GAMMADD02GF,
                          GAMMADD11GF, GAMMADD12GF, GAMMADD22GF};

  int flux_dirn = 0;
  read_from_binary_file_all("flux_source_dirn1.bin", auxevol_gfs);

  calculate_metric_face_values( flux_dirn, 
                                metric_faces_gfs,
                                metric_gfs,
                                auxevol_gfs);

  LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - (NGHOSTS),
              NGHOSTS, Nxx_plus_2NGHOSTS1 - (NGHOSTS),
              NGHOSTS, Nxx_plus_2NGHOSTS2 - (NGHOSTS)) {
    
    int idx  = IDX3S(i0, i1, i2);

    // CCTK_REAL u0_l = Ul[RHOB]*Ul[BX_CENTER]/Ul[VY];

    // CCTK_REAL h_r = Ur[PRESSURE]*Ur[VX]/Ur[VZ];

    auxevol_gfs[IDX4ptS(U4U0GF, idx)] = auxevol_gfs[IDX4ptS(RHOBGF, idx)]*auxevol_gfs[IDX4ptS(BU0GF, idx)] / auxevol_gfs[IDX4ptS(VU1GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU0GF, idx)] = auxevol_gfs[IDX4ptS(RHOB_RGF, idx)]*auxevol_gfs[IDX4ptS(BRU0GF, idx)] / auxevol_gfs[IDX4ptS(VRU1GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU0GF, idx)] = auxevol_gfs[IDX4ptS(RHOB_LGF, idx)]*auxevol_gfs[IDX4ptS(BLU0GF, idx)] / auxevol_gfs[IDX4ptS(VLU1GF, idx)];

    auxevol_gfs[IDX4ptS(HGF, idx)] = auxevol_gfs[IDX4ptS(PGF, idx)]*auxevol_gfs[IDX4ptS(VU0GF, idx)] / auxevol_gfs[IDX4ptS(VU2GF, idx)];
    auxevol_gfs[IDX4ptS(H_RGF, idx)] = auxevol_gfs[IDX4ptS(P_RGF, idx)]*auxevol_gfs[IDX4ptS(VRU0GF, idx)] / auxevol_gfs[IDX4ptS(VRU2GF, idx)];
    auxevol_gfs[IDX4ptS(H_LGF, idx)] = auxevol_gfs[IDX4ptS(P_LGF, idx)]*auxevol_gfs[IDX4ptS(VLU0GF, idx)] / auxevol_gfs[IDX4ptS(VLU2GF, idx)];

    auxevol_gfs[IDX4ptS(U4U1GF, idx)] = auxevol_gfs[IDX4ptS(VU0GF, idx)]*auxevol_gfs[IDX4ptS(U4U0GF, idx)];
    auxevol_gfs[IDX4ptS(U4U2GF, idx)] = auxevol_gfs[IDX4ptS(VU1GF, idx)]*auxevol_gfs[IDX4ptS(U4U0GF, idx)];
    auxevol_gfs[IDX4ptS(U4U3GF, idx)] = auxevol_gfs[IDX4ptS(VU2GF, idx)]*auxevol_gfs[IDX4ptS(U4U0GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU1GF, idx)] = auxevol_gfs[IDX4ptS(VRU0GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU2GF, idx)] = auxevol_gfs[IDX4ptS(VRU1GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU3GF, idx)] = auxevol_gfs[IDX4ptS(VRU2GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU1GF, idx)] = auxevol_gfs[IDX4ptS(VLU0GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU2GF, idx)] = auxevol_gfs[IDX4ptS(VLU1GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU3GF, idx)] = auxevol_gfs[IDX4ptS(VLU2GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];

    initialize_structs(idx, auxevol_gfs, 
                            &prims,
                            &prims_r,
                            &prims_l,
                            &metric,
                            &metric_face,
                            &curv,
                            &metric_derivs);

    calculate_HLLE_fluxes_dirn0(&prims_r, 
                           &prims_l,
                           &eos,
                           &metric_face, 
                           &cons_fluxes);

    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)]  = cons_fluxes.S_x;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)]  = cons_fluxes.S_y;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)]  = cons_fluxes.S_z;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)]  = cons_fluxes.rho;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)]  = cons_fluxes.tau;
  }


  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

    evol_gfs[IDX4ptS(STILDED0GF, idx)]  = invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED1GF, idx)]  = invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED2GF, idx)]  = invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idxp1)]);
    evol_gfs[IDX4ptS(RHO_STARGF, idx)]  = invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idxp1)]);
    evol_gfs[IDX4ptS(TAU_TILDEGF,idx)]  = invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idxp1)]);


    auxevol_gfs[IDX4ptS(ALPHA_DD0GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idxp1)]
                                                    - auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD00GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD10GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD20GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx  )]);


    auxevol_gfs[IDX4ptS(GAMMADD_DD000GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD010GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD020GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD110GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD120GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD220GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx  )]);
  }


  flux_dirn = 1;
  read_from_binary_file_recons("flux_source_dirn2.bin", auxevol_gfs);

  calculate_metric_face_values( flux_dirn, 
                                metric_faces_gfs,
                                metric_gfs,
                                auxevol_gfs);

  
  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int idx  = IDX3S(i0, i1, i2);

    auxevol_gfs[IDX4ptS(U4U0GF, idx)] = auxevol_gfs[IDX4ptS(RHOBGF, idx)]*auxevol_gfs[IDX4ptS(BU0GF, idx)] / auxevol_gfs[IDX4ptS(VU1GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU0GF, idx)] = auxevol_gfs[IDX4ptS(RHOB_RGF, idx)]*auxevol_gfs[IDX4ptS(BRU0GF, idx)] / auxevol_gfs[IDX4ptS(VRU1GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU0GF, idx)] = auxevol_gfs[IDX4ptS(RHOB_LGF, idx)]*auxevol_gfs[IDX4ptS(BLU0GF, idx)] / auxevol_gfs[IDX4ptS(VLU1GF, idx)];

    auxevol_gfs[IDX4ptS(HGF, idx)] = auxevol_gfs[IDX4ptS(PGF, idx)]*auxevol_gfs[IDX4ptS(VU0GF, idx)] / auxevol_gfs[IDX4ptS(VU2GF, idx)];
    auxevol_gfs[IDX4ptS(H_RGF, idx)] = auxevol_gfs[IDX4ptS(P_RGF, idx)]*auxevol_gfs[IDX4ptS(VRU0GF, idx)] / auxevol_gfs[IDX4ptS(VRU2GF, idx)];
    auxevol_gfs[IDX4ptS(H_LGF, idx)] = auxevol_gfs[IDX4ptS(P_LGF, idx)]*auxevol_gfs[IDX4ptS(VLU0GF, idx)] / auxevol_gfs[IDX4ptS(VLU2GF, idx)];


    auxevol_gfs[IDX4ptS(U4RU1GF, idx)] = auxevol_gfs[IDX4ptS(VRU0GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU2GF, idx)] = auxevol_gfs[IDX4ptS(VRU1GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU3GF, idx)] = auxevol_gfs[IDX4ptS(VRU2GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU1GF, idx)] = auxevol_gfs[IDX4ptS(VLU0GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU2GF, idx)] = auxevol_gfs[IDX4ptS(VLU1GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU3GF, idx)] = auxevol_gfs[IDX4ptS(VLU2GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];

    initialize_structs(idx, auxevol_gfs, 
                            &prims,
                            &prims_r,
                            &prims_l,
                            &metric,
                            &metric_face,
                            &curv,
                            &metric_derivs);

    calculate_HLLE_fluxes_dirn1(&prims_r, 
                           &prims_l,
                           &eos,
                           &metric_face, 
                           &cons_fluxes);

    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)]  = cons_fluxes.S_x;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)]  = cons_fluxes.S_y;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)]  = cons_fluxes.S_z;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)]  = cons_fluxes.rho;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)]  = cons_fluxes.tau;
  }


  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

    evol_gfs[IDX4ptS(STILDED0GF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED1GF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED2GF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idxp1)]);
    evol_gfs[IDX4ptS(RHO_STARGF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idxp1)]);
    evol_gfs[IDX4ptS(TAU_TILDEGF,idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idxp1)]);


    auxevol_gfs[IDX4ptS(ALPHA_DD1GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idxp1)]
                                                    - auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD01GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD11GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD21GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx  )]);


    auxevol_gfs[IDX4ptS(GAMMADD_DD001GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD011GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD021GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD111GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD121GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD221GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx  )]);
  }

  flux_dirn = 2;
  read_from_binary_file_recons("flux_source_dirn3.bin", auxevol_gfs);

  calculate_metric_face_values( flux_dirn, 
                                metric_faces_gfs,
                                metric_gfs,
                                auxevol_gfs);

  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {

    int idx  = IDX3S(i0, i1, i2);

    auxevol_gfs[IDX4ptS(U4U0GF, idx)] = auxevol_gfs[IDX4ptS(RHOBGF, idx)]*auxevol_gfs[IDX4ptS(BU0GF, idx)] / auxevol_gfs[IDX4ptS(VU1GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU0GF, idx)] = auxevol_gfs[IDX4ptS(RHOB_RGF, idx)]*auxevol_gfs[IDX4ptS(BRU0GF, idx)] / auxevol_gfs[IDX4ptS(VRU1GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU0GF, idx)] = auxevol_gfs[IDX4ptS(RHOB_LGF, idx)]*auxevol_gfs[IDX4ptS(BLU0GF, idx)] / auxevol_gfs[IDX4ptS(VLU1GF, idx)];

    auxevol_gfs[IDX4ptS(HGF, idx)] = auxevol_gfs[IDX4ptS(PGF, idx)]*auxevol_gfs[IDX4ptS(VU0GF, idx)] / auxevol_gfs[IDX4ptS(VU2GF, idx)];
    auxevol_gfs[IDX4ptS(H_RGF, idx)] = auxevol_gfs[IDX4ptS(P_RGF, idx)]*auxevol_gfs[IDX4ptS(VRU0GF, idx)] / auxevol_gfs[IDX4ptS(VRU2GF, idx)];
    auxevol_gfs[IDX4ptS(H_LGF, idx)] = auxevol_gfs[IDX4ptS(P_LGF, idx)]*auxevol_gfs[IDX4ptS(VLU0GF, idx)] / auxevol_gfs[IDX4ptS(VLU2GF, idx)];


    auxevol_gfs[IDX4ptS(U4RU1GF, idx)] = auxevol_gfs[IDX4ptS(VRU0GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU2GF, idx)] = auxevol_gfs[IDX4ptS(VRU1GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4RU3GF, idx)] = auxevol_gfs[IDX4ptS(VRU2GF, idx)]*auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU1GF, idx)] = auxevol_gfs[IDX4ptS(VLU0GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU2GF, idx)] = auxevol_gfs[IDX4ptS(VLU1GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];
    auxevol_gfs[IDX4ptS(U4LU3GF, idx)] = auxevol_gfs[IDX4ptS(VLU2GF, idx)]*auxevol_gfs[IDX4ptS(U4LU0GF, idx)];

    initialize_structs(idx, auxevol_gfs, 
                            &prims,
                            &prims_r,
                            &prims_l,
                            &metric,
                            &metric_face,
                            &curv,
                            &metric_derivs);

    calculate_HLLE_fluxes_dirn2(&prims_r, 
                           &prims_l,
                           &eos,
                           &metric_face, 
                           &cons_fluxes);

    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)]  = cons_fluxes.S_x;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)]  = cons_fluxes.S_y;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)]  = cons_fluxes.S_z;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)]  = cons_fluxes.rho;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)]  = cons_fluxes.tau;

    }


  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

    evol_gfs[IDX4ptS(STILDED0GF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED1GF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED2GF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idxp1)]);
    evol_gfs[IDX4ptS(RHO_STARGF, idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idxp1)]);
    evol_gfs[IDX4ptS(TAU_TILDEGF,idx)]  += invdx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idxp1)]);


    auxevol_gfs[IDX4ptS(ALPHA_DD2GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idxp1)]
                                                    - auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD02GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD12GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx  )]);

    auxevol_gfs[IDX4ptS(BETAU_DD22GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idxp1)]
                                                     - auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx  )]);


    auxevol_gfs[IDX4ptS(GAMMADD_DD002GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD012GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD022GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD112GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD122GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx  )]);

    auxevol_gfs[IDX4ptS(GAMMADD_DD222GF, idx)] = invdx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idxp1)]
                                                        - auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx  )]);
  }

    LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
                NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
                NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
        int idx  = IDX3S(i0, i1, i2);
    
      initialize_structs(idx, auxevol_gfs, 
                           &prims,
                           &prims_r,
                           &prims_l,
                           &metric,
                           &metric_face,
                           &curv,
                           &metric_derivs);
                           
           
      calculate_all_source_terms(&prims,
                                 &eos,
                                 &metric,
                                 &curv,
                                 &metric_derivs,
                                 &cons_sources);
                                  
      evol_gfs[IDX4ptS(STILDED0GF, idx)]  += cons_sources.S_x;
      evol_gfs[IDX4ptS(STILDED1GF, idx)]  += cons_sources.S_y;
      evol_gfs[IDX4ptS(STILDED2GF, idx)]  += cons_sources.S_z;
      evol_gfs[IDX4ptS(TAU_TILDEGF, idx)] += cons_sources.tau;

      // evol_gfs[IDX4ptS(STILDED0GF, idx)]  = cons_sources.S_x;
      // evol_gfs[IDX4ptS(STILDED1GF, idx)]  = cons_sources.S_y;
      // evol_gfs[IDX4ptS(STILDED2GF, idx)]  = cons_sources.S_z;
      // evol_gfs[IDX4ptS(TAU_TILDEGF, idx)] = cons_sources.tau;
      // evol_gfs[IDX4ptS(RHO_STARGF, idx)]  = 0;
    }

    FILE *infile = fopen("flux_source_output_rhs_data.bin", "rb");
    double rhs_correct_magic_number = 9.524300707856655e-3;
    double rhs_magic_number1, rhs_magic_number2, rhs_magic_number3;
    
    int key;
    
    key  = fread(&rhs_magic_number1, sizeof(double), 1, infile);
    key += fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*RHO_STARGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    key += fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*TAU_TILDEGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    key += fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*STILDED0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    key += fread(&rhs_magic_number2, sizeof(double), 1, infile);
    key += fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*STILDED1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    key += fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*STILDED2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    key += fread(&rhs_magic_number3, sizeof(double), 1, infile);
    
    fclose(infile);
    if(rhs_magic_number1!=rhs_correct_magic_number){ printf("ERROR: magic_number1 does not match"); exit(1);}
    if(rhs_magic_number2!=rhs_correct_magic_number){ printf("ERROR: magic_number2 does not match"); exit(1);}
    if(rhs_magic_number3!=rhs_correct_magic_number){ printf("ERROR: magic_number3 does not match"); exit(1);}

    LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS-1,
                NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS-1,
                NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS-1) {

      int idx  = IDX3S(i0, i1, i2);

      double S_x_rel_diff, S_y_rel_diff, S_z_rel_diff;
      double rho_rel_diff, tau_rel_diff;

      S_x_rel_diff = log10(fabs(0.5*(evol_gfs[IDX4ptS(STILDED0GF, idx)] - etk_evol_gfs[IDX4ptS(STILDED0GF, idx)]) / (evol_gfs[IDX4ptS(STILDED0GF, idx)] + etk_evol_gfs[IDX4ptS(STILDED0GF, idx)]))),
      S_y_rel_diff = log10(fabs(0.5*(evol_gfs[IDX4ptS(STILDED1GF, idx)] - etk_evol_gfs[IDX4ptS(STILDED1GF, idx)]) / (evol_gfs[IDX4ptS(STILDED1GF, idx)] + etk_evol_gfs[IDX4ptS(STILDED1GF, idx)]))),
      S_z_rel_diff = log10(fabs(0.5*(evol_gfs[IDX4ptS(STILDED2GF, idx)] - etk_evol_gfs[IDX4ptS(STILDED2GF, idx)]) / (evol_gfs[IDX4ptS(STILDED2GF, idx)] + etk_evol_gfs[IDX4ptS(STILDED2GF, idx)]))),
      rho_rel_diff = log10(fabs(0.5*(evol_gfs[IDX4ptS(TAU_TILDEGF, idx)] - etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)]) / (evol_gfs[IDX4ptS(TAU_TILDEGF, idx)] + etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)]))),
      tau_rel_diff = log10(fabs(0.5*(evol_gfs[IDX4ptS(RHO_STARGF, idx)] - etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)]) / (evol_gfs[IDX4ptS(RHO_STARGF, idx)] + etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)])));

      if(S_x_rel_diff > -9){printf("ERROR: S_x_rel_diff is too high"); exit(1);}
      if(S_y_rel_diff > -9){printf("ERROR: S_y_rel_diff is too high"); exit(1);}
      if(S_z_rel_diff > -9){printf("ERROR: S_z_rel_diff is too high"); exit(1);}
      if(rho_rel_diff > -9){printf("ERROR: rho_rel_diff is too high"); exit(1);}
      if(tau_rel_diff > -9){printf("ERROR: tau_rel_diff is too high"); exit(1);}

      // printf("%d %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", 
      // i0, i1, i2,
      // evol_gfs[IDX4ptS(STILDED0GF, idx)], etk_evol_gfs[IDX4ptS(STILDED0GF, idx)],
      // evol_gfs[IDX4ptS(STILDED1GF, idx)], etk_evol_gfs[IDX4ptS(STILDED1GF, idx)],
      // evol_gfs[IDX4ptS(STILDED2GF, idx)], etk_evol_gfs[IDX4ptS(STILDED2GF, idx)],
      // evol_gfs[IDX4ptS(TAU_TILDEGF, idx)], etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)],
      // evol_gfs[IDX4ptS(RHO_STARGF, idx)], etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)]);

  }

 
  // Step 4: Free all allocated memory
  free(evol_gfs);
  free(etk_evol_gfs);
  free(auxevol_gfs);
  return 0;
}
