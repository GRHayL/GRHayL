#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Calculate the metric gridfunctions
 */
void calculate_metric_gfs(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict auxevol_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        const REAL xx0 = xx[0][i0];
        {
          /*
           * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
           */
          const double KDD00 = auxevol_gfs[IDX4S(KDD00GF, i0,i1,i2)];
          const double KDD01 = auxevol_gfs[IDX4S(KDD01GF, i0,i1,i2)];
          const double KDD02 = auxevol_gfs[IDX4S(KDD02GF, i0,i1,i2)];
          const double KDD11 = auxevol_gfs[IDX4S(KDD11GF, i0,i1,i2)];
          const double KDD12 = auxevol_gfs[IDX4S(KDD12GF, i0,i1,i2)];
          const double KDD22 = auxevol_gfs[IDX4S(KDD22GF, i0,i1,i2)];
          /*
           * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
           */
          auxevol_gfs[IDX4S(GAMMADD00GF, i0, i1, i2)] = 1.0;
          auxevol_gfs[IDX4S(GAMMADD01GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(GAMMADD02GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(GAMMADD11GF, i0, i1, i2)] = 1.0;
          auxevol_gfs[IDX4S(GAMMADD12GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(GAMMADD22GF, i0, i1, i2)] = 1.0;
          auxevol_gfs[IDX4S(KDD00GF, i0, i1, i2)] = KDD00;
          auxevol_gfs[IDX4S(KDD01GF, i0, i1, i2)] = KDD01;
          auxevol_gfs[IDX4S(KDD02GF, i0, i1, i2)] = KDD02;
          auxevol_gfs[IDX4S(KDD11GF, i0, i1, i2)] = KDD11;
          auxevol_gfs[IDX4S(KDD12GF, i0, i1, i2)] = KDD12;
          auxevol_gfs[IDX4S(KDD22GF, i0, i1, i2)] = KDD22;
          auxevol_gfs[IDX4S(BETAU0GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(BETAU1GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(BETAU2GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(ALPHAGF, i0, i1, i2)] = 1;
        }
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
