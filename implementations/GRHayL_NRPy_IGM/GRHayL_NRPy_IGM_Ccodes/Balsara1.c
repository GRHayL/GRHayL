#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Generate Balsara1 1D initial data
 */
void Balsara1(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs) {
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
           * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
           */
          const double FDPart3_0 = (1.0/(-TINYDOUBLE + xx0));
          const double FDPart3_2 = (1.0/2.0)*xx0 - 1.0/2.0*fabs(xx0);
          const double FDPart3_4 = ((1.0/2.0)*TINYDOUBLE + (1.0/2.0)*xx0 + (1.0/2.0)*fabs(TINYDOUBLE + xx0))/(TINYDOUBLE + xx0);
          evol_gfs[IDX4S(AD0GF, i0, i1, i2)] = (0.5*dxx2 + xx2)*(FDPart3_0*FDPart3_2 - FDPart3_4);
          evol_gfs[IDX4S(AD1GF, i0, i1, i2)] = 0;
          evol_gfs[IDX4S(AD2GF, i0, i1, i2)] = 0.25*dxx1 + (1.0/2.0)*xx1;
          auxevol_gfs[IDX4S(VU0GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(VU1GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(VU2GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(RHOBGF, i0, i1, i2)] = FDPart3_0*FDPart3_2 + (1.0/8.0)*FDPart3_4;
          auxevol_gfs[IDX4S(PGF, i0, i1, i2)] = FDPart3_0*FDPart3_2 + (1.0/10.0)*FDPart3_4;
          evol_gfs[IDX4S(PSI6PHIGF, i0, i1, i2)] = 0;
        }
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
