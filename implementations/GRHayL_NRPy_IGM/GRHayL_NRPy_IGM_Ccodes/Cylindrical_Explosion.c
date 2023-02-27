#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Generate 2D Cylindrical Explosion initial data
 */
void Cylindrical_Explosion(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs) {
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
          const double FDPart3_0 = sqrt(((xx0)*(xx0)) + ((xx1)*(xx1)) + ((xx2)*(xx2)));
          const double FDPart3_1 = FDPart3_0 - 1.0;
          const double FDPart3_5 = ((1.0/2.0)*FDPart3_0 + (1.0/2.0)*TINYDOUBLE + (1.0/2.0)*fabs(FDPart3_1 + TINYDOUBLE) - 0.5)/(FDPart3_1 + TINYDOUBLE);
          const double FDPart3_6 = -FDPart3_0 + TINYDOUBLE;
          const double FDPart3_7 = -FDPart3_6 - 0.80000000000000004;
          const double FDPart3_9 = ((1.0/2.0)*FDPart3_0 - 1.0/2.0*TINYDOUBLE - 1.0/2.0*fabs(FDPart3_7) - 0.40000000000000002)/FDPart3_7;
          const double FDPart3_10 = FDPart3_0 - 0.80000000000000004;
          const double FDPart3_11 = 5.0000000000000009*FDPart3_10;
          const double FDPart3_12 = ((1.0/2.0)*FDPart3_0 - 1.0/2.0*fabs(FDPart3_1) - 0.5)*((1.0/2.0)*FDPart3_0 + (1.0/2.0)*fabs(FDPart3_10) - 0.40000000000000002)/((FDPart3_10 + TINYDOUBLE)*(-FDPart3_6 - 1.0));
          evol_gfs[IDX4S(AD0GF, i0, i1, i2)] = 0;
          evol_gfs[IDX4S(AD1GF, i0, i1, i2)] = 0;
          evol_gfs[IDX4S(AD2GF, i0, i1, i2)] = 0.050000000000000003*dxx1 + (1.0/10.0)*xx1;
          auxevol_gfs[IDX4S(VU0GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(VU1GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(VU2GF, i0, i1, i2)] = 0;
          auxevol_gfs[IDX4S(RHOBGF, i0, i1, i2)] = FDPart3_12*exp(5.0000000000000009*FDPart3_1*log(100) - FDPart3_11*log(10000)) + (1.0/10000.0)*FDPart3_5 + (1.0/100.0)*FDPart3_9;
          auxevol_gfs[IDX4S(PGF, i0, i1, i2)] = FDPart3_12*exp(FDPart3_11*log(3.0/100000.0)) + (3.0/100000.0)*FDPart3_5 + FDPart3_9;
          evol_gfs[IDX4S(PSI6PHIGF, i0, i1, i2)] = 0;
        }
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
