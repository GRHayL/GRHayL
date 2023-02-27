#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"

// Side note: the following values could be used for cell averaged gfs:
//     am2=-1.0/12.0, am1=7.0/12.0, a0=7.0/12.0, a1=-1.0/12.0
// However, since the metric gfs store the grid point values instead of the cell average,
//     the following coefficients should be used:
//     am2 = -1/16, am1 = 9/16, a0 = 9/16, a1 = -1/16
// This will yield the third-order-accurate face values at m-1/2,
//      using values specified at {m-2,m-1,m,m+1}
#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))


const int metric_gfs_list[10] = {GAMMADD00GF,
                                 GAMMADD01GF,
                                 GAMMADD02GF,
                                 GAMMADD11GF,
                                 GAMMADD12GF,
                                 GAMMADD22GF,
                                 BETAU0GF,
                                 BETAU1GF,
                                 BETAU2GF,
                                 ALPHAGF};

const int metric_gfs_face_list[10] = {GAMMA_FACEDD00GF,
                                      GAMMA_FACEDD01GF,
                                      GAMMA_FACEDD02GF,
                                      GAMMA_FACEDD11GF,
                                      GAMMA_FACEDD12GF,
                                      GAMMA_FACEDD22GF,
                                      BETA_FACEU0GF,
                                      BETA_FACEU1GF,
                                      BETA_FACEU2GF,
                                      ALPHA_FACEGF};

const int num_metric_gfs = 10;


/*
 * Interpolate metric gridfunctions to cell faces
 */
void interpolate_metric_gfs_to_cell_faces(const paramstruct *params,REAL *auxevol_gfs,const int flux_dirn) {
#include "./set_Cparameters.h"
    int in_gf,out_gf;
    REAL Qm2,Qm1,Qp0,Qp1;


    
    const int dirn0 = (flux_dirn == 0);
    const int dirn1 = (flux_dirn == 1);
    const int dirn2 = (flux_dirn == 2);
    
    for(int gf = 0;gf < num_metric_gfs;gf++) {
        in_gf  = metric_gfs_list[gf];
        out_gf = metric_gfs_face_list[gf];
        for (int i2 = 2;i2 < Nxx_plus_2NGHOSTS2-1;i2++) {
            for (int i1 = 2;i1 < Nxx_plus_2NGHOSTS1-1;i1++) {
                for (int i0 = 2;i0 < Nxx_plus_2NGHOSTS0-1;i0++) {
                    Qm2 = auxevol_gfs[IDX4S(in_gf,i0-2*dirn0,i1-2*dirn1,i2-2*dirn2)];
                    Qm1 = auxevol_gfs[IDX4S(in_gf,i0-  dirn0,i1-  dirn1,i2-  dirn2)];
                    Qp0 = auxevol_gfs[IDX4S(in_gf,i0,        i1,        i2        )];
                    Qp1 = auxevol_gfs[IDX4S(in_gf,i0+  dirn0,i1+  dirn1,i2+  dirn2)];
                    auxevol_gfs[IDX4S(out_gf,i0,i1,i2)] = COMPUTE_FCVAL(Qm2,Qm1,Qp0,Qp1);

                }
            }
        }
    }
}
