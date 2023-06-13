#include "TestDataID.h"

void TestDataID_1D_tests_magnetic_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestDataID_1D_tests_magnetic_data;
  DECLARE_CCTK_PARAMETERS;

  double Bx_l, By_l, Bz_l;
  double Bx_r, By_r, Bz_r;
  if(CCTK_EQUALS(test_1D_initial_data,"Balsara1")) {
    Bx_l = Bx_r = 0.5;
    By_l = 1.0;
    By_r = -1.0;
    Bz_l = Bz_r = 0.0;   
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara2")) {
    Bx_l = Bx_r = 5.0;
    By_l = Bz_l = 6.0;
    By_r = Bz_r = 0.7;
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara3")) {
    Bx_l = Bx_r = 10.0;
    By_l = Bz_l = 7.0;
    By_r = Bz_r = 0.7;
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara4")) {
    Bx_l = Bx_r = 10.0;
    By_l = Bz_l = 7.0;
    By_r = Bz_r = -7.0; 
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara5")) {
    Bx_l = Bx_r = 2.0;
    By_l = Bz_l = 0.3;
    By_r = -0.7;
    Bz_r = 0.5; 
  } else if(CCTK_EQUALS(test_1D_initial_data,"equilibrium")
         || CCTK_EQUALS(test_1D_initial_data,"sound wave")
         || CCTK_EQUALS(test_1D_initial_data,"shock tube")) {
    Bx_l = By_l = Bz_l = 0.0;
    Bx_r = By_r = Bz_r = 0.0;   
  } else {
    CCTK_VERROR("Parameter test_1D_initial_data is not set "
                "to a valid test name. Something has gone wrong.");
  }

  //Above data assumes that the shock is in the x direction; we just have to rotate
  //the data for it to work in other directions.
  if(CCTK_EQUALS(test_shock_direction, "y")) {
    //x-->y, y-->z, z-->x
    double Bxtmp, Bytmp, Bztmp;
    Bxtmp = Bx_l; Bytmp = By_l; Bztmp = Bz_l;
    By_l = Bxtmp; Bz_l = Bytmp; Bx_l = Bztmp;

    Bxtmp = Bx_r; Bytmp = By_r; Bztmp = Bz_r;
    By_r = Bxtmp; Bz_r = Bytmp; Bx_r = Bztmp;
  } else if(CCTK_EQUALS(test_shock_direction, "z")) {
    //x-->z, y-->x, z-->y
    double Bxtmp, Bytmp, Bztmp;
    Bxtmp = Bx_l; Bytmp = By_l; Bztmp = Bz_l;
    Bz_l = Bxtmp; Bx_l = Bytmp; By_l = Bztmp;

    Bxtmp = Bx_r; Bytmp = By_r; Bztmp = Bz_r;
    Bz_r = Bxtmp; Bx_r = Bytmp; By_r = Bztmp;
  }

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int ind4x = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        const int ind4y = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        const int ind4z = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        if(CCTK_EQUALS(test_shock_direction, "x")) {
          if(x[index] <= discontinuity_position) { 
            Avec[ind4x] = By_l * z[index] - Bz_l * y[index];  
          } else { 
            Avec[ind4x] = By_r * z[index] - Bz_r * y[index];
          }
          Avec[ind4y] = 0.0;
          Avec[ind4z] = Bx_r * y[index];

        } else if(CCTK_EQUALS(test_shock_direction, "y")) {
          if(y[index] <= discontinuity_position) {
            Avec[ind4y] = Bz_l * x[index] - Bx_l * z[index];
          } else {
            Avec[ind4y] = Bz_r * x[index] - Bx_r * z[index];
          }
          Avec[ind4x] = By_r * z[index];
          Avec[ind4z] = 0.0;

        } else if(CCTK_EQUALS(test_shock_direction, "z")) {
          if(z[index] <= discontinuity_position) {
            Avec[ind4z] = Bx_l * y[index] - By_l * x[index];
          } else {
            Avec[ind4z] = Bx_r * y[index] - By_r * x[index];
          }
          Avec[ind4x] = 0.0;
          Avec[ind4y] = Bz_r * x[index];
        }
      }
    }
  }
  CCTK_VINFO("Finished writing magnetic ID for %s test", test_1D_initial_data);
}
