#include "TestDataIDX.h"

extern "C" void TestDataIDX_1D_tests_magnetic_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestDataIDX_1D_tests_magnetic_data;
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

//  grid.loop_all_device<1, 0, 0>(
//      grid.nghostzones,
//      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
//
//    if(CCTK_EQUALS(test_shock_direction, "x")) {
//      if(p.x <= discontinuity_position) { 
//        Avec_x(p.I) = By_l * (p.z) - Bz_l * (p.y);  
//      } else { 
//        Avec_x(p.I) = By_r * (p.z) - Bz_r * (p.y);
//      }
//    } else if(CCTK_EQUALS(test_shock_direction, "y")) {
//      //Avec_x(p.I) = By_r * (p.z);
//    } else if(CCTK_EQUALS(test_shock_direction, "z")) {
//      //Avec_x(p.I) = 0.0;
//    }
//  });

//  grid.loop_all_device<0, 1, 0>(
//      grid.nghostzones,
//      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
//
//    if(CCTK_EQUALS(test_shock_direction, "x")) {
//      //Avec_y(p.I) = 0.0;
//    } else if(CCTK_EQUALS(test_shock_direction, "y")) {
//      //if(p.y <= discontinuity_position) {
//      //  Avec_y(p.I) = Bz_l * (p.x) - Bx_l * (p.z);
//      //} else {
//      //  Avec_y(p.I) = Bz_r * (p.x) - Bx_r * (p.z);
//      //}
//    } else if(CCTK_EQUALS(test_shock_direction, "z")) {
//      //Avec_y(p.I) = Bz_r * (p.x);
//    }
//  });
//
//  grid.loop_all_device<0, 0, 1>(
//      grid.nghostzones,
//      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
//
//    if(CCTK_EQUALS(test_shock_direction, "x")) {
//      //Avec_z(p.I) = Bx_r * (p.y);
//    } else if(CCTK_EQUALS(test_shock_direction, "y")) {
//      //Avec_z(p.I) = 0.0;
//    } else if(CCTK_EQUALS(test_shock_direction, "z")) {
//      //if(p.z <= discontinuity_position) {
//      //  Avec_z(p.I) = Bx_l * (p.y) - By_l * (p.x);
//      //} else {
//      //  Avec_z(p.I) = Bx_r * (p.y) - By_r * (p.x);
//      //}
//    }
//  });
  CCTK_VINFO("Finished writing magnetic ID for %s test", test_1D_initial_data);
}
