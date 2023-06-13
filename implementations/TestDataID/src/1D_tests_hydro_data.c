#include "TestDataID.h"

void TestDataID_1D_tests_hydro_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestDataID_1D_tests_hydro_data;
  DECLARE_CCTK_PARAMETERS;

  double rho_l, rho_r;
  double press_l, press_r;
  double vx_l, vy_l, vz_l;
  double vx_r, vy_r, vz_r;
  if(CCTK_EQUALS(test_1D_initial_data,"Balsara1")) {
    rho_l = press_l = 1.0;
    rho_r = 0.125;
    press_r = 0.1;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara2")) {
    press_l = 30.0;
    rho_l = rho_r = press_r = 1.0;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara3")) {
    rho_l = rho_r = 1.0;
    press_l = 1000.0;
    press_r = 0.1;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara4")) {
    rho_l = rho_r = 1.0;
    press_l = press_r = 0.1;
    vx_l = 0.999;
    vx_r = -0.999;
    vy_l = vz_l = 0.0;
    vy_r = vz_r = 0.0;
  } else if(CCTK_EQUALS(test_1D_initial_data,"Balsara5")) {
    rho_l = 1.08;
    press_l = 0.95;
    rho_r = press_r = 1.0;
    vx_l = 0.4;
    vy_l = 0.3;
    vx_r = -0.45;
    vy_r = -0.2;
    vz_l = vz_r = 0.2;
  } else if(CCTK_EQUALS(test_1D_initial_data,"equilibrium")) {
    rho_l = rho_r = press_l = press_r = 1.0;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  /*
  } else if(CCTK_EQUALS(test_1D_initial_data,"sound wave")) {
    this case is handled in the loop because it isn't a
    step function but a sin() wave
  */
  } else if(CCTK_EQUALS(test_1D_initial_data,"shock tube")) {
    rho_l = press_l = 2.0;
    rho_r = press_r = 1.0;
    vx_l = vy_l = vz_l = 0.0;
    vx_r = vy_r = vz_r = 0.0;
  } else {
    CCTK_VERROR("Parameter test_1D_initial_data is not set "
                "to a valid test name. Something has gone wrong.");
  }

  //Above data assumes that the shock is in the x direction; we just have to rotate
  //the data for it to work in other directions.
  if(CCTK_EQUALS(test_shock_direction, "y")) {
    //x-->y, y-->z, z-->x
    double vxtmp, vytmp, vztmp;
    vxtmp = vx_l; vytmp = vy_l; vztmp = vz_l;
    vy_l = vxtmp; vz_l = vytmp; vx_l = vztmp;

    vxtmp = vx_r; vytmp = vy_r; vztmp = vz_r;
    vy_r = vxtmp; vz_r = vytmp; vx_r = vztmp;
  } else if(CCTK_EQUALS(test_shock_direction, "z")) {
    //x-->z, y-->x, z-->y
    double vxtmp, vytmp, vztmp;
    vxtmp = vx_l; vytmp = vy_l; vztmp = vz_l;
    vz_l = vxtmp; vx_l = vytmp; vy_l = vztmp;

    vxtmp = vx_r; vytmp = vy_r; vztmp = vz_r;
    vz_r = vxtmp; vx_r = vytmp; vy_r = vztmp;
  }

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        int ind4x = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        int ind4y = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        int ind4z = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        double step = x[index];
        if(CCTK_EQUALS(test_shock_direction, "y")) {
          step = y[index];
        } else if(CCTK_EQUALS(test_shock_direction, "z")) {
          step = z[index];
        }
    
        if(CCTK_EQUALS(test_1D_initial_data,"sound wave")) {
          rho[index]   = 1.0;
          press[index] = 1.0; // should add kinetic energy here
          vel[ind4x]   = test_wave_amplitude * sin(M_PI * step);
          vel[ind4y]   = 0.0;
          vel[ind4z]   = 0.0;
        } else if(step <= discontinuity_position) {
          rho[index]   = rho_l;
          press[index] = press_l;
          vel[ind4x]   = vx_l;
          vel[ind4y]   = vy_l;
          vel[ind4z]   = vz_l;
        } else {
          rho[index]   = rho_r;
          press[index] = press_r;
          vel[ind4x]   = vx_r;
          vel[ind4y]   = vy_r;
          vel[ind4z]   = vz_r;
        }
        double Gamma = ghl_eos->Gamma_ppoly[
                                 ghl_hybrid_find_polytropic_index(
                                             ghl_eos, rho[index])];
        eps[index] = press[index]/( rho[index]*(Gamma-1) );
      }
    }
  }
  CCTK_VINFO("Finished writing hydrodynamic ID for %s test", test_1D_initial_data);
}
