#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Method of Lines (MoL) for "RK4" method: Step forward one full timestep.
 */
void MoL_step_forward_in_time(griddata_struct *restrict griddata, const REAL dt) {

  // C code implementation of -={ RK4 }=- Method of Lines timestepping.

  // -={ START k1 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    REAL * xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata->xx[ww];
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    //rhs_eval(params, xx, y_n_gfs, k_odd_gfs);
    if(params->outer_bc_type == RADIATION_OUTER_BCS)
      apply_bcs_outerradiation_and_inner(params, bcstruct, griddata->xx,
                                         gridfunctions_wavespeed,gridfunctions_f_infinity,
                                         y_n_gfs, k_odd_gfs);
    LOOP_ALL_GFS_GPS(i) {
      const REAL k_odd_gfsL = k_odd_gfs[i];
      const REAL y_n_gfsL = y_n_gfs[i];
      y_nplus1_running_total_gfs[i] = (1.0/6.0)*dt*k_odd_gfsL;
      k_odd_gfs[i] = (1.0/2.0)*dt*k_odd_gfsL + y_n_gfsL;
    }
    if(params->outer_bc_type == EXTRAPOLATION_OUTER_BCS)
      apply_bcs_outerextrap_and_inner(params, bcstruct, k_odd_gfs);
  }
  // -={ END k1 substep }=-

  // -={ START k2 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    REAL * xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata->xx[ww];
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    //rhs_eval(params, xx, k_odd_gfs, k_even_gfs);
    if(params->outer_bc_type == RADIATION_OUTER_BCS)
      apply_bcs_outerradiation_and_inner(params, bcstruct, griddata->xx,
                                         gridfunctions_wavespeed,gridfunctions_f_infinity,
                                         k_odd_gfs, k_even_gfs);
    LOOP_ALL_GFS_GPS(i) {
      const REAL k_even_gfsL = k_even_gfs[i];
      const REAL y_nplus1_running_total_gfsL = y_nplus1_running_total_gfs[i];
      const REAL y_n_gfsL = y_n_gfs[i];
      y_nplus1_running_total_gfs[i] = (1.0/3.0)*dt*k_even_gfsL + y_nplus1_running_total_gfsL;
      k_even_gfs[i] = (1.0/2.0)*dt*k_even_gfsL + y_n_gfsL;
    }
    if(params->outer_bc_type == EXTRAPOLATION_OUTER_BCS)
      apply_bcs_outerextrap_and_inner(params, bcstruct, k_even_gfs);
  }
  // -={ END k2 substep }=-

  // -={ START k3 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    REAL * xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata->xx[ww];
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    //rhs_eval(params, xx, k_even_gfs, k_odd_gfs);
    if(params->outer_bc_type == RADIATION_OUTER_BCS)
      apply_bcs_outerradiation_and_inner(params, bcstruct, griddata->xx,
                                         gridfunctions_wavespeed,gridfunctions_f_infinity,
                                         k_even_gfs, k_odd_gfs);
    LOOP_ALL_GFS_GPS(i) {
      const REAL k_odd_gfsL = k_odd_gfs[i];
      const REAL y_nplus1_running_total_gfsL = y_nplus1_running_total_gfs[i];
      const REAL y_n_gfsL = y_n_gfs[i];
      y_nplus1_running_total_gfs[i] = (1.0/3.0)*dt*k_odd_gfsL + y_nplus1_running_total_gfsL;
      k_odd_gfs[i] = dt*k_odd_gfsL + y_n_gfsL;
    }
    if(params->outer_bc_type == EXTRAPOLATION_OUTER_BCS)
      apply_bcs_outerextrap_and_inner(params, bcstruct, k_odd_gfs);
  }
  // -={ END k3 substep }=-

  // -={ START k4 substep }=-
  {
    // Set gridfunction aliases from gridfuncs struct
    REAL *restrict y_n_gfs = griddata->gridfuncs.y_n_gfs;  // y_n gridfunctions
    // Temporary timelevel & AUXEVOL gridfunctions:
    REAL *restrict y_nplus1_running_total_gfs = griddata->gridfuncs.y_nplus1_running_total_gfs;
    REAL *restrict k_odd_gfs = griddata->gridfuncs.k_odd_gfs;
    REAL *restrict k_even_gfs = griddata->gridfuncs.k_even_gfs;
    REAL *restrict auxevol_gfs = griddata->gridfuncs.auxevol_gfs;
    paramstruct *restrict params = &griddata->params;
    REAL * xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata->xx[ww];
    const bc_struct *restrict bcstruct = &griddata->bcstruct;
    const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;
    //rhs_eval(params, xx, k_odd_gfs, k_even_gfs);
    if(params->outer_bc_type == RADIATION_OUTER_BCS)
      apply_bcs_outerradiation_and_inner(params, bcstruct, griddata->xx,
                                         gridfunctions_wavespeed,gridfunctions_f_infinity,
                                         k_odd_gfs, k_even_gfs);
    LOOP_ALL_GFS_GPS(i) {
      const REAL k_even_gfsL = k_even_gfs[i];
      const REAL y_n_gfsL = y_n_gfs[i];
      const REAL y_nplus1_running_total_gfsL = y_nplus1_running_total_gfs[i];
      y_n_gfs[i] = (1.0/6.0)*dt*k_even_gfsL + y_n_gfsL + y_nplus1_running_total_gfsL;
    }
    if(params->outer_bc_type == EXTRAPOLATION_OUTER_BCS)
      apply_bcs_outerextrap_and_inner(params, bcstruct, y_n_gfs);
  }
  // -={ END k4 substep }=-

}
