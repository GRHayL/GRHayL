#include "reconstruction.h"

void grhayl_ppm_Ur_Ul_no_rho_P(
      const double pressure[5],
      const double var_data[][5],
      const int num_vars,
      const double v_flux_dirn[5],
      const double Gamma_eff,
      double *restrict varsr,
      double *restrict varsl);

/* Function    : grhayl_simple_ppm_no_rho_P()
 * Description : reconstructs variables at the points
 *                   Ur(i) = U(i-1/2+epsilon)
 *                   Ul(i) = U(i-1/2-epsilon)
 *               input stencils require primitives in stencil (i-3,i-2,i-1,i,i+1,i+2)
 *               along reconstruction direction; much is identical to grhayl_simple_ppm(), but
 *               rho_b has extra function call; this function is a simpler version for when
 *               rho_b and pressure do not need to be reconstructed
 *
 * Inputs      : pressure[6]    - stencil for pressure; used to compute ftilde for flattening
 *             : var_data[][6]  - array of pointers to stencils for reconstructed
 *                                variables; has length num_vars
 *             : num_vars       - number of variables to reconstruct
 *             : v_flux_dirn[6] - stencil of velocity along reconstruction direction
 *             : Gamma_eff      - effective Gamma = (partial P / partial rho0)_s /(P/rho0)
 *
 * Outputs     : var_datar      - array of reconstructed variables on right face
 *             : var_datal      - array of reconstructed variables on left face
 */
void grhayl_simple_ppm_no_rho_P(
      const double pressure[6],
      const double var_data[][6],
      const int num_vars,
      const double v_flux_dirn[6],
      const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
      double *restrict var_datar,
      double *restrict var_datal) {

  // For the variables array, the data slices need to be
  // explicitly extracted to pass to the PPM routine.
  double tmp_data[num_vars][5];
  for(int var=0; var<num_vars; var++)
    for(int i=0; i<5; i++)
    tmp_data[var][i] = var_data[var][i+1];
  double tmp_varsr[num_vars], tmp_varsl[num_vars];

  // grhayl_ppm_Ur_Ul evaluates
  //  * tmp_Ur[PLUS_0] = U(i+1/2)
  //  * tmp_Ul[PLUS_0] = U(i-1/2)
  // However, we want
  //  * (STEP 1) Ur[PLUS_0] = U(i-1/2+epsilon) = tmp_Ul[PLUS_0]
  //  AND 
  //  * (STEP 2) Ul[PLUS_0] = U(i-1/2-epsilon) = tmp_Ur[MINUS1]

  // STEP 1: Evaluate Ur[PLUS_0] and Ul[PLUS_0],
  //         which depend on U[1],U[2],U[3],U[4],U[5],
  //         hence the passing of the address U[1]
  //         as the lower bound of each U array.
  grhayl_ppm_Ur_Ul_no_rho_P(&pressure[1], tmp_data, num_vars,
            &v_flux_dirn[1], Gamma_eff,
            tmp_varsr, tmp_varsl);

  // tmp_Ul[PLUS_0] is Ur[PLUS0], so set that now:
  for(int var=0;var<num_vars;var++) {
    var_datar[var] = tmp_varsl[var];
    // Setting the data slice for Step 2
    for(int i=0; i<5; i++)
    tmp_data[var][i] = var_data[var][i];
  }

  // STEP 2: Evaluate Ur[MINUS1] and Ul[MINUS1],
  //         which depend on U[0],U[1],U[2],U[3],U[4]
  grhayl_ppm_Ur_Ul_no_rho_P(pressure, tmp_data, num_vars,
            v_flux_dirn, Gamma_eff,
            tmp_varsr, tmp_varsl);

  // tmp_Ur[MINUS1] is Ul[PLUS0], so set that now:
  for(int var=0;var<num_vars;var++)
    var_datal[var] = tmp_varsr[var];
}

// Inputs: Primitives U at *five* locations: i-2,i-1,i,i+1,i+2,
//         where i-2 == MINUS2; i-1 == MINUS1; i == PLUS_0, etc.
// Outputs: tmp_Ur[PLUS_0] = U(i+1/2)
//          tmp_Ul[PLUS_0] = U(i-1/2)
void grhayl_ppm_Ur_Ul_no_rho_P(
      const double pressure[5],
      const double var_data[][5],
      const int num_vars,
      const double v_flux_dirn[5],
      const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
      double *restrict varsr,
      double *restrict varsl) {

  // Interpolate primitives to faces.
  for(int var=0;var<num_vars;var++) {
    grhayl_compute_UrUl_onevar(var_data[var], &varsr[var], &varsl[var]);
  }

  // Flatten all variables
  const double ftilde = grhayl_shock_detection_ftilde(pressure, v_flux_dirn);
  for(int var=0;var<num_vars;var++)
    grhayl_flatten_and_monotonize_Ur_and_Ul(var_data[var][PLUS_0], ftilde, &varsr[var], &varsl[var]);
}
