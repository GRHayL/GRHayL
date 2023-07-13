#include "reconstruction.h"

void ghl_ppm_Ur_Ul(
      const double rho[5],
      const double pressure[5],
      const double var_data[][5],
      const int num_vars,
      const double v_flux_dirn[5],
      const double Gamma_eff,
      double *restrict rhor,
      double *restrict rhol,
      double *restrict pressr,
      double *restrict pressl,
      double *restrict varsr,
      double *restrict varsl);

/* Function    : ghl_ppm()
 * Description : reconstructs rho_b, pressure, and other variables at the points
 *                   Ur(i) = U(i-1/2+epsilon)
 *                   Ul(i) = U(i-1/2-epsilon)
 *               input stencils require primitives in stencil (i-3,i-2,i-1,i,i+1,i+2)
 *               along reconstruction direction
 *
 * Inputs      : rho[6]         - stencil for baryonic density
 *             : pressure[6]    - stencil for baryonic density
 *             : var_data[][6]  - array of pointers to stencils for other reconstructed
 *                                variables; has length num_vars
 *             : num_vars       - number of other variables to reconstruct
 *             : v_flux_dirn[6] - stencil of velocity along reconstruction direction
 *             : Gamma_eff      - effective Gamma = (partial P / partial rho0)_s /(P/rho0)
 *
 * Outputs     : rhor           - reconstructed rho_b on right face
 *             : rhol           - reconstructed rho_b on left face
 *             : pressr         - reconstructed pressure on right face
 *             : pressl         - reconstructed pressure on left face
 *             : var_datar      - array of reconstructed other variables on right face
 *             : var_datal      - array of reconstructed other variables on left face
 */

void ghl_ppm(
      const double rho[6],
      const double pressure[6],
      const double var_data[][6],
      const int num_vars,
      const double v_flux_dirn[6],
      const double Gamma_eff,
      double *restrict rhor,
      double *restrict rhol,
      double *restrict pressr,
      double *restrict pressl,
      double *restrict var_datar,
      double *restrict var_datal) {

  // For the variables array, the data slices need to be
  // explicitly extracted to pass to the PPM routine.
  double tmp_data[num_vars][5];
  for(int var=0; var<num_vars; var++)
    for(int i=0; i<5; i++)
    tmp_data[var][i] = var_data[var][i+1];
  double tmp_rhor, tmp_rhol;
  double tmp_pressr, tmp_pressl;
  double tmp_varsr[num_vars], tmp_varsl[num_vars];

  // ghl_ppm_Ur_Ul evaluates
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
  ghl_ppm_Ur_Ul(&rho[1], &pressure[1], tmp_data, num_vars,
            &v_flux_dirn[1], Gamma_eff,
            &tmp_rhor, &tmp_rhol, &tmp_pressr, &tmp_pressl,
            tmp_varsr, tmp_varsl);

  // tmp_Ul[PLUS_0] is Ur[PLUS0], so set that now:
  *rhor = tmp_rhol;
  *pressr = tmp_pressl;
  for(int var=0;var<num_vars;var++) {
    var_datar[var] = tmp_varsl[var];
    // Setting the data slice for Step 2
    for(int i=0; i<5; i++)
    tmp_data[var][i] = var_data[var][i];
  }

  // STEP 2: Evaluate Ur[MINUS1] and Ul[MINUS1],
  //         which depend on U[0],U[1],U[2],U[3],U[4]
  ghl_ppm_Ur_Ul(rho, pressure, tmp_data, num_vars,
            v_flux_dirn, Gamma_eff,
            &tmp_rhor, &tmp_rhol, &tmp_pressr, &tmp_pressl,
            tmp_varsr, tmp_varsl);

  // tmp_Ur[MINUS1] is Ul[PLUS0], so set that now:
  *rhol = tmp_rhor;
  *pressl = tmp_pressr;
  for(int var=0;var<num_vars;var++)
    var_datal[var] = tmp_varsr[var];
}

// Inputs: Primitives U at *five* locations: i-2,i-1,i,i+1,i+2,
//         where i-2 == MINUS2; i-1 == MINUS1; i == PLUS_0, etc.
// Outputs: tmp_Ur[PLUS_0] = U(i+1/2)
//          tmp_Ul[PLUS_0] = U(i-1/2)
void ghl_ppm_Ur_Ul(
      const double rho[5],
      const double pressure[5],
      const double var_data[][5],
      const int num_vars,
      const double v_flux_dirn[5],
      const double Gamma_eff, // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
      double *restrict rhor,
      double *restrict rhol,
      double *restrict pressr,
      double *restrict pressl,
      double *restrict varsr,
      double *restrict varsl) {

  // Interpolate primitives to faces with a slope limiter.
  ghl_compute_UrUl_onevar(rho, rhor, rhol);
  ghl_compute_UrUl_onevar(pressure, pressr, pressl);
  for(int var=0;var<num_vars;var++) {
    ghl_compute_UrUl_onevar(var_data[var], &varsr[var], &varsl[var]);
  }

  // Steepen rhol and rhor
  ghl_steepen_rhor_rhol(rho, pressure, Gamma_eff, rhor, rhol);

  // Flatten all variables
  const double ftilde = ghl_shock_detection_ftilde(pressure, v_flux_dirn);
  ghl_flatten_and_monotonize_Ur_and_Ul(rho[PLUS_0], ftilde, rhor, rhol);
  ghl_flatten_and_monotonize_Ur_and_Ul(pressure[PLUS_0], ftilde, pressr, pressl);
  for(int var=0;var<num_vars;var++)
    ghl_flatten_and_monotonize_Ur_and_Ul(var_data[var][PLUS_0], ftilde, &varsr[var], &varsl[var]);
}
