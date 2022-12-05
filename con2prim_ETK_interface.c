#include "cctk.h"
#include "con2prim_header.h"
void con2prim_ETK_kernel( CCTK_POINTER_TO_CONST void_params, CCTK_POINTER_TO_CONST void_eos,
                          CCTK_POINTER void_metric, CCTK_POINTER void_cons,
                          CCTK_POINTER void_prims, CCTK_POINTER void_diagnostics, CCTK_POINTER void_Tmunu) {

  const GRMHD_parameters *params = (const GRMHD_parameters *)void_params;
  const eos_parameters *eos = (const eos_parameters *)void_eos;
  metric_quantities *metric = (metric_quantities *)void_metric;
  conservative_quantities *cons = (conservative_quantities *)void_cons;
  primitive_quantities *prims = (primitive_quantities *)void_prims;
  con2prim_diagnostics *diagnostics = (con2prim_diagnostics *)void_diagnostics;
  stress_energy *Tmunu = (stress_energy *)void_Tmunu;

  con2prim_loop_kernel(params, eos, metric, cons, prims, diagnostics, Tmunu);
}

void c2p_initialize_conservatives(
             const CCTK_REAL rho, const CCTK_REAL tau,
             const CCTK_REAL S_x, const CCTK_REAL S_y, const CCTK_REAL S_z,
             const CCTK_REAL Y_e, const CCTK_REAL entropy,
             CCTK_POINTER void_cons) {

  conservative_quantities *cons = (conservative_quantities *)void_cons;

  initialize_conservatives(rho, tau, S_x, S_y, S_z, Y_e, entropy, cons);
}

void c2p_initialize_diagnostics(CCTK_POINTER void_diagnostics) {

  con2prim_diagnostics *diagnostics = (con2prim_diagnostics *)void_diagnostics;
  initialize_diagnostics(diagnostics);
}

void c2p_initialize_general_eos(const int type,
                const CCTK_REAL tau_atm, const CCTK_REAL W_max,
                const CCTK_REAL entropy_atm, const CCTK_REAL entropy_min, const CCTK_REAL entropy_max,
                const CCTK_REAL rho_atm, const CCTK_REAL rho_min, const CCTK_REAL rho_max,
                CCTK_POINTER void_eos) {

  eos_parameters *eos = (eos_parameters *)void_eos;

  initialize_general_eos(type, tau_atm, W_max,
                         entropy_atm, entropy_min, entropy_max,
                         rho_atm, rho_min, rho_max, eos);
}

void c2p_initialize_hybrid_eos(const int neos,
                const CCTK_REAL rho_ppoly[], const CCTK_REAL Gamma_ppoly[],
                const CCTK_REAL K_ppoly,
                const CCTK_REAL Gamma_th, CCTK_POINTER void_eos) {

  eos_parameters *eos = (eos_parameters *)void_eos;

  initialize_hybrid_eos(neos, rho_ppoly, Gamma_ppoly, K_ppoly, Gamma_th, eos);
}

void c2p_initialize_metric(const CCTK_REAL lapse,
                const CCTK_REAL gxx, const CCTK_REAL gxy, const CCTK_REAL gxz,
                const CCTK_REAL gyy, const CCTK_REAL gyz, const CCTK_REAL gzz,
                const CCTK_REAL betax, const CCTK_REAL betay, const CCTK_REAL betaz,
                CCTK_POINTER void_metric) {

  metric_quantities *metric = (metric_quantities *)void_metric;

  initialize_metric(lapse, gxx, gxy, gxz,
                    gyy, gyz, gzz, betax, betay, betaz, metric);

}

void c2p_initialize_parameters(const int main, const int backup[3],
                const int evolve_entropy, const int evolve_temp, const int calc_prim_guess,
                const CCTK_REAL psi6threshold, const int update_Tmunu, const int Cupp_Fix,
                CCTK_POINTER void_params) {

  GRMHD_parameters *params = (GRMHD_parameters *)void_params;

initialize_parameters(main, backup,
                evolve_entropy, evolve_temp, calc_prim_guess,
                psi6threshold, update_Tmunu, Cupp_Fix, params);

}

void c2p_initialize_primitives(
             const CCTK_REAL rho, const CCTK_REAL press, const CCTK_REAL epsilon,
             const CCTK_REAL vx, const CCTK_REAL vy, const CCTK_REAL vz,
             const CCTK_REAL Bx, const CCTK_REAL By, const CCTK_REAL Bz,
             const CCTK_REAL entropy, const CCTK_REAL Y_e, const CCTK_REAL temp,
             CCTK_POINTER void_prims) {

  primitive_quantities *prims = (primitive_quantities *)void_prims;

  initialize_primitives(rho, press, epsilon, vx, vy, vz,
                        Bx, By, Bz, entropy, Y_e, temp, prims);
}

void c2p_return_conservatives(
             CCTK_POINTER_TO_CONST void_cons,
             CCTK_REAL *restrict rho, CCTK_REAL *restrict tau,
             CCTK_REAL *restrict S_x, CCTK_REAL *restrict S_y, CCTK_REAL *restrict S_z,
             CCTK_REAL *restrict Y_e, CCTK_REAL *restrict entropy) {

  const conservative_quantities *cons = (const conservative_quantities *)void_cons;

  return_conservatives(cons, rho, tau,
                       S_x, S_y, S_z, Y_e, entropy);

}
void c2p_return_primitives(
             CCTK_POINTER_TO_CONST void_prims,
             CCTK_REAL *restrict rho, CCTK_REAL *restrict press, CCTK_REAL *restrict epsilon,
             CCTK_REAL *restrict vx, CCTK_REAL *restrict vy, CCTK_REAL *restrict vz,
             CCTK_REAL *restrict Bx, CCTK_REAL *restrict By, CCTK_REAL *restrict Bz,
             CCTK_REAL *restrict entropy, CCTK_REAL *restrict Y_e, CCTK_REAL *restrict temp) {

  const primitive_quantities *prims = (const primitive_quantities *)void_prims;

  return_primitives(prims, rho, press, epsilon, vx, vy, vz,
                    Bx, By, Bz, entropy, Y_e, temp);
}
