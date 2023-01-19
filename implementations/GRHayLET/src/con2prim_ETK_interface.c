#include "cctk.h"
#include "con2prim.h"

void c2p_initialize_conservatives(
             const CCTK_REAL rho, const CCTK_REAL tau,
             const CCTK_REAL S_x, const CCTK_REAL S_y, const CCTK_REAL S_z,
             const CCTK_REAL Y_e, const CCTK_REAL entropy,
             CCTK_POINTER cons) {

  initialize_conservatives(rho, tau, S_x, S_y, S_z, Y_e, entropy,
                           (conservative_quantities *)cons);
}

void c2p_initialize_diagnostics(CCTK_POINTER diagnostics) {

  initialize_diagnostics((con2prim_diagnostics *)diagnostics);
}

void c2p_initialize_general_eos(const int type,
                const CCTK_REAL W_max,
                const CCTK_REAL rho_atm, const CCTK_REAL rho_min, const CCTK_REAL rho_max,
                CCTK_POINTER eos) {
  initialize_general_eos(type, W_max, rho_atm, rho_min, rho_max, (eos_parameters *)eos);
}

void c2p_initialize_hybrid_eos(const int neos,
                const CCTK_REAL rho_ppoly[], const CCTK_REAL Gamma_ppoly[],
                const CCTK_REAL K_ppoly,
                const CCTK_REAL Gamma_th, CCTK_POINTER eos) {

  initialize_hybrid_eos(neos, rho_ppoly, Gamma_ppoly, K_ppoly, Gamma_th,
                        (eos_parameters *)eos);
}

void c2p_initialize_metric(const CCTK_REAL lapse,
                const CCTK_REAL gxx, const CCTK_REAL gxy, const CCTK_REAL gxz,
                const CCTK_REAL gyy, const CCTK_REAL gyz, const CCTK_REAL gzz,
                const CCTK_REAL betax, const CCTK_REAL betay, const CCTK_REAL betaz,
                CCTK_POINTER metric) {

  initialize_metric(lapse, gxx, gxy, gxz,
                    gyy, gyz, gzz, betax, betay, betaz,
                    (metric_quantities *)metric);

}

void c2p_initialize_parameters(const int main, const int *restrict backup,
                const int evolve_entropy, const int evolve_temp, const int calc_prim_guess,
                const CCTK_REAL psi6threshold, const int update_Tmunu, const int Cupp_Fix,
                CCTK_POINTER params) {

initialize_GRHayL(main, backup,
                evolve_entropy, evolve_temp, calc_prim_guess,
                psi6threshold, update_Tmunu, Cupp_Fix,
                (GRHayL_parameters *)params);

}

void c2p_initialize_primitives(
             const CCTK_REAL rho, const CCTK_REAL press, const CCTK_REAL epsilon,
             const CCTK_REAL vx, const CCTK_REAL vy, const CCTK_REAL vz,
             const CCTK_REAL Bx, const CCTK_REAL By, const CCTK_REAL Bz,
             const CCTK_REAL entropy, const CCTK_REAL Y_e, const CCTK_REAL temp,
             CCTK_POINTER prims) {

  initialize_primitives(rho, press, epsilon, vx, vy, vz,
                        Bx, By, Bz, entropy, Y_e, temp,
                        (primitive_quantities *)prims);
}

void c2p_return_conservatives(
             CCTK_POINTER_TO_CONST cons,
             CCTK_REAL *restrict rho, CCTK_REAL *restrict tau,
             CCTK_REAL *restrict S_x, CCTK_REAL *restrict S_y, CCTK_REAL *restrict S_z,
             CCTK_REAL *restrict Y_e, CCTK_REAL *restrict entropy) {

  return_conservatives( (const conservative_quantities *)cons,
                        rho, tau, S_x, S_y, S_z, Y_e, entropy );

}
void c2p_return_primitives(
             CCTK_POINTER_TO_CONST prims,
             CCTK_REAL *restrict rho, CCTK_REAL *restrict press, CCTK_REAL *restrict epsilon,
             CCTK_REAL *restrict vx, CCTK_REAL *restrict vy, CCTK_REAL *restrict vz,
             CCTK_REAL *restrict Bx, CCTK_REAL *restrict By, CCTK_REAL *restrict Bz,
             CCTK_REAL *restrict entropy, CCTK_REAL *restrict Y_e, CCTK_REAL *restrict temp) {

  return_primitives( (const primitive_quantities *)prims,
                     rho, press, epsilon, vx, vy, vz,
                     Bx, By, Bz, entropy, Y_e, temp );
}
