#include "cctk.h"
#include "con2prim_header.h"
void con2prim_ETK_kernel( CCTK_POINTER_TO_CONST const void_params, CCTK_POINTER_TO_CONST const void_eos,
                          CCTK_POINTER *void_metric, CCTK_POINTER *void_cons,
                          CCTK_POINTER *void_prims, CCTK_POINTER *void_diagnostics) {

  const GRMHD_parameters *params = (const GRMHD_parameters *)void_params;
  const eos_parameters *eos = (const eos_parameters *)void_eos;
  metric_quantities *metric = *( (metric_quantities **)void_metric );
  conservative_quantities *cons = *( (conservative_quantities **)void_cons );
  primitive_quantities *prims = *( (primitive_quantities **)void_prims );
  con2prim_diagnostics *diagnostics = *( (con2prim_diagnostics **)void_diagnostics );

  con2prim_loop_kernel(params, eos, metric, cons, prims, diagnostics);
}

void c2p_initialize_conservatives(
             CCTK_POINTER_TO_CONST const void_params,
             CCTK_POINTER_TO_CONST const void_eos,
             const double rho, const double tau,
             const double S_x, const double S_y, const double S_z,
             const double Y_e, const double entropy,
             CCTK_POINTER *void_cons) {

  const GRMHD_parameters *params = (const GRMHD_parameters *)void_params;
  const eos_parameters *eos = (const eos_parameters *)void_eos;
  conservative_quantities *cons = *( (conservative_quantities **)void_cons );

  initialize_conservatives(params, eos, rho, tau, S_x, S_y, S_z, Y_e, entropy, cons);
}

void c2p_initialize_diagnostics(CCTK_POINTER *void_diagnostics) {

  con2prim_diagnostics *diagnostics = *( (con2prim_diagnostics **)void_diagnostics );
  initialize_diagnostics(diagnostics);
}

void c2p_initialize_general_eos(CCTK_POINTER *void_eos, const int type,
                const double tau_atm, const double W_max,
                const double eps_atm, const double eps_min, const double eps_max,
                const double press_atm, const double press_min, const double press_max,
                const double entropy_atm, const double entropy_min, const double entropy_max,
                const double rho_atm, const double rho_min, const double rho_max) {

  eos_parameters *eos = *( (eos_parameters **)void_eos );

  initialize_general_eos(eos, type, tau_atm, W_max, eps_atm, eps_min, eps_max,
                             press_atm, press_min, press_max, entropy_atm, entropy_min,
                             entropy_max, rho_atm, rho_min, rho_max);
}

void c2p_initialize_hybrid_eos(CCTK_POINTER *void_eos, const int neos,
                const double rho_ppoly_tab[], const double Gamma_ppoly_tab[],
                const double K_ppoly_tab[], const double eps_integ_const[],
                const double Gamma_th) {

  eos_parameters *eos = *( (eos_parameters **)void_eos );

  initialize_hybrid_eos(eos, neos, rho_ppoly_tab, Gamma_ppoly_tab, K_ppoly_tab, eps_integ_const, Gamma_th);
}

void c2p_initialize_metric(CCTK_POINTER *void_metric, 
                const double phi, const double psi, const double lapse,
                const double gxx, const double gxy, const double gxz,
                const double gyy, const double gyz, const double gzz,
                const double betax, const double betay, const double betaz) {

  metric_quantities *metric = *( (metric_quantities **)void_metric );

  initialize_metric(metric, phi, psi, lapse, gxx, gxy, gxz,
                    gyy, gyz, gzz, betax, betay, betaz);

}

void c2p_initialize_parameters(CCTK_POINTER *void_params, const int main, const int backup[3],
                const int evolve_entropy, const int evolve_temp, const int calc_prim_guess,
                const double psi6threshold, const int update_Tmunu) {

  GRMHD_parameters *params = * ( (GRMHD_parameters **)void_params );

initialize_parameters(params, main, backup,
                evolve_entropy, evolve_temp, calc_prim_guess,
                psi6threshold, update_Tmunu);

}

void c2p_initialize_primitives(
             CCTK_POINTER_TO_CONST const void_eos,
             CCTK_POINTER_TO_CONST const void_metric,
             const double rho, const double press, const double epsilon,
             const double vx, const double vy, const double vz,
             const double Bx, const double By, const double Bz,
             const double entropy, const double Y_e, const double temp,
             CCTK_POINTER *void_prims) {

  const eos_parameters *eos = (const eos_parameters *)void_eos;
  const metric_quantities *metric = (metric_quantities *)void_metric;
  primitive_quantities *prims = *( (primitive_quantities **)void_prims );

  initialize_primitives(eos, metric, rho,  press, epsilon, vx, vy, vz,
                        Bx, By, Bz, entropy, Y_e, temp, prims);
}

void c2p_return_conservatives(
             CCTK_POINTER_TO_CONST const void_params,
             CCTK_POINTER_TO_CONST const void_eos,
             CCTK_POINTER_TO_CONST const void_cons,
             double *restrict rho, double *restrict tau,
             double *restrict S_x, double *restrict S_y, double *restrict S_z,
             double *restrict Y_e, double *restrict entropy) {

  const GRMHD_parameters *params = (const GRMHD_parameters *)void_params;
  const eos_parameters *eos = (const eos_parameters *)void_eos;
  const conservative_quantities *cons = (const conservative_quantities *)void_cons;


  return_conservatives(params, eos, cons, rho, tau,
                       S_x, S_y, S_z, Y_e, entropy);

}
void c2p_return_primitives(
             CCTK_POINTER_TO_CONST const void_eos,
             CCTK_POINTER_TO_CONST const void_prims,
             double *restrict rho, double *restrict press, double *restrict epsilon,
             double *restrict vx, double *restrict vy, double *restrict vz,
             double *restrict Bx, double *restrict By, double *restrict Bz,
             double *restrict entropy, double *restrict Y_e, double *restrict temp) {

  const eos_parameters *eos = (const eos_parameters *)void_eos;
  const primitive_quantities *prims = (const primitive_quantities *)void_prims;

  return_primitives(eos, prims, rho, press, epsilon, vx, vy, vz,
                    Bx, By, Bz, entropy, Y_e, temp);
}
