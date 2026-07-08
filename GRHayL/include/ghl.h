#ifndef GHL_H_
#define GHL_H_

#include <stdbool.h>
#include <math.h>
#include "ghl_io.h"
#include "ghl_metric_helpers.h"

#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

#ifndef ONE_OVER_SQRT_4PI
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
extern "C" {
#endif

typedef enum {
  ghl_success,
  ghl_error_unknown_eos_type,
  ghl_error_invalid_c2p_key,
  ghl_error_neg_rho,
  ghl_error_neg_pressure,
  ghl_error_neg_vsq,
  ghl_error_c2p_max_iter,
  ghl_error_c2p_singular,
  ghl_error_root_not_bracketed,
  ghl_error_table_max_rho,
  ghl_error_table_min_rho,
  ghl_error_table_max_ye,
  ghl_error_table_min_ye,
  ghl_error_table_max_T,
  ghl_error_table_min_T,
  ghl_error_exceed_table_vars,
  ghl_error_table_neg_energy,
  ghl_error_table_bisection,
  ghl_error_u0_singular,
  ghl_error_invalid_utsq,
  ghl_error_invalid_Z,
  ghl_error_newman_invalid_discriminant,
  ghl_error_used_disabled_hdf5,
  ghl_error_out_of_memory,
  ghl_error_eos_struct_is_null,
  ghl_error_invalid_eos_type,
  ghl_error_invalid_eos_table_type,
  ghl_error_could_not_open_file,
  ghl_error_hdf5_dataset_could_not_open,
  ghl_error_hdf5_dataset_could_not_read,
  ghl_error_hdf5_dataset_invalid_ndims,
  ghl_error_hdf5_dataset_size_mismatch,
  ghl_error_invalid_rho_atm,
  ghl_error_rho_min_gt_rho_max,
  ghl_error_invalid_press_atm,
  ghl_error_press_min_gt_press_max,
  ghl_error_invalid_Y_e_atm,
  ghl_error_Y_e_min_gt_Y_e_max,
  ghl_error_invalid_T_atm,
  ghl_error_T_min_gt_T_max,
  ghl_error_invalid_fermi_dirac_integral_key,
  ghl_error_nn_c2p_model_is_null,
  ghl_error_nn_c2p_invalid_dimensions,
  ghl_error_nn_c2p_invalid_input_index,
  ghl_error_nn_c2p_missing_array,
  ghl_error_nn_c2p_invalid_kind,
  ghl_error_nn_c2p_invalid_number,
} ghl_error_codes_t;

typedef enum {
  ghl_con2prim_id_None = -1,
  ghl_con2prim_id_Noble2D,
  ghl_con2prim_id_Noble1D,
  ghl_con2prim_id_Noble1D_entropy,
  ghl_con2prim_id_Noble1D_entropy2,
  ghl_con2prim_id_Font1D,
  ghl_con2prim_id_Palenzuela1D,
  ghl_con2prim_id_Palenzuela1D_entropy,
  ghl_con2prim_id_Newman1D,
  ghl_con2prim_id_Newman1D_entropy,
} ghl_con2prim_id_t;

typedef enum {
  ghl_eos_simple,
  ghl_eos_hybrid,
  ghl_eos_tabulated,
} ghl_eos_t;

typedef enum {
  ghl_eos_table_unknown,
  ghl_eos_table_stellarcollapse,
  ghl_eos_table_types,
} ghl_eos_table_t;

typedef struct ghl_c2p_nn_model ghl_c2p_nn_model;

/*
 * Struct        : ghl_parameters
 * Description   : stores basic GRHayL parameters
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_parameters
*/
typedef struct ghl_parameters {
  ghl_con2prim_id_t main_routine, backup_routine[3];
  bool evolve_entropy;
  bool evolve_temp;
  bool calc_prim_guess;
  double max_Lorentz_factor;
  double inv_sq_max_Lorentz_factor;
  double psi6threshold;
  double Lorenz_damping_factor;

  // Con2Prim parameters
  int con2prim_max_iterations;
  double con2prim_solver_tolerance;

  // PPM parameters
  double ppm_flattening_epsilon;
  double ppm_flattening_omega1;
  double ppm_flattening_omega2;
  double ppm_shock_k0;
  double ppm_shock_eta1;
  double ppm_shock_eta2;
  double ppm_shock_epsilon;
} ghl_parameters;

/*
 * Struct        : ghl_primitive_quantities
 * Description   : stores pointwise information about the primitive variables
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_primitive_quantities
*/
typedef struct ghl_primitive_quantities {
  double rho, press, eps;
  double u0, vU[3];
  double BU[3];
  double Y_e, temperature, entropy;
} ghl_primitive_quantities;

/*
 * Struct        : ghl_conservative_quantities
 * Description   : stores pointwise information about the conservative variables
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_conservative_quantities
*/
typedef struct ghl_conservative_quantities {
  double rho, tau, Y_e;
  double SD[3];
  double entropy;
} ghl_conservative_quantities;

/*
 * Struct        : ghl_metric_quantities
 * Description   : stores pointwise information about the spacetime
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_metric_quantities
*/
typedef struct ghl_metric_quantities {
  double lapse, lapseinv, lapseinv2;
  double detgamma, sqrt_detgamma;
  double betaU[3];
  double gammaDD[3][3];
  double gammaUU[3][3];
} ghl_metric_quantities;

/*
 * Struct        : ghl_ADM_aux_quantities
 * Description   : stores auxiliary information based on ADM metric quantities
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_ADM_aux_quantities
*/
typedef struct ghl_ADM_aux_quantities {
  double g4DD[4][4], g4UU[4][4];
} ghl_ADM_aux_quantities;

/*
 * Struct        : ghl_extrinsic_curvature
 * Description   : stores pointwise information about the extrinsic curvature
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_extrinsic_curvature
*/
typedef struct ghl_extrinsic_curvature {
  double K[3][3];
} ghl_extrinsic_curvature;

/*
 * Struct        : ghl_stress_energy
 * Description   : stores pointwise information about the stress energy tensor
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_stress_energy
*/
typedef struct ghl_stress_energy {
  double T4[4][4];
} ghl_stress_energy;

// Maxmimum number of polytropic EOS pieces
#define MAX_EOS_PARAMS (10)
/*
   The struct ghl_eos_parameters contains information about the eos being used
   by the simulation. The struct elements are detailed below:

 --type: selects the type of EOS (hybrid or tabulated) where
   Hybrid = 0, Tabulated = 1.

 --rho_atm, tau_atm, press_atm, Ye_atm, temp_atm, eps_atm, entropy_atm: all
   variables marked by "_atm" are the values for those quantities at
   atmosphere.

 --rho_min, tau_min, press_min, Ye_min, temp_min, eps_min, entropy_min: all
   variables marked by "_min" are the minimum value for these quantities.
   This is often just the atmospheric value. If the simulation does not have
   a separate minimum value, simply pass the atmospheric value for these
   parameters.

 --rho_max, tau_max, press_max, Ye_max, temp_max, eps_max, entropy_max:
   all variables marked by "_max" are the maximum value for these quantities.

           ----------- Hybrid Equation of State -----------
 --neos: sets the number of polytropic pieces for the hybrid EOS.
   The maximum number of polytropic pieces is controlled by MAX_EOS_PARAMS below.

 --rho_ppoly: array of the density values which divide the polytropic pieces

 --Gamma_ppoly: array of the polytropic indices

 --K_ppoly: array of the adiabatic constants

 --eps_integ_const: array of the integration constants for specific internal energy

 --p_ppoly: array of pressure values which divide the polytropic pieces

 --Gamma_th: thermal adiabatic index

         ----------- Tabulated Equation of State -----------
 --root_finding_precision: root-finding precision for table inversions
*/

typedef struct ghl_eos_parameters {

  //-------------- General parameters --------------
  ghl_eos_t eos_type;
  ghl_eos_table_t table_type;
  bool clean_sound_speed;
  double rho_atm, rho_min, rho_max;
  double tau_atm;
  double press_atm, press_min, press_max;
  //------------------------------------------------

  //----------- Hybrid Equation of State -----------
  int neos;
  double rho_ppoly[MAX_EOS_PARAMS - 1];
  double Gamma_ppoly[MAX_EOS_PARAMS];
  double K_ppoly[MAX_EOS_PARAMS];
  double eps_integ_const[MAX_EOS_PARAMS];
  double p_ppoly[MAX_EOS_PARAMS];
  double Gamma_th;
  //------------------------------------------------

  //---------- Tabulated Equation of State ---------
  double Y_e_atm, Y_e_min, Y_e_max;
  double T_atm, T_min, T_max;
  double eps_atm, eps_min, eps_max;
  double entropy_atm, entropy_min, entropy_max;
  double root_finding_precision;

  // Table size
  int N_rho, N_T, N_Ye;

  // Tabulated quantities
  double *restrict table_all;
  double *restrict table_logrho;
  double *restrict table_logT;
  double *restrict table_Y_e;
  double *restrict table_eps;
  double *restrict table_logh;

  // Table bounds
  double table_rho_min, table_rho_max;
  double table_T_min, table_T_max;
  double table_Y_e_min, table_Y_e_max;
  double table_P_min, table_P_max;
  double table_eps_min, table_eps_max;
  double table_ent_min, table_ent_max;

  // Auxiliary variables
  double energy_shift;
  double dtempi;
  double drhoi;
  double dyei;
  double drhotempi;
  double drhoyei;
  double dtempyei;
  double drhotempyei;

  // These are used for beta-equilibrium
  double *lp_of_lr, *le_of_lr, *Ye_of_lr, *lh_of_lr;

  bool enable_neural_net_c2p;
  ghl_c2p_nn_model *restrict c2p_nn;
  //------------------------------------------------
} ghl_eos_parameters;

const char *ghl_get_con2prim_routine_name(const ghl_con2prim_id_t key);

void ghl_initialize_eos_functions(
    const ghl_eos_t eos_type);

ghl_error_codes_t ghl_initialize_simple_eos(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *restrict eos);

ghl_error_codes_t ghl_initialize_hybrid_eos(
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *restrict eos);

ghl_error_codes_t ghl_initialize_tabulated_eos(
      const char *table_path,
      const ghl_eos_table_t table_type,
      const bool clean_sound_speed,
      const bool enable_neural_net_c2p,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      ghl_eos_parameters *restrict eos);

ghl_error_codes_t ghl_initialize_simple_eos_functions_and_params(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *restrict eos);

ghl_error_codes_t ghl_initialize_hybrid_eos_functions_and_params(
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *restrict eos);

ghl_error_codes_t ghl_initialize_tabulated_eos_functions_and_params(
      const char *table_path,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      ghl_eos_parameters *restrict eos);

//---- Basic struct packing/unpacking functions ----
void ghl_initialize_params(
      const ghl_con2prim_id_t main_routine,
      const ghl_con2prim_id_t backup_routine[3],
      const bool evolve_entropy,
      const bool evolve_temp,
      const bool calc_prim_guess,
      const double psi6threshold,
      const double max_Lorentz_factor,
      const double Lorenz_damping_factor,
      ghl_parameters *restrict params);

void ghl_initialize_primitives(
      const double rho,
      const double press,
      const double epsilon,
      const double vx,
      const double vy,
      const double vz,
      const double Bx,
      const double By,
      const double Bz,
      const double entropy,
      const double Y_e,
      const double temp,
      ghl_primitive_quantities *restrict prims);

void ghl_initialize_conservatives(
      const double rho,
      const double tau,
      const double S_x,
      const double S_y,
      const double S_z,
      const double entropy,
      const double Y_e,
      ghl_conservative_quantities *restrict cons);

void ghl_return_primitives(
      const ghl_primitive_quantities *restrict prims,
      double *restrict rho,
      double *restrict press,
      double *restrict epsilon,
      double *restrict vx,
      double *restrict vy,
      double *restrict vz,
      double *restrict Bx,
      double *restrict By,
      double *restrict Bz,
      double *restrict entropy,
      double *restrict Y_e,
      double *restrict temp);

void ghl_return_conservatives(
      const ghl_conservative_quantities *restrict cons,
      double *restrict rho,
      double *restrict tau,
      double *restrict S_x,
      double *restrict S_y,
      double *restrict S_z,
      double *restrict entropy,
      double *restrict Y_e);

void ghl_initialize_metric(
      const double lapse,
      const double betax,
      const double betay,
      const double betaz,
      const double gxx,
      const double gxy,
      const double gxz,
      const double gyy,
      const double gyz,
      const double gzz,
      ghl_metric_quantities *restrict metric);

void ghl_compute_ADM_auxiliaries(
      const ghl_metric_quantities *restrict metric_adm,
      ghl_ADM_aux_quantities *restrict metric_aux);

void ghl_enforce_detgtij_and_initialize_ADM_metric(
      const double lapse,
      const double betax,
      const double betay,
      const double betaz,
      const double gxx,
      const double gxy,
      const double gxz,
      const double gyy,
      const double gyz,
      const double gzz,
      ghl_metric_quantities *restrict metric_adm);

void ghl_initialize_extrinsic_curvature(
      const double Kxx,
      const double Kxy,
      const double Kxz,
      const double Kyy,
      const double Kyz,
      const double Kzz,
      ghl_extrinsic_curvature *restrict curv);

void ghl_initialize_stress_energy(
      const double Ttt,
      const double Ttx,
      const double Tty,
      const double Ttz,
      const double Txx,
      const double Txy,
      const double Txz,
      const double Tyy,
      const double Tyz,
      const double Tzz,
      ghl_stress_energy *restrict Tmunu);

void ghl_return_stress_energy(
      const ghl_stress_energy *restrict Tmunu,
      double *restrict Ttt,
      double *restrict Ttx,
      double *restrict Tty,
      double *restrict Ttz,
      double *restrict Txx,
      double *restrict Txy,
      double *restrict Txz,
      double *restrict Tyy,
      double *restrict Tyz,
      double *restrict Tzz);

ghl_error_codes_t ghl_limit_v_and_compute_u0(
      const ghl_parameters *restrict params,
      const ghl_metric_quantities *restrict metric_adm,
      ghl_primitive_quantities *restrict prims,
      bool *restrict speed_limited);

void ghl_compute_TDNmunu(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_stress_energy *restrict Tmunu);

void ghl_compute_TUPmunu(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_stress_energy *restrict Tmunu);

void ghl_compute_smallb_and_b2(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_primitive_quantities *restrict prims,
      const double uDN[4],
      double smallb[4],
      double *restrict smallb2);

void ghl_compute_SU_Bsq_Ssq_BdotS(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      double *restrict SU,
      double *restrict Bsq,
      double *restrict Ssq,
      double *restrict BdotS);

void ghl_abort_if_error(const ghl_error_codes_t error);

int ghl_imin(int a, int b);

int ghl_imax(int a, int b);

int ghl_iclamp(int x, int x_min, int x_max);

double ghl_clamp(double x, double x_min, double x_max);

#ifdef __cplusplus
}
#endif

#include "ghl_eos_functions.h"
#include "ghl_debug.h"

#endif // GHL_H
