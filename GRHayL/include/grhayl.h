#ifndef GRHAYL_H_
#define GRHAYL_H_

#include <stdbool.h>
#include <math.h>
#include "grhayl_io.h"

#ifndef MIN
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#endif
#ifndef MAX
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
#endif
#define SQR(x) ((x) * (x))
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

#ifndef GRHAYL_DISABLE_HDF5
#define GRHAYL_USE_HDF5
#endif

/*
   The struct grhayl_parameters contains parameters for controlling
   the behavior of the GRHayL gems. The struct elements are detailed below:

 --main_routine: selects which con2prim routine to use. The
   available con2prim routines are given by the enum con2prim_routines.

 --backup_routine[3]: stores up to three backup routines using the values
   in con2prim_routines. The "None" option is provided if backups are not
   desired or if there are less than 3 backup routines.

 --evolve_entropy: tells the code whether or not entropy is being evolved.
   If true, then the code will set the conservative entropy when needed.

 --calc_prim_guess: currently does nothing. However, some codes might be able to
   directly provide guesses for the initial primitive values for con2prim. If they
   can, this should be supported.
*/

typedef struct grhayl_parameters {
  int main_routine, backup_routine[3];
  bool evolve_entropy;
  bool evolve_temp;
  bool calc_prim_guess;
  double psi6threshold;
  bool Cupp_Fix;
  double Lorenz_damping_factor;
} grhayl_parameters;

void ghl_initialize_params(
      const int main,
      const int backup[3],
      const bool evolve_entropy,
      const bool evolve_temp,
      const bool calc_prim_guess,
      const double psi6threshold,
      const bool Cupp_Fix,
      const double Lorenz_damping_factor,
      grhayl_parameters *restrict params);

//--------------------------------------------------

//--------------------- EOS facets -----------------

// Maxmimum number of polytropic EOS pieces
#define MAX_EOS_PARAMS (10)

/*
   The struct primitive_quantities contains variables for storing the (point-wise)
   primitive variable data. The struct elements are detailed below:

 --rho: the baryonic density rho_b

 --press: the pressure P

 --u0: the zeroth component of the fluid four-velocity

 --v*: the 3-velocity v^i used in IllinoisGRMHD. This is defined as u^i/u^0. The other
   commonly used choice is the Valencia 3-velocity Vv^i defined as
   Vv^i = u^i/W + beta^i/lapse.


 --B*: the magnetic field TODO: give specific B definition

 --entropy: the entropy S

 --Y_e: the electron fraction Y_e

 --temp: the temperature T
*/

typedef struct primitive_quantities {
  double rho, press, eps;
  double u0, vU[3];
  double BU[3];
  double Y_e, temperature, entropy;
} primitive_quantities;

/*
   The struct eos_parameters contains information about the eos being used
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

 --rho_max, tau_max, press_max, Ye_max, temp_max, eps_max, entropy_max, W_max:
   all variables marked by "_max" are the maximum value for these quantities.

           ----------- Hybrid Equation of State -----------
 --neos: sets the number of polytropic pieces for the hybrid EOS.
   The maximum number of polytropic pieces is controlled by MAX_EOS_PARAMS below.

 --rho_ppoly: array of the density values which divide the polytropic pieces

 --Gamma_ppoly: array of the polytropic indices

 --K_ppoly: array of the adiabatic constants

 --eps_integ_const: array of the integration constants for specific internal energy

 --Gamma_th: thermal adiabatic index

         ----------- Tabulated Equation of State -----------
 --root_finding_precision: root-finding precision for table inversions

 --depsdT_threshold: this threshold is used by the Palenzuela con2prim routine
*/

typedef enum {grhayl_eos_hybrid, grhayl_eos_tabulated} grhayl_eos_t;

typedef struct eos_parameters {

  //-------------- General parameters --------------
  grhayl_eos_t eos_type;
  double rho_atm, rho_min, rho_max;
  double tau_atm;
  double press_atm, press_min, press_max;
  double W_max, inv_W_max_squared;
  //------------------------------------------------

  //------- Functions available for all EOSs -------
  void (*compute_h_and_cs2)(
        struct eos_parameters const *restrict eos,
        primitive_quantities const *restrict prims,
        double *restrict h,
        double *restrict cs2);

  //----------- Hybrid Equation of State -----------
  int neos;
  double rho_ppoly[MAX_EOS_PARAMS-1];
  double Gamma_ppoly[MAX_EOS_PARAMS];
  double K_ppoly[MAX_EOS_PARAMS];
  double eps_integ_const[MAX_EOS_PARAMS];
  double Gamma_th;

  // Function prototypes
  int  (*hybrid_find_polytropic_index)(
        const struct eos_parameters *restrict eos,
        const double rho_in);

  void (*hybrid_get_K_and_Gamma)(
        const struct eos_parameters *restrict eos,
        const double rho_in,
        double *restrict K,
        double *restrict Gamma);

  void (*hybrid_set_K_ppoly_and_eps_integ_consts)(struct eos_parameters *restrict eos);

  void (*hybrid_compute_P_cold)(
        const struct eos_parameters *restrict eos,
        const double rho_in,
        double *restrict P_cold_ptr);

  void (*hybrid_compute_P_cold_and_eps_cold)(
        const struct eos_parameters *restrict eos,
        const double rho_in,
        double *restrict P_cold_ptr,
        double *restrict eps_cold_ptr);

  void (*hybrid_compute_entropy_function)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double P,
        double *restrict S );
  //------------------------------------------------

  //---------- Tabulated Equation of State ---------
  double Y_e_atm, Y_e_min, Y_e_max;
  double T_atm, T_min, T_max;
  double eps_atm, eps_min, eps_max;
  double entropy_atm, entropy_min, entropy_max;
  double root_finding_precision;
  double depsdT_threshold;

  // Table size
  int N_rho, N_T, N_Ye;

  // Tabulated quantities
  double *restrict table_all;
  double *restrict table_logrho;
  double *restrict table_logT;
  double *restrict table_Y_e;
  double *restrict table_eps;

  // Table bounds
  double table_rho_min, table_rho_max;
  double table_T_min  , table_T_max;
  double table_Y_e_min , table_Y_e_max;
  double table_P_min  , table_P_max;
  double table_eps_min, table_eps_max;
  double table_ent_min, table_ent_max;

  // Auxiliary variables
  double energy_shift;
  double temp0, temp1;
  double dlintemp, dlintempi;
  double drholintempi;
  double dlintempyei;
  double drholintempyei;
  double dtemp, dtempi;
  double drho, drhoi;
  double dye, dyei;
  double drhotempi;
  double drhoyei;
  double dtempyei;
  double drhotempyei;

  // Function prototypes
  void (*tabulated_read_table_set_EOS_params)(
        const char *nuceos_table_name,
        struct eos_parameters *restrict eos);

  void (*tabulated_free_memory)(struct eos_parameters *restrict eos);

  void (*tabulated_compute_P_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict P);

  void (*tabulated_compute_eps_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict eps);

  void (*tabulated_compute_P_eps_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict P,
        double *restrict eps);

  void (*tabulated_compute_P_eps_S_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict P,
        double *restrict eps,
        double *restrict S);

  void (*tabulated_compute_P_eps_cs2_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict P,
        double *restrict eps,
        double *restrict cs2);

  void (*tabulated_compute_P_eps_S_cs2_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict P,
        double *restrict eps,
        double *restrict S,
        double *restrict cs2);

  void (*tabulated_compute_P_eps_depsdT_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict P,
        double *restrict eps,
        double *restrict depsdT);

  void (*tabulated_compute_P_eps_muhat_mue_mup_mun_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict P,
        double *restrict eps,
        double *restrict muhat,
        double *restrict mu_e,
        double *restrict mu_p,
        double *restrict mu_n);

  void (*tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double T,
        double *restrict muhat,
        double *restrict mu_e,
        double *restrict mu_p,
        double *restrict mu_n,
        double *restrict X_n,
        double *restrict X_p);

  void (*tabulated_compute_T_from_eps)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double eps,
        double *restrict T);

  void (*tabulated_compute_P_T_from_eps)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double eps,
        double *restrict P,
        double *restrict T);

  void (*tabulated_compute_P_cs2_T_from_eps)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double eps,
        double *restrict P,
        double *restrict cs2,
        double *restrict T);

  void (*tabulated_compute_eps_cs2_T_from_P)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double P,
        double *restrict eps,
        double *restrict cs2,
        double *restrict T);

  void (*tabulated_compute_P_T_from_S)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double S,
        double *restrict P,
        double *restrict T);

  void (*tabulated_compute_P_S_depsdT_T_from_eps)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double eps,
        double *restrict P,
        double *restrict S,
        double *restrict depsdT,
        double *restrict T);

  void (*tabulated_compute_eps_S_T_from_P)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double P,
        double *restrict eps,
        double *restrict S,
        double *restrict T);

  void (*tabulated_compute_P_eps_T_from_S)(
        const struct eos_parameters *restrict eos,
        const double rho,
        const double Y_e,
        const double S,
        double *restrict P,
        double *restrict eps,
        double *restrict T);
  //------------------------------------------------

} eos_parameters;

#ifdef __cplusplus
extern "C" {
#endif

void ghl_initialize_eos_functions(
    grhayl_eos_t const eos_type,
    eos_parameters *restrict eos );

void ghl_initialize_hybrid_eos(
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      eos_parameters *restrict eos );

void ghl_initialize_tabulated_eos(
      const char *table_path,
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      eos_parameters *restrict eos );

void ghl_initialize_hybrid_eos_functions_and_params(
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      eos_parameters *restrict eos );

void ghl_initialize_tabulated_eos_functions_and_params(
      const char *table_path,
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      eos_parameters *restrict eos );

#ifdef __cplusplus
}
#endif

//--------------------------------------------------

//--------------- Con2Prim facets ------------------

/*
   The struc conservative_quantities contains variables for storing the (point-wise)
   conservative variable data. Since most of these variables are densitized, let's
   define dens = sqrt(gamma). Then, the struct elements are detailed below:

 --rho: the densitized baryonic density \tilde{D} = rho_star = dens * lapse * rho_b * u^0

 --tau: the densitized energy density variable \tilde{tau} = dens * tau = dens * lapse^2 * T^{00} - rho_star

 --S*: the densitized momentum density \tilde{S}_i = dens * S_i = dens * lapse * T_i^0

 --Y_e: the densitized electron fraction \tilde{Y}_e = TODO

 --entropy: the densitized entropy \tilde{S} = TODO
*/

typedef struct conservative_quantities {
  double rho, tau, Y_e;
  double SD[3];
  double entropy;
} conservative_quantities;

#ifdef __cplusplus
extern "C" {
#endif

void ghl_initialize_primitives(
      const double rho, const double press, const double epsilon,
      const double vx, const double vy, const double vz,
      const double Bx, const double By, const double Bz,
      const double entropy, const double Y_e, const double temp,
      primitive_quantities *restrict prims);

void ghl_initialize_conservatives(
      const double rho, const double tau,
      const double S_x, const double S_y, const double S_z,
      const double Y_e, const double entropy,
      conservative_quantities *restrict cons);

void ghl_return_primitives(
      const primitive_quantities *restrict prims,
      double *restrict rho, double *restrict press, double *restrict epsilon,
      double *restrict vx, double *restrict vy, double *restrict vz,
      double *restrict Bx, double *restrict By, double *restrict Bz,
      double *restrict entropy, double *restrict Y_e, double *restrict temp);

void ghl_return_conservatives(
      const conservative_quantities *restrict cons,
      double *restrict rho, double *restrict tau,
      double *restrict S_x, double *restrict S_y, double *restrict S_z,
      double *restrict Y_e, double *restrict entropy);

#ifdef __cplusplus
}
#endif

//--------------------------------------------------

//-------------- Space-time facets -----------------

/* The struct metric_quantities contains variables for storing the (point-wise)
   metric data. The struct elements are detailed below:

 --bssn_phi: TODO

 --bssn_psi: TODO

 --bssn_gij: the BSSN conformal metric g_{i j}. These variables are set via the
   ghl_initialize_metric function.

 --betai: the shift vector (beta^x, beta^y, beta^z) from the 3+1
   decomposition. These variables are set via the ghl_initialize_metric function.

 --lapse: the lapse from the 3+1 decomposition. These variables are set via the
   ghl_initialize_metric function.

   All of the following varibles are calculated from the previous quantities in the ghl_initialize_metric
   function.

 --gij: the input metric g_{i j}

 --gupij: the metric inverse g^{i j}. They are calculated from
     gupxx =   ( gyy * gzz - gyz * gyz )/detgij;
     gupxy = - ( gxy * gzz - gyz * gxz )/detgij;
     gupxz =   ( gxy * gyz - gyy * gxz )/detgij;
     gupyy =   ( gxx * gzz - gxz * gxz )/detgij;
     gupyz = - ( gxx * gyz - gxy * gxz )/detgij;
     gupzz =   ( gxx * gyy - gxy * gxy )/detgij;

 --lapseinv: this variable stores the inverse of the lapse.

 --lapm1: this variable stores the quantity (lapse-1).

 --psi2: this variable stores the quantity exp(2.0*metric.bssn_phi).

 --psi4: this variable stores the quantity (psi2^2).

 --psi6: this variable stores the quantity (psi4*psi2).

 --psi4inv: this variable stores the inverse of psi4.

 --lapseinv2: this variables stores the quantity (lapse^(-2))

 --g4DD: the 4-metric g_{\mu \nu}. This quantity is needed for computing T_{\mu \nu} and T^{\mu \nu}
   and the HARM con2prim lowlevel functions.

 --g4UU: the 4-metric inverse g^{\mu \nu}. This quantity is needed for computing T_{\mu \nu} and T^{\mu \nu}
   and the HARM con2prim lowlevel functions.
*/

typedef struct metric_quantities {
  double lapse, gijdet;
  double betaU[3];
  double gammaDD[3][3];
  double gammaUU[3][3];
  double lapseinv, lapseinv2;
} metric_quantities;

typedef struct ADM_aux_quantities {
  double phi;
  double psi2, psi6;
  double psi4, psi4inv;
  double g4DD[4][4],g4UU[4][4];
} ADM_aux_quantities;

typedef struct extrinsic_curvature {
  double K[3][3];
} extrinsic_curvature;

typedef struct stress_energy {
  double T4[4][4];
} stress_energy;

#ifdef __cplusplus
extern "C" {
#endif

void ghl_initialize_metric(
      const double lapse,
      const double betax, const double betay, const double betaz,
      const double gxx, const double gxy, const double gxz,
      const double gyy, const double gyz, const double gzz,
      metric_quantities *restrict metric);

void ghl_compute_ADM_auxiliaries(
      const metric_quantities *restrict ADM_metric,
      ADM_aux_quantities *restrict metric_aux);

void ghl_enforce_detgtij_and_initialize_ADM_metric(
      const double lapse,
      const double betax, const double betay, const double betaz,
      const double gxx, const double gxy, const double gxz,
      const double gyy, const double gyz, const double gzz,
      metric_quantities *restrict ADM_metric);

void ghl_initialize_extrinsic_curvature(
      const double Kxx, const double Kxy, const double Kxz,
      const double Kyy, const double Kyz, const double Kzz,
      extrinsic_curvature *restrict curv);

void ghl_initialize_stress_energy(
      const double Ttt,
      const double Ttx, const double Tty, const double Ttz,
      const double Txx, const double Txy, const double Txz,
      const double Tyy, const double Tyz, const double Tzz,
      stress_energy *restrict Tmunu);

void ghl_return_stress_energy(
      const stress_energy *restrict Tmunu,
      double *restrict Ttt, double *restrict Ttx, double *restrict Tty,
      double *restrict Ttz, double *restrict Txx, double *restrict Txy,
      double *restrict Txz, double *restrict Tyy, double *restrict Tyz,
      double *restrict Tzz);

void ghl_limit_v_and_compute_u0(
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      primitive_quantities *restrict prims,
      int *restrict speed_limit);

void ghl_raise_vector_4D(
      const double g4UU[4][4],
      const double vecD[4],
      double vecU[4]);

void ghl_raise_vector_3D(
      const double gammaUU[3][3],
      const double vecD[3],
      double vecU[3]);

void ghl_lower_vector_4D(
      const double g4DD[4][4],
      const double vecU[4],
      double vecD[4]);

void ghl_lower_vector_3D(
      const double gammaDD[3][3],
      const double vecU[3],
      double vecD[3]);

double ghl_compute_vec2_from_vecD(
      const double gammaUU[3][3],
      const double *restrict vecD);

double ghl_compute_vec2_from_vecU(
      const double gammaDD[3][3],
      const double *restrict vecU);

void ghl_compute_TDNmunu(
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const primitive_quantities *restrict prims,
      stress_energy *restrict Tmunu);

void ghl_compute_TUPmunu(
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const primitive_quantities *restrict prims,
      stress_energy *restrict Tmunu);

void ghl_compute_smallb_and_b2(
      const metric_quantities *restrict ADM_metric,
      const primitive_quantities *restrict prims,
      const double uDN[4],
      double smallb[4],
      double *restrict smallb2);

#ifdef __cplusplus
}
#endif

#endif // GRHayL_H
