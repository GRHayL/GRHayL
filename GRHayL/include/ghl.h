#ifndef GHL_H_
#define GHL_H_

#include <stdbool.h>
#include <math.h>
#include "ghl_io.h"
#include "ghl_metric_helpers.h"

/** @cond DOXYGEN_IGNORE */
#ifndef MIN
#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#endif
#ifndef MAX
#define MAX(A, B) (((A) > (B)) ? (A) : (B))
#endif
#define SQR(x) ((x) * (x))
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

#ifndef GRHAYL_DISABLE_HDF5
#define GRHAYL_USE_HDF5
#endif

#ifdef __cplusplus
#ifndef restrict
#define restrict __restrict__
#endif
#endif

static inline int ghl_iclamp(
      int a,
      int lo,
      int hi) {
  return MAX( MIN(a, hi), lo);
}
static inline float ghl_fclamp(
      float a,
      float lo,
      float hi) {
  return fmaxf( fminf(a, hi), lo);
}
static inline double ghl_clamp(
      double a,
      double lo,
      double hi) {
  return fmax( fmin(a, hi), lo);
}
/** @endcond */

/**
 * @ingroup GRHayL_Core
 * @enum ghl_error_codes_t
 * @brief Integer constants to track error types
 *
 * @todo
 * S.Cupp: The tabulated error codes should be reviewed by Leo Werneck for validity
 * and to check if any should be split or merged with other error codes. E.g. does
 * ghl_error_table_bisection actually only represent a single error, or should it
 * be divided into two?
 */
typedef enum {
  ghl_success,                   /**< Success/no error found */
  ghl_error_unknown_eos_type,   /**< Found an invalid @ref ghl_eos_t key */
  ghl_error_invalid_c2p_key,    /**< Found an invalid @ref ghl_con2prim_method_t key */
  ghl_error_neg_rho,            /**< Negative density from @ref Con2Prim routine */
  ghl_error_neg_pressure,       /**< Negative pressure from @ref Con2Prim routine */
  ghl_error_neg_vsq,            /**< Imaginary velocity from @ref Con2Prim routine */
  ghl_error_c2p_max_iter,       /**< Maximum iterations reached in @ref Con2Prim routine */
  ghl_error_c2p_singular,       /**< Singular value found in @ref Con2Prim routine */
  ghl_error_root_not_bracketed, /**< Value not contained by bounds in @ref Con2Prim root finder routine */
  ghl_error_table_max_rho,      /**< Density above table bounds in @ref tab_eos routine */
  ghl_error_table_min_rho,      /**< Density below table bounds in @ref tab_eos routine */
  ghl_error_table_max_ye,       /**< Electron fraction above table bounds in @ref tab_eos routine */
  ghl_error_table_min_ye,       /**< Electron fraction below table bounds in @ref tab_eos routine */
  ghl_error_table_max_T,        /**< Temperature above table bounds in @ref tab_eos routine */
  ghl_error_table_min_T,        /**< Temperature below table bounds in @ref tab_eos routine */
  ghl_error_exceed_table_vars,  /**< Requested more output variables than exist in the @ref tab_eos table */
  ghl_error_table_neg_energy,   /**< Negative energy found after energy shift in @ref tab_eos routine */
  ghl_error_table_bisection,    /**< Failure to find solution via bisection in @ref tab_eos routine */
  ghl_error_u0_singular         /**< Singular \f$ u^0 \f$ while computing velocities */
} ghl_error_codes_t;

/**
 * @ingroup Con2Prim
 * @enum ghl_con2prim_method_t
 * @brief Integer constants to specify conservative-to-primitive method
 */
typedef enum {
  None = -1,            /**< No method (for disabling backups) */
  Noble2D,              /**< The @ref ghl_hybrid_Noble2D routine */
  Noble1D,              /**< The @ref ghl_hybrid_Noble1D routine */
  Noble1D_entropy,      /**< The @ref ghl_hybrid_Noble1D routine */
  Font1D,               /**< The @ref ghl_hybrid_Font1D routine */
  Palenzuela1D,         /**< The @ref ghl_hybrid_Palenzuela1D_energy or @ref ghl_tabulated_Palenzuela1D_energy routine */
  Palenzuela1D_entropy, /**< The @ref ghl_hybrid_Palenzuela1D_entropy or @ref ghl_tabulated_Palenzuela1D_entropy routine */
  Newman1D,             /**< The @ref ghl_tabulated_Newman1D_energy routine */
  Newman1D_entropy      /**< The @ref ghl_tabulated_Newman1D_entropy routine */
} ghl_con2prim_method_t;

/**
 * @ingroup EOS
 * @enum ghl_eos_t
 * @brief Integer constants to track EOS type
 */
typedef enum {
  ghl_eos_simple,   /**< Simple EOS */ 
  ghl_eos_hybrid,   /**< Hybrid EOS */ 
  ghl_eos_tabulated /**< Tabulated EOS */ 
} ghl_eos_t;

/**
 * @ingroup GRHayL_Core
 * @struct ghl_parameters
 * @brief Stores basic GRHayL parameters
 */
typedef struct ghl_parameters {
  /** Primary @ref Con2Prim method */
  ghl_con2prim_method_t main_routine;
  /** Secondary @ref Con2Prim methods */
  ghl_con2prim_method_t backup_routine[3];
  /** Whether entropy is evolved (true) or not (false) */
  bool evolve_entropy;
  /** Whether temperature is evolved (true) or not (false) */
  bool evolve_temp;
  /** Maximum Lorentz factor for velocity limiters */
  double max_Lorentz_factor;
  /** Pre-computed value of the inverse square of max_Lorenz_factor */
  double inv_sq_max_Lorentz_factor;
  /** Threshold of \f$ \psi^6 \f$ above which limits on primitives and
      conservatives change. Approximates whether the point is inside the horizon. */
  double psi6threshold;
  /** Damping factor for the generalized Lorenz gauge in the @ref Induction */
  double Lorenz_damping_factor;

  // Con2Prim parameters
  /** Whether @ref ghl_con2prim_multi_method should calculate an initial
      primitive variable guess (true) or use user input (false) */
  bool calc_prim_guess;
  /** Maximum number of iterations allowed in @ref Con2Prim solvers */
  int con2prim_max_iterations;
  /** Convergence tolerance in @ref Con2Prim solvers */
  double con2prim_solver_tolerance;

  // PPM parameters
  /** Parameter controlling the flattening algorithm in @ref ghl_shock_detection_ftilde */
  double ppm_flattening_epsilon;
  /** Parameter controlling the flattening algorithm in @ref ghl_shock_detection_ftilde */
  double ppm_flattening_omega1;
  /** Parameter controlling the flattening algorithm in @ref ghl_shock_detection_ftilde */
  double ppm_flattening_omega2;
  /** Parameter controlling the steepening algorithm in @ref ghl_steepen_var */
  double ppm_shock_k0;
  /** Parameter controlling the steepening algorithm in @ref ghl_steepen_var */
  double ppm_shock_eta1;
  /** Parameter controlling the steepening algorithm in @ref ghl_steepen_var */
  double ppm_shock_eta2;
  /** Parameter controlling the steepening algorithm in @ref ghl_steepen_var */
  double ppm_shock_epsilon;
} ghl_parameters;

/**
 * @ingroup GRHayL_Core
 * @struct ghl_primitive_quantities
 * @brief Stores pointwise data for hydrodynamic primitive variables
 */
typedef struct ghl_primitive_quantities {
  /** Baryonic density \f$ \rho \f$ of the fluid */
  double rho;
  /** Pressure \f$ P \f$ of the fluid */
  double press;
  /** Specific internal energy \f$ \epsilon \f$ of the fluid */
  double eps;
  /** 0th component of fluid 4-velocity \f$ u^0 \f$ */
  double u0;
  /** Fluid 3-velocity \f$ \frac{u^i}{u^0} \f$ */
  double vU[3];
  /** Magnetic field \f$ B^i \f$ */
  double BU[3];
  /** Entropy \f$ S \f$ */
  double entropy;
  /** Electron fraction \f$ Y_e \f$ */
  double Y_e;
  /** Temperature \f$ T \f$ */
  double temperature;
} ghl_primitive_quantities;

/**
 * @ingroup GRHayL_Core
 * @struct ghl_conservative_quantities
 * @brief Stores pointwise data for hydrodynamic conservative variables
 *
 * @details
 * This struct can be used to store both the densitized and
 * undensitized versions of these variables. Densitized is the standard
 * expectation. Functions expecting or returning undensitized variants
 * will explicitly note this in the documentation and variable name.
 */
typedef struct ghl_conservative_quantities {
  /** Conservative density variable \f$ \tilde{D} \f$ (also called \f$ \rho_* \f$) */
  double rho;
  /** Conservative energy variable \f$ \tilde{\tau} \f$ */
  double tau;
  /** Conservative electron fraction variable \f$ \tilde{Y_e} \f$ */
  double Y_e;
  /** Conservative momentum variable \f$ \tilde{S}_i \f$ */
  double SD[3];
  /** Conservative entropy \f$ \tilde{S} \f$ */
  double entropy;
} ghl_conservative_quantities;

/**
 * @ingroup GRHayL_Core
 * @struct ghl_metric_quantities
 * @brief Stores pointwise data for spacetime quantities
 *
 * @details
 * This metric can be used to store any 3-metric, but
 * most @grhayl functions expect the ADM metric. Any exceptions
 * (such as taking the BSSN metric as an input) will be
 * explicitly noted in the documentation. The auxiliary quantities
 * can be automatically computed by @grhayl, so we recommend using
 * @ref ghl_initialize_metric to initialize this struct.
 */
typedef struct ghl_metric_quantities {
  /** The lapse \f$ \alpha \f$ */
  double lapse;
  /** The inverse lapse \f$ \frac{1}{\alpha} \f$ */
  double lapseinv;
  /** The inverse squared lapse \f$ \frac{1}{\alpha^2} \f$ */
  double lapseinv2;
  /** The determinant of the metric \f$ \left| \gamma \right| \f$ */
  double detgamma;
  /** The square root of the determinant of the metric \f$ \sqrt{\left| \gamma \right|} \f$ */
  double sqrt_detgamma;
  /** The shift \f$ \beta^i \f$ */
  double betaU[3];
  /** The 3-metric \f$ \gamma_{ij} \f$ */
  double gammaDD[3][3];
  /** The inverse 3-metric \f$ \gamma^{ij} \f$ */
  double gammaUU[3][3];
} ghl_metric_quantities;

/**
 * @ingroup GRHayL_Core
 * @struct ghl_ADM_aux_quantities
 * @brief Stores auxiliary data derived from ADM metric quantities
 *
 * @details
 * This struct is usually used in concert with an instance of ghl_metric_quantities.
 * As such, using @ref ghl_compute_ADM_auxiliaries is highly recommended.
 */
typedef struct ghl_ADM_aux_quantities {
  /** 4-metric \f$ g_{\mu\nu} \f$ */
  double g4DD[4][4];
  /** Inverse 4-metric \f$ g^{\mu\nu} \f$ */
  double g4UU[4][4];
} ghl_ADM_aux_quantities;

/**
 * @ingroup GRHayL_Core
 * @struct ghl_extrinsic_curvature
 * @brief Stores pointwise data for extrinsic curvature
 *
 * @details
 * This array is used for both the raised and lowered variant
 * of this tensor, and individual functions specify their output or
 * expected input.
 */
typedef struct ghl_extrinsic_curvature {
  /** Extrinsic curvature array \f$ K^{ij} \f$ or \f$ K_{ij} \f$ */
  double K[3][3];
} ghl_extrinsic_curvature;

/**
 * @ingroup GRHayL_Core
 * @struct ghl_stress_energy
 * @brief Stores pointwise data for stress-energy tensor
 *
 * @details
 * This array is used for both the raised and lowered variant
 * of this tensor, and individual functions specify their output or
 * expected input.
 */
typedef struct ghl_stress_energy {
  /** Stress-energy tensor array \f$ T^{\mu\nu} \f$ or \f$ T_{\mu\nu} \f$ */
  double T4[4][4];
} ghl_stress_energy;

/**
 * @ingroup EOS
 * @struct ghl_eos_parameters
 * @brief Stores parameters controlling the @ref EOS
 *
 * @details
 *
 * @todo
 * Tabulated variables need descriptions
 */
typedef struct ghl_eos_parameters {
  /** @name General Parameters */
  /** @ref ghl_eos_t setting the EOS type */
  ghl_eos_t eos_type;
  /** Atmospheric value of the density \f$ \rho \f$ */
  double rho_atm;
  /** Minimum allowed density \f$ \rho \f$ */
  double rho_min;
  /** Maximum allowed density \f$ \rho \f$ */
  double rho_max;
  /** Atmospheric value of the conservative energy variable \f$ \tau \f$ */
  double tau_atm;
  /** Atmospheric value of the eps \f$ \epsilon \f$ */
  double eps_atm;
  /** Minimum allowed eps \f$ \epsilon \f$ */
  double eps_min;
  /** Maximum allowed eps \f$ \epsilon \f$ */
  double eps_max;
  /** Atmospheric value of the pressure \f$ P \f$ */
  double press_atm;
  /** Minimum allowed pressure \f$ P \f$ */
  double press_min;
  /** Maximum allowed pressure \f$ P \f$ */
  double press_max;
  /** Atmospheric value of the entropy \f$ S \f$ */
  double entropy_atm;
  /** Minimum allowed entropy \f$ S \f$ */
  double entropy_min;
  /** Maximum allowed entropy \f$ S \f$ */
  double entropy_max;

  /** @name Hybrid EOS Parameters */
  /** Maximum allowed number of polytropic pieces */
#define MAX_EOS_PARAMS (10)
  /** Number of polytropic pieces (must be \f$ \eq \f$ @ref MAX_EOS_PARAMS) */
  int neos;
  /** Density boundaries for each polytropic piece */
  double rho_ppoly[MAX_EOS_PARAMS - 1];
  /** Polytropic indices \f$ \Gamma \f$ for each polytropic piece */
  double Gamma_ppoly[MAX_EOS_PARAMS];
  /** Polytropic indices \f$ \Gamma \f$ for each polytropic piece */
  double K_ppoly[MAX_EOS_PARAMS];
  /** Integration constants for specific internal energy \f$ \epsilon \f$ for each polytropic piece */
  double eps_integ_const[MAX_EOS_PARAMS];
  /** Adiabatic index\f$ \Gamma_th \f$ defining the thermal contribution to the EOS */
  double Gamma_th;

  /** @name Tabulated EOS Parameters */
  /** Atmospheric value of the electron fraction \f$ Y_e \f$ */
  double Y_e_atm;
  /** Minimum allowed electron fraction \f$ Y_e \f$ */
  double Y_e_min;
  /** Maximum allowed electron fraction \f$ Y_e \f$ */
  double Y_e_max;
  /** Atmospheric value of the temperature \f$ T \f$ */
  double T_atm;
  /** Minimum allowed temperature \f$ T \f$ */
  double T_min;
  /** Maximum allowed temperature \f$ T \f$ */
  double T_max;
  /** root_finding_precision: root-finding precision for table inversions */
  double root_finding_precision;

  /** @name Internal Tabulated EOS Variables */
  /** Number of density points in the EOS table */
  int N_rho;
  /** Number of temperature points in the EOS table */
  int N_T;
  /** Number of electron fraction points in the EOS table */
  int N_Ye;

  /** Pointer to full EOS data table */
  double *restrict table_all;
  /** @todo needs documentation */
  double *restrict table_logrho;
  /** @todo needs documentation */
  double *restrict table_logT;
  /** @todo needs documentation */
  double *restrict table_Y_e;
  /** @todo needs documentation */
  double *restrict table_eps;

  /** Minimum density in the table */
  double table_rho_min;
  /** Maximum density in the table */
  double table_rho_max;
  /** Minimum temperature in the table */
  double table_T_min;
  /** Maximum temperature in the table */
  double table_T_max;
  /** Minimum electron fraction in the table */
  double table_Y_e_min;
  /** Maximum electron fraction in the table */
  double table_Y_e_max;
  /** Minimum pressure in the table */
  double table_P_min;
  /** Maximum pressure in the table */
  double table_P_max;
  /** Minimum specific internal energy in the table */
  double table_eps_min;
  /** Maximum specific internal energy in the table */
  double table_eps_max;
  /** Minimum entropy in the table */
  double table_ent_min;
  /** Maximum entropy in the table */
  double table_ent_max;

  /** The energy shift \f$ \epsilon_0 \f$ such that \f$ \epsilon + \epsilon_0 \f$ is always positive */
  double energy_shift;
  /** Inverse of the temperature difference \f$ \Delta T = \log{T_1} - \log{T_0} \f$ for the table */
  double dtempi;
  /** Inverse of the density difference \f$ \Delta \rho = \log{\rho_1} - \log{\rho_0} \f$ for the table */
  double drhoi;
  /** Inverse of the electron fraction difference \f$ \Delta Y_e = {Y_e}_1 - {Y_e}_0 \f$ for the table */
  double dyei;
  /** Value of ghl_eos_parameters::drhoi * ghl_eos_parameters::dtempi */
  double drhotempi;
  /** Value of ghl_eos_parameters::drhoi * ghl_eos_parameters::dyei */
  double drhoyei;
  /** Value of ghl_eos_parameters::dtempi * ghl_eos_parameters::dyei */
  double dtempyei;
  /** Value of ghl_eos_parameters::drhoi * ghl_eos_parameters::dtempi * ghl_eos_parameters::dyei */
  double drhotempyei;

  // These are used for beta-equilibrium
  /** @todo needs documentation */
  double *lp_of_lr;
  /** @todo needs documentation */
  double *le_of_lr;
  /** @todo needs documentation */
  double *Ye_of_lr;

} ghl_eos_parameters;

#ifdef __cplusplus
extern "C" {
#endif

char *ghl_get_con2prim_routine_name(const ghl_con2prim_method_t key);

void ghl_initialize_eos_functions(
    const ghl_eos_t eos_type);

void ghl_initialize_simple_eos(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *restrict eos);

void ghl_initialize_hybrid_eos(
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *restrict eos);

void ghl_initialize_tabulated_eos(
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

void ghl_initialize_simple_eos_functions_and_params(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *restrict eos);

void ghl_initialize_hybrid_eos_functions_and_params(
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *restrict eos);

void ghl_initialize_tabulated_eos_functions_and_params(
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
      const ghl_con2prim_method_t main_routine,
      const ghl_con2prim_method_t backup_routine[3],
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
      const ghl_metric_quantities *restrict ADM_metric,
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
      ghl_metric_quantities *restrict ADM_metric);

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
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims,
      bool *restrict speed_limited);

void ghl_compute_TDNmunu(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_stress_energy *restrict Tmunu);

void ghl_compute_TUPmunu(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_stress_energy *restrict Tmunu);

void ghl_compute_smallb_and_b2(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_primitive_quantities *restrict prims,
      const double uDN[4],
      double smallb[4],
      double *restrict smallb2);

void ghl_read_error_codes(
      const ghl_error_codes_t error);

#ifdef __cplusplus
}
#endif

#include "ghl_eos_functions.h"
#include "ghl_debug.h"

#endif // GHL_H
