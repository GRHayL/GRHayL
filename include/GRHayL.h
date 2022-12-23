#ifndef GRHayL_H_
#define GRHayL_H_

#include <math.h>
#include <stdbool.h>

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SQR(x) ((x) * (x))
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780

/*
   The struct GRHayL_parameters contains parameters for controlling
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

typedef struct GRHayL_parameters {
  int main_routine, backup_routine[3];
  bool evolve_entropy;
  bool evolve_temp;
  bool calc_prim_guess;
  double psi6threshold;
  bool update_Tmunu;
  bool Cupp_Fix;
  double Lorenz_damping_factor;
} GRHayL_parameters;

void initialize_GRHayL(
      const int main,
      const int backup[3],
      const int evolve_entropy,
      const int evolve_temp,
      const int calc_prim_guess,
      const double psi6threshold,
      const int update_Tmunu,
      const int Cupp_Fix,
      GRHayL_parameters *restrict params);

//--------------------------------------------------

//--------------------- EOS facets -----------------

// Maxmimum number of polytropic EOS pieces
#define MAX_EOS_PARAMS (10)

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

typedef struct eos_parameters {

  //-------------- General parameters --------------
  int eos_type;
  double rho_atm, rho_min, rho_max;
  double tau_atm;
  double press_atm, press_min, press_max;
  double W_max, inv_W_max_squared;
  //------------------------------------------------

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
  double Ye_atm, Ye_min, Ye_max;
  double T_atm, T_min, T_max;
  double eps_atm, eps_min, eps_max;
  double entropy_atm, entropy_min, entropy_max;
  double root_finding_precision;
  double depsdT_threshold;

  // Table size
  int N_rho, N_T, N_Ye;

  // Tabulated quantities
  double *restrict table_logrho;
  double *restrict table_logT;
  double *restrict table_Ye;
  double *restrict table_all;
  double *restrict table_eps;

  // Table bounds
  double table_rho_max, table_rho_min;
  double table_T_min  , table_T_max;
  double table_Ye_min , table_Ye_max;

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
  //------------------------------------------------

} eos_parameters;

#include "GRHayL_EOS_helpers.h"

void initialize_general_eos(
      const int type,
      const double tau_atm,
      const double W_max,
      const double entropy_atm,
      const double entropy_min,
      const double entropy_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      eos_parameters *restrict eos);

// Leo says: I think this function should
//           be called by the general one.
void initialize_hybrid_eos(
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly,
      const double Gamma_th,
      eos_parameters *restrict eos);

// Leo says: Same comment.
void initialize_tabulated_eos(
      const double precision,
      const double threshold,
      const double temp_atm,
      const double temp_min,
      const double temp_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      eos_parameters *restrict eos);

//--------------------------------------------------

//--------------- Con2Prim facets ------------------

/*
   The struct primitive_quantities contains variables for storing the (point-wise)
   primitive variable data. The struct elements are detailed below:

 --rho: the baryonic density rho_b

 --press: the pressure P

 --v*: the 3-velocity v^i used in IllinoisGRMHD. This is defined as u^i/u^0. The other
   commonly used choice is the Valencia 3-velocity Vv^i defined as
   Vv^i = u^i/W + beta^i/lapse.


 --B*: the magnetic field TODO: give specific B definition

 --entropy: the entropy S

 --Y_e: the electron fraction Y_e

 --temp: the temperature T

 --u: this is the energy variable needed by HARM-type routines. It is currently not
   possible to set this via the initialize_primitives function, as most code do not
   have this. A guess is automatically generated for this quantity from the other
   values.
*/

typedef struct primitive_quantities {
  double rho, press, eps;
  double vx, vy, vz;
  double Bx, By, Bz;
  double entropy, Y_e, temperature;
} primitive_quantities;

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
  double S_x, S_y, S_z;
  double entropy;
} conservative_quantities;

void initialize_primitives(
             const double rho, const double press, const double epsilon,
             const double vx, const double vy, const double vz,
             const double Bx, const double By, const double Bz,
             const double entropy, const double Y_e, const double temp,
             primitive_quantities *restrict prims);

void initialize_conservatives(
             const double rho, const double tau,
             const double S_x, const double S_y, const double S_z,
             const double Y_e, const double entropy,
             conservative_quantities *restrict cons);

void return_primitives(
             const primitive_quantities *restrict prims,
             double *restrict rho, double *restrict press, double *restrict epsilon,
             double *restrict vx, double *restrict vy, double *restrict vz,
             double *restrict Bx, double *restrict By, double *restrict Bz,
             double *restrict entropy, double *restrict Y_e, double *restrict temp);

void return_conservatives(
             const conservative_quantities *restrict cons,
             double *restrict rho, double *restrict tau,
             double *restrict S_x, double *restrict S_y, double *restrict S_z,
             double *restrict Y_e, double *restrict entropy);

//--------------------------------------------------

//-------------- Space-time facets -----------------

/* The struct metric_quantities contains variables for storing the (point-wise)
   metric data. The struct elements are detailed below:

 --bssn_phi: TODO

 --bssn_psi: TODO

 --bssn_gij: the BSSN conformal metric g_{i j}. These variables are set via the
   initialize_metric function.

 --betai: the shift vector (beta^x, beta^y, beta^z) from the 3+1
   decomposition. These variables are set via the initialize_metric function.

 --lapse: the lapse from the 3+1 decomposition. These variables are set via the
   initialize_metric function.

   All of the following varibles are calculated from the previous quantities in the initialize_metric
   function.

 --bssn_gupij: the BSSN conformal metric inverse g^{i j}. They are calculated from
   the bssn_gij variables using the relations
     gupxx =   ( gyy * gzz - gyz * gyz );
     gupxy = - ( gxy * gzz - gyz * gxz );
     gupxz =   ( gxy * gyz - gyy * gxz );
     gupyy =   ( gxx * gzz - gxz * gxz );
     gupyz = - ( gxx * gyz - gxy * gxz );
     gupzz =   ( gxx * gyy - gxy * gxy );

 --lapseinv: this variable stores the inverse of the lapse.

 --lapm1: this variable stores the quantity (lapse-1).

 --psi2: this variable stores the quantity exp(2.0*metric.bssn_phi).

 --psi4: this variable stores the quantity (psi2^2).

 --psi6: this variable stores the quantity (psi4*psi2).

 --psi4inv: this variable stores the inverse of psi4.

 --lapseinv2: this variables stores the quantity (lapse^(-2))

 --adm_gij: the ADM metric g_{i j} They are computed by the initialize_metric
   function using the relation adm_gij = psi4*bssn_gij.

 --adm_gupij: the ADM metric inverse g^{i j}. They are computed by the
   initialize_metric function using the relation adm_gupij = psi4inv*bssn_gupij.

 --g4dn: the 4-metric g_{\mu \nu}. This quantity is needed for computing T_{\mu \nu} and T^{\mu \nu}
   and the HARM con2prim lowlevel functions.

 --g4up: the 4-metric inverse g^{\mu \nu}. This quantity is needed for computing T_{\mu \nu} and T^{\mu \nu}
   and the HARM con2prim lowlevel functions.
*/

typedef struct metric_quantities {
  double adm_gxx, adm_gxy, adm_gxz;
  double adm_gyy, adm_gyz, adm_gzz;
  double adm_gupxx, adm_gupxy, adm_gupxz;
  double adm_gupyy, adm_gupyz, adm_gupzz;
  double betax, betay, betaz;
  double lapse, lapseinv;
  double psi2, psi4, psi6;
  double psi4inv, lapseinv2;
  double g4dn[4][4],g4up[4][4];
} metric_quantities;


//TODO: comment/add to this
typedef struct stress_energy {
  double Ttt, Ttx, Tty, Ttz;
  double Txx, Txy, Txz;
  double Tyy, Tyz, Tzz;
} stress_energy;

void initialize_metric(const double lapse,
             const double gxx, const double gxy, const double gxz,
             const double gyy, const double gyz, const double gzz,
             const double betax, const double betay, const double betaz,
             metric_quantities *restrict metric);

//--------------------------------------------------

#endif // GRHayL_H
