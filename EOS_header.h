// Thorn      : IllinoisGRMHD
// File       : EOS_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains function prototypes for the
//              EOS functions available in IllinoisGRMHD

#ifndef EOS_HEADER_H
#define EOS_HEADER_H

#include "GRMHD_header.h"

//------------------- EOS struct -----------------

// This struct defines all the necessary quantities
// that IllinoisGRMHD needs in order to use the
// EOS functions correctly

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

 --rho_ppoly_tab: array of the density values which divide the polytropic pieces

 --Gamma_ppoly_tab: array of the polytropic indices

 --K_ppoly_tab: array of the adiabatic constants

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
  double rho_ppoly_tab[MAX_EOS_PARAMS-1];
  double Gamma_ppoly_tab[MAX_EOS_PARAMS];
  double K_ppoly_tab[MAX_EOS_PARAMS];
  double eps_integ_const[MAX_EOS_PARAMS];
  double Gamma_th;
  //------------------------------------------------

  //---------- Tabulated Equation of State ---------
  double Ye_atm, Ye_min, Ye_max;
  double temp_atm, temp_min, temp_max;
  double eps_atm, eps_min, eps_max;
  double entropy_atm, entropy_min, entropy_max;
  double root_finding_precision;
  double depsdT_threshold;
  //------------------------------------------------
} eos_parameters;

//---------- Initialization routines ---------------
void initialize_general_eos(eos_parameters *restrict eos, const int type,
             const double tau_atm, const double W_max, const double inv_W_max_squared,
             const double eps_atm, const double eps_min, const double eps_max,
             const double press_atm, const double press_min, const double press_max,
             const double entropy_atm, const double entropy_min, const double entropy_max,
             const double rho_atm, const double rho_min, const double rho_max);

void initialize_hybrid_eos(eos_parameters *restrict eos, const int neos, const double rho_ppoly_tab[],
             const double Gamma_ppoly_tab[], const double K_ppoly_tab[],
             const double eps_integ_const[], const double Gamma_th);

void initialize_tabulated_eos(eos_parameters *restrict eos, const double precision, const double threshold,
             const double temp_atm, const double temp_min, const double temp_max,
             const double Ye_atm, const double Ye_min, const double Ye_max);

//------------------------------------------------

void compute_P_cold__eps_cold(const eos_parameters *restrict eos, const double rho_in,
                              double *restrict P_cold_ptr, double *restrict eps_cold_ptr);

#endif // EOS_HEADER_H


