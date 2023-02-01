#ifndef FLUX_SOURCE_H_
#define FLUX_SOURCE_H_

#include "GRHayL.h"

/*
   The struct reconstructed_prims_struct contains variables for storing the (point-wise)
   reconstructed primitive variables. The struct elements are detailed below:
   
 --rhob: the baryonic density
 
 --P: the pressure
 
 --u4U*: the contravariant fluid 4-velocity
   
 --BU*: the contravariant magnetic field
 
 --h: the enthalpy
 
 --Gamma_th: an EOS quantity
 
 --epsilon_th: an EOS quantity
 
 --dPcold_drhob: an EOS quantity
*/

typedef struct reconstructed_prims_struct {
  double u4U0, u4U1, u4U2, u4U3;
  double BU0, BU1, BU2;
  double rhob;
  double P, h, Gamma_th, epsilon_th, dPcold_drhob;
} reconstructed_prims_struct;

/*
   The struct conservative_fluxes_struct contains variables for HLLE fluxes that the user
   will finite difference later. Note that we do not store the fluxes for all three
   directions. The struct elements are detailed below:
   
 --cmin* / cmax*: the minimum and maximum characteristic speeds for all three flux directions
 
 --HLLE_flux_StildeD*: the HLLE flux for the momentum term term, all three components
 
 --HLLE_flux_rho_star: the HLLE flux for the baryonic density
   
 --HLLE_flux_tau_tilde: the HLLE flux for the energy
*/
typedef struct conservative_fluxes_struct {
  double cmin_dirn0, cmin_dirn1, cmin_dirn2;
  double cmax_dirn0, cmax_dirn1, cmax_dirn2;
  double HLLE_flux_StildeD0, HLLE_flux_StildeD1, HLLE_flux_StildeD2;
  double HLLE_flux_rho_star, HLLE_flux_tau_tilde;
} conservative_fluxes_struct;

/*
   The struct conservative_sources_struct contains stores the point-wise values 
   of the source terms for the momentum and energy equations. The struct 
   elements are detailed below:
   
 --StildeD*_src: momentum equations source terms
 
 --tau_tilde_src: tau equation source term
*/
typedef struct conservative_sources_struct {
  double StildeD0_src, StildeD1_src, StildeD2_src;
  double tau_tilde_src;
} conservative_sources_struct;

/*
   The struct metric_quantities_struct contains variables for storing the (point-wise)
   metric variables. The struct elements are detailed below:
   
 --alpha: the lapse
 
 --betaU*: the shift vector
   
 --gammaDD**: the physical 3-metric
 
 --KDD*: the extrinsic curvature
*/
typedef struct metric_quantities2 {
  double alpha;
  double betaU0, betaU1, betaU2;
  double gammaDD01, gammaDD02, gammaDD12, gammaDD00, gammaDD11, gammaDD22;
  double KDD01, KDD02, KDD12, KDD00, KDD11, KDD22;
} metric_quantities2;

/*
   The struct metric_quantities_derivatives_struct contains variables for
   precomputed derivatives of the metric quantities.
*/

typedef struct metric_quantities_derivatives {
  double alpha_dD1;
  double gammaDD_dD111;
  double alpha_dD2;
  double betaU_dD12;
  double gammaDD_dD020;
  double gammaDD_dD122;
  double gammaDD_dD001;
  double gammaDD_dD022;
  double gammaDD_dD121;
  double alpha_dD0;
  double betaU_dD10;
  double betaU_dD01;
  double gammaDD_dD110;
  double betaU_dD00;
  double betaU_dD22;
  double gammaDD_dD221;
  double gammaDD_dD112;
  double gammaDD_dD021;
  double betaU_dD21;
  double betaU_dD20;
  double gammaDD_dD012;
  double betaU_dD02;
  double gammaDD_dD222;
  double gammaDD_dD220;
  double gammaDD_dD010;
  double betaU_dD11;
  double gammaDD_dD002;
  double gammaDD_dD000;
  double gammaDD_dD120;
  double gammaDD_dD011;
} metric_quantities_derivatives;

#endif //FLUX_SOURCE_H_
