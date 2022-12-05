#ifndef __HARM_U2P_UTIL__C__
#define __HARM_U2P_UTIL__C__
/*
  -------------------------------------------------------------------------------
  Copyright 2005 Scott C. Noble, Charles F. Gammie,
  Jonathan C. McKinney, and Luca Del Zanna


  This file is part of PVS-GRMHD.

  PVS-GRMHD is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  PVS-GRMHD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with PVS-GRMHD; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------------------
*/

#include "GRMHD_header.h"
#include <stdio.h>
static bool print;

// HARM uses lots of globals. These auxiliary variables
// allow us to pass useful quantities to the con2prim
// functions without having a large number of function
// arguments.
typedef struct _harm_auxiliary_vars_ {
  double Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,gamma,gamma_times_S,ye,T_guess;
  bool use_entropy;
} harm_aux_vars_struct;

/**********************************************************************
    raise_g():

         -- calculates the contravariant form of a covariant tensor,
            using the inverse of the metric;
***********************************************************************/
static inline void raise_g(const double vcov[NDIM], const double gcon[NDIM][NDIM], double vcon[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcon[i] = 0. ;
    for(j=0;j<NDIM;j++)
      vcon[i] += gcon[i][j]*vcov[j];
  }

  return ;
}


/**********************************************************************
     lower_g():

          -- calculates the ocvariant form of a contravariant tensor
             using the metric;
***********************************************************************/
static inline void lower_g(const double vcon[NDIM], const double gcov[NDIM][NDIM], double vcov[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcov[i] = 0. ;
    for(j=0;j<NDIM;j++)
      vcov[i] += gcov[i][j]*vcon[j] ;
  }

  return ;
}


/**********************************************************************
     ncov_calc():

         -- calculates the covariant form of the normal vector to our
            spacelike hypersurfaces ala the ADM formalism.

         -- requires the inverse metric;
***********************************************************************/
static inline void ncov_calc(const double gcon[NDIM][NDIM],double ncov[NDIM])
{
  double lapse ;
  int i;

  lapse = sqrt(-1./gcon[0][0]) ;

  ncov[0] = -lapse ;
  for( i = 1; i < NDIM; i++) {
    ncov[i] = 0. ;
  }

  return ;
}

/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/*
pressure as a function of rho0 and u
this is used by primtoU and Utoprim_?D
*/
static inline double pressure_rho0_u(const eos_parameters *restrict eos, const double rho0, const double u)
{

  // Compute P_cold, eps_cold
  double P_cold, eps_cold;
  compute_P_cold_and_eps_cold(eos, rho0, &P_cold, &eps_cold);

  /* Compute the pressure as a function of rho_b (rho0) and
   * u = rho_b * eps, using our hybrid EOS:
   * .-------------------------------------------------------------.
   * | p(rho_b,u) = P_cold + (Gamma_th - 1)*(u - rho_b * eps_cold) |
   * .-------------------------------------------------------------.
   */
//  return( P_cold + (eos->Gamma_th - 1.0)*(u - rho0*eps_cold) );
  return((eos->Gamma_th - 1.0)*u) ;

}


/*
   pressure as a function of rho0 and w = rho0 + u + p
   this is used by primtoU and Utoprim_1D
*/
static inline double pressure_rho0_w(const eos_parameters *restrict eos, const double rho0, const double w)
{

  // Compute P_cold, eps_cold
  double P_cold, eps_cold;
  compute_P_cold_and_eps_cold(eos,rho0, &P_cold, &eps_cold);

  /* Compute the pressure as a function of rho_b (rho0) and
   * w = u + rho_b + p, using our hybrid EOS:
   *  ----------------------------------------------------------------------------
   * | p(rho_b,w) = ( P_cold + (Gamma_th-1)*( w - rho_b*(1+eps_cold) ) )/Gamma_th |
   *  ----------------------------------------------------------------------------
   */
//  return( (P_cold + (eos->Gamma_th-1.0)*( w - rho0*(1.0 + eps_cold) ) )/eos->Gamma_th );
  return((eos->Gamma_th-1.0)*(w - rho0)/eos->Gamma_th ) ;
}

#endif
