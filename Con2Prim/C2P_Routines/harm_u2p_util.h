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

#include "con2prim.h"

/************* This version of vars are from Leo ****************/
//static const int MAX_NEWT_ITER       = 50;     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
//static const double NEWT_TOL      = 5e-9;    /* Min. of tolerance allowed for Newton-Raphson iterations */
//static const double MIN_NEWT_TOL  = 5e-9;    /* Max. of tolerance allowed for Newton-Raphson iterations */
/****************************************************************/
static const int NPR =8;
static const int NDIM=4;

static const int MAX_NEWT_ITER    = 30;     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
//#define MAX_NEWT_ITER 300     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
static const double NEWT_TOL      = 1.0e-10;    /* Min. of tolerance allowed for Newton-Raphson iterations */
static const double MIN_NEWT_TOL  = 1.0e-10;    /* Max. of tolerance allowed for Newton-Raphson iterations */
static const int EXTRA_NEWT_ITER  = 0; /* ZACH SAYS: Original value = 2. But I don't think this parameter > 0 is warranted. Just slows the code for no reason, since our tolerances are fine. */

static const double NEWT_TOL2     = 1.0e-15;      /* TOL of new 1D^*_{v^2} gnr2 method */
static const double MIN_NEWT_TOL2 = 1.0e-10;  /* TOL of new 1D^*_{v^2} gnr2 method */

static const double Z_TOO_BIG    = 1.e20;    /* \gamma^2 (\rho_0 + u + p) is assumed
                                                  to always be smaller than this.  This
                                                  is used to detect solver failures */
static const double UTSQ_TOO_BIG = 1.e20;    /* \tilde{u}^2 is assumed to be smaller
                                                  than this.  Used to detect solver
                                                  failures */

static const double FAIL_VAL     = 1.e30;    /* Generic value to which we set variables when a problem arises */

static const double NUMEPSILON   = 2.2204460492503131e-16;

// HARM uses lots of globals. These auxiliary variables
// allow us to pass useful quantities to the con2prim
// functions without having a large number of function
// arguments.
typedef struct _harm_auxiliary_vars_ {
  double Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,W,W_times_S,ye,T_guess;
  bool use_entropy;
} harm_aux_vars_struct;

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
  eos->hybrid_compute_P_cold_and_eps_cold(eos, rho0, &P_cold, &eps_cold);

  /* Compute the pressure as a function of rho_b (rho0) and
   * u = rho_b * eps, using our hybrid EOS:
   * .-------------------------------------------------------------.
   * | p(rho_b,u) = P_cold + (Gamma_th - 1)*(u - rho_b * eps_cold) |
   * .-------------------------------------------------------------.
   */
  return( P_cold + (eos->Gamma_th - 1.0)*(u - rho0*eps_cold) );
}


/*
   pressure as a function of rho0 and w = rho0 + u + p
   this is used by primtoU and Utoprim_1D
*/
static inline double pressure_rho0_w(const eos_parameters *restrict eos, const double rho0, const double w)
{

  // Compute P_cold, eps_cold
  double P_cold, eps_cold;
  eos->hybrid_compute_P_cold_and_eps_cold(eos,rho0, &P_cold, &eps_cold);

  /* Compute the pressure as a function of rho_b (rho0) and
   * w = u + rho_b + p, using our hybrid EOS:
   *  ----------------------------------------------------------------------------
   * | p(rho_b,w) = ( P_cold + (Gamma_th-1)*( w - rho_b*(1+eps_cold) ) )/Gamma_th |
   *  ----------------------------------------------------------------------------
   */
  return( (P_cold + (eos->Gamma_th-1.0)*( w - rho0*(1.0 + eps_cold) ) )/eos->Gamma_th );
}

/****************************************************************************
   vsq_calc():

      -- evaluate v^2 (spatial, normalized velocity) from
            W = \gamma^2 w

****************************************************************************/
static inline double vsq_calc(const harm_aux_vars_struct *restrict harm_aux, const double W)
{
  const double Wsq = W*W ;
  const double Xsq = (harm_aux->Bsq + W) * (harm_aux->Bsq + W);

  return(  ( Wsq * harm_aux->Qtsq  + harm_aux->QdotBsq * (harm_aux->Bsq + 2.*W)) / (Wsq*Xsq) );
}

/**********************************************************************/
/**********************************************************************
  pressure_W_vsq():

        -- Hybrid single and piecewise polytropic equation of state;
        -- pressure as a function of P_cold, eps_cold, W, vsq, and D:
**********************************************************************/
static inline double pressure_W_vsq(const eos_parameters *restrict eos, const double W, const double vsq, const double D)
{

  // Compute gamma^{-2} = 1 - v^{2} and gamma^{-1}
  const double inv_gammasq = 1.0 - vsq;
  const double inv_gamma   = sqrt(inv_gammasq);

  // Compute rho_b = D / gamma
  const double rho_b = D*inv_gamma;

  // Compute P_cold and eps_cold
  double P_cold, eps_cold;
  eos->hybrid_compute_P_cold_and_eps_cold(eos,rho_b, &P_cold, &eps_cold);

  // Compute p = P_{cold} + P_{th}
  return( ( P_cold + (eos->Gamma_th - 1.0)*( W*inv_gammasq - D*inv_gamma*( 1.0 + eps_cold ) ) )/eos->Gamma_th );
}

/**********************************************************************/
/**********************************************************************
  dpdvsq_calc():

      -- partial derivative of pressure with respect to vsq
**********************************************************************/
static inline double dpdvsq_calc(const eos_parameters *restrict eos, const double W, const double vsq, const double D)
{

  // Set gamma and rho
  const double gamma = 1.0/sqrt(1.0 - vsq);
  const double rho_b = D/gamma;

  // Compute P_cold and eps_cold
  double P_cold, eps_cold;
  eos->hybrid_compute_P_cold_and_eps_cold(eos, rho_b, &P_cold, &eps_cold);

  // Set basic polytropic quantities
  const int polytropic_index = eos->hybrid_find_polytropic_index(eos, rho_b);
  const double Gamma_ppoly = eos->Gamma_ppoly[polytropic_index];


  /* Now we implement the derivative of P_cold with respect
   * to v^{2}, given by
   *  ----------------------------------------------------
   * | dP_cold/dvsq = gamma^{2 + Gamma_{poly}/2} P_{cold} |
   *  ----------------------------------------------------
   */
  const double dPcold_dvsq = P_cold * pow(gamma,2.0 + 0.5*Gamma_ppoly);


  /* Now we implement the derivative of eps_cold with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------------
   * | deps_cold/dvsq = gamma/(D*(Gamma_ppoly-1)) * (dP_cold/dvsq + gamma^{2} P_cold / 2) |
   *  -----------------------------------------------------------------------------------
   */
  const double depscold_dvsq = ( gamma/(D*(Gamma_ppoly-1.0)) ) * ( dPcold_dvsq + 0.5*gamma*gamma*P_cold );

  /* Now we implement the derivative of p_hybrid with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------
   * | dp/dvsq = Gamma_th^{-1}( dP_cold/dvsq                                       |
   * |                          + (Gamma_{th}-1)*(-W                               |
   * |                                            + D gamma (1 + eps_cold)/2       |
   * |                                            - (D/gamma) * deps_cold/dvsq) )  |
   *  -----------------------------------------------------------------------------
   */
  return( ( dPcold_dvsq + (eos->Gamma_th-1.0)*( -W + D*gamma*(1+eps_cold)/2.0 - D*depscold_dvsq/gamma ) )/eos->Gamma_th );
}

int general_newton_raphson(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double indep_var_in,
      int *restrict n_iter_ptr,
      double x[],
      void (*funcd)(const eos_parameters *restrict, const harm_aux_vars_struct *restrict, const int, const double, const double [], double [],
                    double [], double [][ndim], double *restrict, double *restrict, int *restrict));

int newton_raphson_1d(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double indep_var_in,
      int *restrict n_iter_ptr,
      double x[],
      void (*funcd)(const eos_parameters *restrict, const harm_aux_vars_struct *restrict, const int, const double, const double [], double [],
                    double [], double [][ndim], double *restrict, double *restrict, int *restrict));

void func_vsq(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double dummy,
      const double x[],
      double dx[],
      double resid[],
      double jac[][2],
      double *restrict f,
      double *restrict df,
      int *restrict n_iter);

void func_1d_orig(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double dummy,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df,
      int *restrict n_iter);

void func_W(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double rho_in,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df,
      int *restrict n_iter);

void func_rho(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double W_in,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df,
      int *restrict n_iter);

#endif
