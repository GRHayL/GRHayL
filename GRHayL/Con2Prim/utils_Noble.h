#ifndef __HARM_U2P_UTIL__C__
#define __HARM_U2P_UTIL__C__

#include "con2prim.h"

static const int NPR =8;
static const int NDIM=4;

/*
   \gamma^2 (\rho_0 + u + p) is assumed
   to always be smaller than this.  This
   is used to detect solver failures
*/
static const double Z_TOO_BIG    = 1e20;

// HARM uses lots of globals. These auxiliary variables
// allow us to pass useful quantities to the con2prim
// functions without having a large number of function
// arguments.
typedef struct _harm_auxiliary_vars_ {
  double BbarU[3];
  double QU[4];
  double Bsq, QdotBsq, QdotB;
  double Qtsq, Qdotn;
  double D, W_times_S;
  int n_iter;
  int max_iterations;
  int solver_tolerance;
} harm_aux_vars_struct;

/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

static inline void ghl_compute_func_auxiliaries(
      const ghl_eos_parameters *restrict eos,
      const double Z,
      const double vsq,
      const double D,
      double *restrict pressure,
      double *restrict dPdvsq,
      double *restrict dPdZ) {

  // Compute W^{-2} = 1 - v^{2} and W^{-1}
  const double inv_Wsq = 1.0 - vsq;
  const double inv_W   = sqrt(inv_Wsq);
  const double W = 1.0/inv_W;

  // Compute rho_b = D / W
  const double rho_b = D*inv_W;

  // Compute P_cold and eps_cold
  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(eos,rho_b, &P_cold, &eps_cold);

  // Set up some intermediate quantities
  const double Gth_m1          = eos->Gamma_th - 1.0;
  const double eps_p1          = eps_cold + 1.0;
  const double Gth_m1_Z_invWsq = Gth_m1*Z*inv_Wsq;
  const double Gth_m1_rho_eps  = Gth_m1*rho_b*eps_p1;

  // Compute P = P_{cold} + P_{th}
  *pressure = (P_cold + Gth_m1_Z_invWsq - Gth_m1_rho_eps)/eos->Gamma_th;

  // Set basic polytropic quantities
  const double Gamma_ppoly = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos, rho_b)];

  // Note that for the following derivatives, we divide by W^2 to simplify the expressions

  /* Now we implement the derivative of P_cold with respect to v^{2}, given by
   *  -------------------------------------------------
   * | dP_cold/dvsq = W^{2 + Gamma_{poly}/2} P_{cold} |
   * | dP_cold/dvsq/W^2 = W^{Gamma_{poly}/2} P_{cold} |
   *  -------------------------------------------------
   */
  const double dPcold_dvsq_Wsq = P_cold * pow(W, 0.5*Gamma_ppoly);

  /* Now we implement the derivative of eps_cold with respect to v^{2}, given by
   *  -----------------------------------------------------------------------------
   * | deps_cold/dvsq = (dP_cold/dvsq + W^{2} P_cold/2) * W/(D*(Gamma_ppoly-1))   |
   * | deps_cold/dvsq/W^2 = (dP_cold/dvsq/W^2 + P_cold/2) / (W*D*(Gamma_ppoly-1)) |
   *  -----------------------------------------------------------------------------
   */
  const double depscold_dvsq_Wsq = ( dPcold_dvsq_Wsq + 0.5*P_cold ) / (rho_b*(Gamma_ppoly-1.0));

  /* Now we implement the derivative of p_hybrid with respect to v^{2}, given by
   *  --------------------------------------------------------------------------
   * | dp/dvsq = Gamma_th^{-1}( dP_cold/dvsq                                   |
   * |                          + (Gamma_{th}-1)*(-Z                           |
   * |                                            + D W (1 + eps_cold)/2       |
   * |                                            - (D/W) * deps_cold/dvsq) )  |
   * |         = W^2/Gamma_th * (dP_cold/dvsq/W^2                              |
   * |                         + (Gamma_{th}-1)*(-Z/W^2                        |
   * |                                           + rho (1 + eps_cold)/2        |
   * |                                           - rho * deps_cold/dvsq/W^2) ) |
   *  --------------------------------------------------------------------------
   */
  *dPdvsq = W*W*(dPcold_dvsq_Wsq - Gth_m1_Z_invWsq + Gth_m1_rho_eps/2.0 - Gth_m1*rho_b*depscold_dvsq_Wsq)/eos->Gamma_th;

  *dPdZ = Gth_m1*inv_Wsq/eos->Gamma_th;
}

/*
   pressure as a function of rho0 and w = rho0 + u + p
   this is used by primtoU and Utoprim_1D
*/
static inline double ghl_pressure_rho0_w(
      const ghl_eos_parameters *restrict eos,
      const double rho0,
      const double w) {

  // Compute P_cold, eps_cold
  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(eos,rho0, &P_cold, &eps_cold);

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
            Z = \gamma^2 w

****************************************************************************/
static inline double ghl_vsq_calc(
      const harm_aux_vars_struct *restrict harm_aux,
      const double Z) {
  const double Zsq = Z*Z ;
  const double Xsq = (harm_aux->Bsq + Z) * (harm_aux->Bsq + Z);

  return(  ( Zsq * harm_aux->Qtsq  + harm_aux->QdotBsq * (harm_aux->Bsq + 2.*Z)) / (Zsq*Xsq) );
}

int ghl_general_newton_raphson(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double indep_var_in,
      double x[],
      void (*validate_x)(
            const double [],
            double []),
      void (*funcd)(
            const ghl_parameters *restrict,
            const ghl_eos_parameters *restrict,
            harm_aux_vars_struct *restrict,
            const double,
            const double [],
            double [],
            double [],
            double [][ndim],
            double *restrict,
            double *restrict));

void ghl_func_Z(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double rho_in,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df);

void ghl_func_rho(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double Z_in,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df);

void ghl_func_rho2(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double Z_in,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df);

void ghl_func_1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double dummy,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df);

void ghl_func_2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double dummy,
      const double x[],
      double dx[],
      double resid[],
      double jac[][2],
      double *restrict f,
      double *restrict df);

int ghl_initialize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      harm_aux_vars_struct *restrict harm_aux,
      double *restrict Z_ptr);

int ghl_initialize_Noble_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      harm_aux_vars_struct *restrict harm_aux,
      double *restrict rho_ptr,
      double *restrict Z_ptr);

int ghl_finalize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const harm_aux_vars_struct *restrict harm_aux,
      const double Z,
      const double vsq,
      ghl_primitive_quantities *restrict prims);

int ghl_finalize_Noble_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const harm_aux_vars_struct *restrict harm_aux,
      const double Z,
      const double W,
      ghl_primitive_quantities *restrict prims);

void ghl_validate_1D(
      const double x0[1],
      double x[1]);

void ghl_validate_2D(
      const double x0[2],
      double x[2]);

#endif
