#ifndef PALENZUELA_H_
#define PALENZUELA_H_

#include "ghl_con2prim.h"
#include "roots.h"

/**
 * @ingroup c2p_internal
 * @brief Compute rest-mass density and Lorentz factor from a Palenzuela iterate.
 * @author Leo Werneck
 *
 * Computes rho and W from x and the conservative variables using Eq. (42)
 * of https://arxiv.org/pdf/1712.07538.pdf and rho = D/W.
 *
 * @param[in] x Trial value of h W
 * @param[in] params Parameters supplying the maximum Lorentz factor
 * @param[in] cons_undens Undensitized conservatives supplying D = rho W
 * @param[in] fparams Normalized contractions r, s, and t from roots.h
 * @param[out] rho_ptr Rest-mass density D/W
 * @param[out] W_ptr Lorentz factor, limited to [1, max_Lorentz_factor]
 */
static inline void
compute_rho_W_from_x_and_conservatives(
  const double x,
  const ghl_parameters *restrict params,
  const ghl_conservative_quantities *restrict cons_undens,
  const fparams_struct *restrict fparams,
  double *restrict rho_ptr,
  double *restrict W_ptr) {

  // Step 1: Unpack the fparams struct
  const double r = fparams->r;
  const double s = fparams->s;
  const double t = fparams->t;

  // Step 2: Compute W
  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t )/ (x*x*(x+s)*(x+s));
  Wminus2        = ghl_clamp(Wminus2, params->inv_sq_max_Lorentz_factor, 1.0);
  const double W = pow(Wminus2, -0.5);

  // Step 3: Compute rho
  const double rho = cons_undens->rho/W;

  // Step 4: Set the output
  *rho_ptr = rho;
  *W_ptr   = W;
}

/*
 * Function prototypes
 */
ghl_error_codes_t ghl_tabulated_Palenzuela1D(
      void compute_rho_P_eps_T_W(
            const double x,
            const ghl_parameters *restrict params,
            const ghl_eos_parameters *restrict eos,
            const ghl_conservative_quantities *restrict cons_undens,
            fparams_struct *restrict fparams,
            ghl_primitive_quantities *restrict prims,
            double *restrict W_ptr),
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

ghl_error_codes_t ghl_hybrid_Palenzuela1D(
      void compute_rho_P_eps_W(
            const double x,
            const ghl_parameters *restrict params,
            const ghl_eos_parameters *restrict eos,
            const ghl_conservative_quantities *restrict cons_undens,
            fparams_struct *restrict fparams,
            ghl_primitive_quantities *restrict prims,
            double *restrict W_ptr),
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

#endif // PALENZUELA_H_
