#ifndef PALENZUELA_H_
#define PALENZUELA_H_

#include "ghl_con2prim.h"
#include "roots.h"

/**
 * @ingroup c2p_internal
 * @brief Computes rho and W from x and the conservative variables using Eq. (42)
 *        of https://arxiv.org/pdf/1712.07538.pdf and rho = D/W.
 * @todo
 * Add in/out information to params
 *
 * @param x:       pointer to fparams_struct.
 *
 * @param fparams: array containing v_{i}.
 *
 * @param rho_ptr: stores the output of rho.
 *
 * @param W_ptr:   stores the output of W.
 *
 * @returns void
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
  Wminus2        = fmin(fmax(Wminus2, params->inv_sq_max_Lorentz_factor ), 1.0);
  const double W = pow(Wminus2, -0.5);

  // Step 3: Compute rho
  const double rho = cons_undens->rho/W;

  // Step 4: Set the output
  *rho_ptr = rho;
  *W_ptr   = W;
}

/**
 * @ingroup c2p_internal
 * @brief Computes auxiliary quantities of \f$ B^i \f$ and \f$ S_i \f$
 *
 * @details
 * This function computes the quantities \f$ S^i \f$, \f$ B^2 \f$, \f$ S^2 \f$,
 * and \f$ B\cdot S \f$ using the functions @ref ghl_raise_lower_vector_3D and
 * @ref ghl_compute_vec2_from_vec3D. It also enforces that
 *
 * \f[
 * S^2 < \left| \gamma \right| \left( \tilde{\tau} + \rho_* \right)^2
 * \f]
 *
 * @todo
 * I (SCupp) compared this limit to @ref ghl_apply_conservative_limits, and the
 * two seem to disagree. Should this \f$ S_i \f$ limiter just do the same as
 * @ref ghl_apply_conservative_limits? It seems to me like one of them should be
 * wrong unless someone can explain why they would differ between EOS. This
 * should be looked at and, if true, that function should be properly extended
 * to work for any EOS and the two limits should converge to the 'correct'
 * version. It also references eq. A5 of [1], but I don't know what [1]
 * represents here.
 *
 * @param[in] ADM_metric:  pointer to ghl_metric_quantities containing the ADM metric
 *
 * @param[in] cons_undens: pointer to ghl_conservative_quantities containing
 *                         **undensitized** conservative variables
 *
 * @param[in] prims:       pointer to ghl_primitive_quantities (only magnetic field is used)
 *
 * @param[out] SU:         returned quantity \f$ S^i \f$
 *
 * @param[out] Bsq:        returned quantity \f$ B^2 \f$
 *
 * @param[out] Ssq:        returned quantity \f$ S^2 \f$
 *
 * @param[out] BdotS:      returned quantity \f$ B \cdot S \f$
 *
 * @returns void
 */
static inline void
compute_SU_Bsq_Ssq_BdotS(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      double SU[restrict 3],
      double *restrict Bsq,
      double *restrict Ssq,
      double *restrict BdotS) {

  // Step 1: Compute S^{2} = gamma^{ij}S_{i}S_{j}
  double SD[3] = {cons_undens->SD[0], cons_undens->SD[1], cons_undens->SD[2]};
  double S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);

  // Step 2: Enforce ceiling on S^{2} (Eq. A5 of [1])
  // Step 2.1: Compute maximum allowed value for S^{2}
  const double S_squared_max = SQR(cons_undens->tau + cons_undens->rho);
  if(S_squared > S_squared_max) {
    // Step 2.2: Rescale S_{i}
    const double rescale_factor = sqrt(0.9999*S_squared_max/S_squared);
    for(int i=0;i<3;i++) {
      SD[i] *= rescale_factor;
    }

    // Step 2.3: Recompute S^{2}
    S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);
  }

  *Ssq = S_squared;

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  *Bsq = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims->BU);

  // Step 4: Compute B.S = B^{i}S_{i}
  *BdotS = prims->BU[0]*SD[0] + prims->BU[1]*SD[1] + prims->BU[2]*SD[2];

  // Step 5: Compute S^{i}
  ghl_raise_lower_vector_3D(ADM_metric->gammaUU, SD, SU);
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
      const ghl_metric_quantities *restrict ADM_metric,
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
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics);

#endif // PALENZUELA_H_
