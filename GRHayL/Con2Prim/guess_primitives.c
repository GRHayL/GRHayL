#include "ghl_con2prim.h"

/**
 * @ingroup Con2Prim
 * @brief Computes an initial guess for the @ref Con2Prim solvers
 *
 * @details
 * This function sets a default initial guess for when a better guess is not
 * available. It sets
 * 
 * \f[
 * \begin{aligned}
 * \rho &= \frac{\rho_*}{\sqrt{|\gamma|}} \\
 * u^0 &= 1 \\
 * v^i &= -\beta^i \\
 * Y_e &= \frac{\tilde{Y_e}}{\rho_*} \\
 * T &= T_\mathrm{max}
 * \end{aligned}
 * \f]
 *
 * Note that this value of \f$ v^i \f$ implies that \f$ u^i=0 \f$, which leads
 * to the given value of \f$ \rho \f$. Also, the quantity \f$ T_\mathrm{max} \f$
 * uses ghl_eos_parameters::T_max. Finally, for hybrid or simple EOS,
 * we set the pressure and specific internal energy \f$ \epsilon \f$ to the
 * cold values.
 *
 * @param[in] eos: pointer to ghl_eos_parameters struct
 *
 * @param[in] ADM_metric: pointer to ghl_metric_quantities struct with ADM metric
 *
 * @param[in] cons_undens: pointer to ghl_conservative_quantities struct with
 *                         **undensitized** conservative variables
 *
 * @param[out] prims: pointer to ghl_primitive_quantities containing initial guess
 *
 * @returns void
 */
void ghl_guess_primitives(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims) {

  // Use atmosphere as initial guess:
  prims->rho         = cons_undens->rho;
  prims->u0          = 1.0;
  prims->vU[0]       = -ADM_metric->betaU[0];
  prims->vU[1]       = -ADM_metric->betaU[1];
  prims->vU[2]       = -ADM_metric->betaU[2];
  prims->Y_e         = cons_undens->Y_e/cons_undens->rho;
  prims->temperature = eos->T_max;

  // If using hybrid or ideal fluid EOS, compute P_cold and eps_cold
  if( eos->eos_type == ghl_eos_hybrid || eos->eos_type == ghl_eos_simple )
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &prims->press, &prims->eps);
}
