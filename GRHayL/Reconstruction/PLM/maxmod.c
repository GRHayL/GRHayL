#include "ghl_reconstruction.h"

/**
 * @ingroup recon_internal
 * @brief Evaluates \f$ \Delta U \f$ using maxmod.
 *
 * @details
 * This function computes the value of \f$ \Delta U \f$  for use in
 * @ref ghl_superbee_reconstruction. It uses the logic
 * \f[
 * result = 
 * \begin{cases}
 *   a & \text{if } |a| > |b| \text{ and } a b > 0 \\
 *   b & \text{if } |a| < |b| \text{ and } a b > 0 \\
 *   0 & \text{if } a b \leq 0
 * \end{cases}
 * \f]
 * to determine what should be used to approximate \f$ \Delta U \f$.
 *
 * @param[in] a: first approximation of \f$ \Delta U \f$
 *
 * @param[in] b: second approximation of \f$ \Delta U \f$
 *
 * @returns approximate \f$ \Delta U \f$ value
 */
double ghl_maxmod(
      const double a,
      const double b) {

  const double ab = a*b;

  if(ab > 0) {
    if(fabs(a) > fabs(b)) {
      return a;
    } else {
      return b;
    }
  } else {
    return 0.0;
  }
}
