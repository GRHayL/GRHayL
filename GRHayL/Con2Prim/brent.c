#include "roots.h"

/**
 * @ingroup c2p_internal
 * @brief Ensures $\f |f(b)| < |f(a)| \f$ while finding root of \f$ f \f$.
 *
 * @param[in,out] a:  first evaluated point f(x)
 * @param[in,out] b:  second evaluated point of f(x)
 * @param[in,out] fa: f(a)
 * @param[in,out] fb: f(b)
 *
 * @returns void
 */
static inline void
ensure_b_is_closest_to_root(
    double *restrict a,
    double *restrict b,
    double *restrict fa,
    double *restrict fb) {

  if( fabs(*fa) < fabs(*fb) ) {
    swap(a, b);
    swap(fa, fb);
  }
}

/*
 * @ingroup c2p_internal
 * @brief Keeps best guesses of a,b,c for root-finding.
 *
 * @details
 * This cycles a,b = b,c. It also resets c = b.
 *
 * @param[in,out] a: bounding point to replace
 * @param[in,out] b: bounding point to keep
 * @param[in,out] c: new bounding point
 *
 * @returns void
 */
static inline void
cycle(
    double *restrict a,
    double *restrict b,
    double *restrict c) {

  *a = *b;
  *b = *c;
  *c = *a;
}

/**
 * @ingroup c2p_internal
 * @brief Used to initialize/check root-finding variables.
 *
 * @details
 * This function performs the following tasks:
 *
 * - Check if either a or b are roots of f;
 * - Check if the root is in the interval [a,b];
 * - Ensure |f(b)| < |f(a)| by swapping a and b if necessary.
 *
 * Parameters : f        - Function for which the root is computed.
 *            : fparams  - Object containing all parameters needed by the
 *                         function f other than the variable x.
 *            : a        - Lower limit of the initial interval.
 *            : b        - Upper limit of the initial interval.
 *            : fa       - f(a)
 *            : fb       - f(b)
 *            : r        - Pointer to Roots parameters.
 *
 */
static inline ghl_error_codes_t
check_a_b_compute_fa_fb(
      double f(
            const double x,
            const ghl_parameters *restrict params,
            const ghl_eos_parameters *restrict eos,
            const ghl_conservative_quantities *restrict cons_undens,
            fparams_struct *restrict fparams,
            ghl_primitive_quantities *restrict prims),
    const ghl_parameters *restrict params,
    const ghl_eos_parameters *restrict eos,
    const ghl_conservative_quantities *restrict cons_undens,
    fparams_struct *restrict fparams,
    ghl_primitive_quantities *restrict prims,
    double *restrict a,
    double *restrict b,
    double *restrict fa,
    double *restrict fb,
    roots_params *restrict r) {

  // Step 1: Compute fa; check if a is the root.
  *fa = f(*a, params, eos, cons_undens, fparams, prims);
  if(*fa == 0.0) {
    r->root     = *a;
    r->residual = *fa;
    return ghl_success;
  }

  // Step 2: Compute fb; check if b is the root.
  *fb = f(*b, params, eos, cons_undens, fparams, prims);
  if( *fb == 0.0 ) {
    r->root     = *b;
    r->residual = *fb;
    return ghl_success;
  }

  // Step 3: Ensure the root is in [a,b]
  if( (*fa)*(*fb) > 0 )
    return ghl_error_root_not_bracketed;

  // Step 4: Ensure b contains the best approximation to the root
  ensure_b_is_closest_to_root(a, b, fa, fb);

  // Step 5: If [a,b] is too small, return b
  if( fabs(*a - *b) < r->tol ) {
    r->root     = *b;
    r->residual = *fb;
    return ghl_success;
  }

  // Step 6: Root not found.
  return ghl_error_c2p_max_iter;
}

/**
 * @ingroup c2p_internal
 * @brief Find the root of f(x) in the interval [a,b] using Brent's method
 *
 * @details
 * References
 * - Press et al., Numerical Recipes, Ch. 9.3. Freely available at
 *   http://numerical.recipes/book/book.html
 * - Brent, Algorithms for Minimization Without Derivatives (1973)
 *
 * param[in] f:           function for which the root is computed
 *
 * param[in] params:      pointer to ghl_parameters
 *
 * param[in] eos:         pointer to ghl_eos_parameters
 *
 * param[in] cons_undens: pointer to ghl_conservative_quantities containing
 *                        **undensitized** conservative variables
 *
 * param[in,out] fparams: object containing all parameters needed by f
 *                        other than the variable x
 *
 * param[in,out] a:       lower limit of the initial interval
 *
 * param[in,out] b:       upper limit of the initial interval
 *
 * param[in,out] r:       pointer to roots_params
 *
 * @returns a ghl_error_codes_t code; possible codes:
 * - ghl_success if the root is found
 * - ghl_error_root_not_bracketed if the interval [a,b] does not bracket a root of f(x)
 * - ghl_error_c2p_max_iter if the maximum allowed number of iterations is exceeded
 */
ghl_error_codes_t
ghl_brent(
      double f(
            const double x,
            const ghl_parameters *restrict params,
            const ghl_eos_parameters *restrict eos,
            const ghl_conservative_quantities *restrict cons_undens,
            fparams_struct *restrict fparams,
           ghl_primitive_quantities *restrict prims),
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims,
      double a,
      double b,
      roots_params *restrict r) {

  // Step 0: Set basic info to the roots_params struct
  r->a = a;
  r->b = b;

  // Step 1: Check whether a or b is the root; compute fa and fb
  double fa, fb;
  ghl_error_codes_t error = check_a_b_compute_fa_fb(f, params, eos, cons_undens, fparams, prims, &a, &b, &fa, &fb, r);
  if(error != ghl_error_c2p_max_iter)
    return error;

  // Step 2: Declare/initialize auxiliary variables
  double c  = a;
  double fc = fa;
  double d  = b-a;
  double e  = d;
  double tol, m, P, Q, R, S;

  // Step 3: Brent's algorithm
  for(r->n_iters=1;r->n_iters<=r->max_iters;r->n_iters++) {

    // Step 3.a: Keep the bracket in [b,c]
    if(fb*fc > 0) {
      c  = a;
      fc = fa;
      d  = e = b-a;
    }

    // Step 3.b: Keep the best guess in b
    if( fabs(fc) < fabs(fb) ) {
      cycle(&a , &b , &c);
      cycle(&fa, &fb, &fc);
    }

    // Step 3.d: Set the tolerance for this iteration
    tol = 2*DBL_EPSILON*fabs(b) + 0.5*r->tol;

    // Step 3.e: Compute midpoint
    m = 0.5*(c-b);

    // Step 3.f: Check for convergence
    if(fabs(m) < tol || fb == 0.0) {
      r->root     = b;
      r->residual = fb;
      return ghl_success;
    }

    // Step 3.g: Check whether to bisect or interpolate
    if( fabs(e) < tol || fabs(fa) <= fabs(fb) )
      e = d = m; // Step 3.g.1: bisect
    else {
      // Attempt interpolation
      S = fb/fa;
      if(a == c) {
        // Step 3.g.2: Linear interpolation
        P = 2*m*S;
        Q = 1-S;
      }
      else {
        // Step 3.g.3: Inverse quadratic interpolation
        Q = fa/fc;
        R = fb/fc;
        P = S*(2*m*Q*(Q-R)-(b-a)*(R-1));
        Q = (Q-1)*(R-1)*(S-1);
      }
      if(P > 0)
        Q = -Q;
      else
        P = -P;

      // Step 3.g.3: Accept interpolation?
      if( 2*P < 3*m*Q - fabs(tol*Q) && 2*P < fabs(e*Q) ) {
        // Yes
        e = d;
        d = P/Q;
      }
      else
        e = d = m; // Interpolation failed; do a bisection
    }
    a  = b;
    fa = fb;
    if(fabs(d) > tol)
      b += d;
    else
      b += m>0 ? tol : -tol;
    fb = f(b, params, eos, cons_undens, fparams, prims);
  }

  // Step 4: The only way to get here is if we have exceeded the maximum number
  //         of iterations allowed.
  return ghl_error_c2p_max_iter;
}
