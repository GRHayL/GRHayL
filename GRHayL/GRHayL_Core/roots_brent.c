#include "ghl_roots.h"

/*
 * Function   : ensure_b_is_closest_to_root
 * Author     : Leo Werneck
 *
 * Given a, b, f(a), f(b), ensures |f(b)| < |f(a)|.
 *
 * Parameters : a        - First point where f(x) is evaluated at.
 *            : b        - Second point where f(x) is evaluated at.
 *            : fa       - f(a).
 *            : fb       - f(b).
 *
 * Returns    : Nothing.
 */
static inline void ensure_b_is_closest_to_root(
      double *restrict a,
      double *restrict b,
      double *restrict fa,
      double *restrict fb) {

  if(fabs(*fa) < fabs(*fb)) {
    swap(a, b);
    swap(fa, fb);
  }
}

/*
 * Function   : check_a_b_compute_fa_fb
 * Author     : Leo Werneck
 *
 * This function is used at the beginning of the root-finding methods. It
 * performs the following tasks:
 *
 *   1. Check if either a or b are roots of f;
 *   2. Check if the root is in the interval [a,b];
 *   3. Ensure |f(b)| < |f(a)| by swapping a and b if necessary.
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
    double f(const double, void *restrict),
    void *restrict fparams,
    double *restrict a,
    double *restrict b,
    double *restrict fa,
    double *restrict fb,
    roots_params *restrict r) {

  // Step 1: Compute fa; check if a is the root.
  *fa = f(*a, fparams);
  if(*fa == 0.0) {
    r->root = *a;
    r->residual = *fa;
    return ghl_success;
  }

  // Step 2: Compute fb; check if b is the root.
  *fb = f(*b, fparams);
  if( fabs(*fb) < 1e-15 ) {
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
  if(fabs(*a - *b) < r->tol) {
    r->root = *b;
    r->residual = *fb;
    return ghl_success;
  }

  // Step 6: Root not found.
  return ghl_error_c2p_max_iter;
}

/*
 * Function   : ghl_brent
 * Author     : Leo Werneck
 *
 * Find the root of f(x) in the interval [a,b] using Brent's method.
 *
 * Parameters : f        - Function for which the root is computed.
 *            : fparams  - Object containing all parameters needed by the
 *                         function f other than the variable x.
 *            : a        - Lower limit of the initial interval.
 *            : b        - Upper limit of the initial interval.
 *            : r        - Pointer to roots parameters (see ghl_roots.h).
 *                         The root is stored in r->root.
 *
 * Returns    : One the following error keys:
 *                 - roots_success if the root is found
 *                 - roots_error_root_not_bracketed if the interval [a,b]
 *                   does not bracket a root of f(x)
 *                 - roots_error_max_iter if the maximum allowed number of
 *                   iterations is exceeded
 *
 * References : Press et al., Numerical Recipes, Ch. 9.3
 *              Freely available at: http://numerical.recipes/book/book.html
 *            : Brent, Algorithms for Minimization Without Derivatives (1973)
 */
ghl_error_codes_t
ghl_brent(
      double f(const double, void *restrict),
      void *restrict fparams,
      double a,
      double b,
      roots_params *restrict r) {

  // Step 0: Set basic info to the roots_params struct
  sprintf(r->routine_name, __func__);
  r->a = a;
  r->b = b;

  // Step 1: Check whether a or b is the root; compute fa and fb
  double fa, fb;
  ghl_error_codes_t error = check_a_b_compute_fa_fb(f, fparams, &a, &b, &fa, &fb, r);
  if(error != ghl_error_c2p_max_iter)
    return error;

  // Step 2: Declare/initialize auxiliary variables
  double c = a;
  double fc = fa;
  double d = b - a;
  double e = d;
  double tol, m, P, Q, R, S;
  // Step 3: Brent's algorithm
  for(r->n_iters = 1; r->n_iters <= r->max_iters; r->n_iters++) {
    // Step 3.a: Keep the bracket in [b,c]
    if(fb * fc > 0) {
      c = a;
      fc = fa;
      d = e = b - a;
    }

    // Step 3.b: Keep the best guess in b
    if(fabs(fc) < fabs(fb)) {
      swap(&b, &c);
      swap(&fb, &fc);
      a = c;
      fa = fc;
    }

    // Step 3.c: Set the tolerance for this iteration
    tol = 2 * DBL_EPSILON * fabs(b) + 0.5 * r->tol;

    // Step 3.d: Compute midpoint
    m = 0.5 * (c - b);

    // Step 3.e: Check for convergence
    if(fabs(m) < tol || fb == 0.0) {
      r->root = b;
      r->residual = fb;
      return ghl_success;
    }

    // Step 3.f: Check whether to bisect or interpolate
    if(fabs(e) < tol || fabs(fa) <= fabs(fb)) {
      e = d = m; // Step 3.f.1: bisect
    }
    else {
      // Attempt interpolation
      S = fb / fa;
      if(a == c) {
        // Step 3.f.2: Linear interpolation
        P = 2 * m * S;
        Q = 1 - S;
      }
      else {
        // Step 3.f.3: Inverse quadratic interpolation
        Q = fa / fc;
        R = fb / fc;
        P = S * (2 * m * Q * (Q - R) - (b - a) * (R - 1));
        Q = (Q - 1) * (R - 1) * (S - 1);
      }
      if(P > 0) {
        Q = -Q;
      }
      else {
        P = -P;
      }

      // Step 3.f.4: Accept interpolation?
      if(2 * P < 3 * m * Q - fabs(tol * Q) && 2 * P < fabs(e * Q)) {
        // Yes
        e = d;
        d = P / Q;
      }
      else {
        e = d = m; // Interpolation failed; do a bisection
      }
    }
    a = b;
    fa = fb;
    if(fabs(d) > tol) {
      b += d;
    }
    else {
      b += m>0 ? tol : -tol;
    }
    fb = f(b, fparams);
  }

  // Step 4: The only way to get here is if we have exceeded the maximum number
  //         of iterations allowed.
  return ghl_error_c2p_max_iter;
}
