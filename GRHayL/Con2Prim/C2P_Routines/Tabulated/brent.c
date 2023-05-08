#include "roots.h"

/*
 * Function   : brent
 * Author     : Leo Werneck
 *
 * Find the root of f(x) in the interval [a,b] using Brent's method.
 *
 * Parameters : f        - Function for which the root is computed.
 *            : fparams  - Object containing all parameters needed by the
 *                         function f other than the variable x.
 *            : a        - Lower limit of the initial interval.
 *            : b        - Upper limit of the initial interval.
 *            : r        - Pointer to roots parameters (see roots.h).
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
roots_error_t
brent(
    double f(double const, void *restrict),
    void *restrict fparams,
    double a,
    double b,
    roots_params *restrict r ) {

  // Step 0: Set basic info to the roots_params struct
  sprintf(r->routine_name, __func__);
  r->a = a;
  r->b = b;

  // Step 1: Check whether a or b is the root; compute fa and fb
  double fa, fb;
  if( check_a_b_compute_fa_fb(f, fparams, &a, &b, &fa, &fb, r) >= roots_success )
    return r->error_key;

  // Step 2: Declare/initialize auxiliary variables
  double c  = a;
  double fc = fa;
  double d  = b-a;
  double e  = d;
  double tol, m, P, Q, R, S;

  // Step 3: Brent's algorithm
  for(r->n_iters=1;r->n_iters<=r->max_iters;r->n_iters++) {

    // Step 3.a: Keep the bracket in [b,c]
    if( fb*fc > 0 ) {
      c  = a;
      fc = fa;
      d  = e = b-a;
    }

    // Step 3.b: Keep the best guess in b
    if( fabs(fc) < fabs(fb) ) {
      cicle(&a , &b , &c );
      cicle(&fa, &fb, &fc);
    }

    // Step 3.d: Set the tolerance for this iteration
    tol = 2*DBL_EPSILON*fabs(b) + 0.5*r->tol;

    // Step 3.e: Compute midpoint
    m = 0.5*(c-b);

    // Step 3.f: Check for convergence
    if( fabs(m) < tol || fb == 0.0 ) {
      r->root     = b;
      r->residual = fb;
      return (r->error_key = roots_success);
    }

    // Step 3.g: Check whether to bisect or interpolate
    if( fabs(e) < tol || fabs(fa) <= fabs(fb) )
      e = d = m; // Step 3.g.1: bisect
    else {
      // Attempt interpolation
      S = fb/fa;
      if( a == c ) {
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
      if( P > 0 )
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
    if( fabs(d) > tol )
      b += d;
    else
      b += m>0 ? tol : -tol;
    fb = f(b, fparams);
  }

  // Step 4: The only way to get here is if we have exceeded the maximum number
  //         of iterations allowed.
  return (r->error_key = roots_error_max_iter);
}
