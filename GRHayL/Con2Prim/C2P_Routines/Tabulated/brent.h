#ifndef BRENT_H_
#define BRENT_H_

#include <stdio.h>
#include <math.h>

#ifdef DBL_EPSILON
#undef DBL_EPSILON
#endif
#define DBL_EPSILON 2.220446049250313e-16

/*
 * Error keys for Brent root-finding routine
 */
typedef enum {
  brent_continue=-1,
  brent_success,
  brent_error_root_not_bracketed,
  brent_error_max_iter
} brent_error_t;

/*
 * Parameters for Brent root-finding routine
 *
 * error_key : Status of the root-finding algorithm.
 * n_iters   : Number of iterations used to find the root.
 * max_iters : Maximum allowed iterations to find the root.
 * a         : Initial interval is [a,b].
 * b         : Initial interval is [a,b].
 * root      : The root.
 * residual  : f(root).
 * tol       : Stop when fabs(b_n-a_n) < tol.
 */
typedef struct brent_params {
  brent_error_t error_key;
  unsigned int n_iters, max_iters;
  double a, b, root, residual, tol;
} brent_params;

/*
 * Function   : brent_info
 * Author     : Leo Werneck
 *
 * Prints information about the root-finding process.
 *
 * Parameters : r        - Pointer to Brent parameters.
 *
 * Returns    : Nothing.
 */
static inline void
brent_info(const brent_params *restrict r) {

  // Step 1: Print basic message to the user
  printf("(brent) Root-finding information:\n");
  printf("(brent)   %16s : Brent's\n", "Method");
  printf("(brent)   %16s : [%c%21.15e, %c%21.15e]\n",
         "Initial interval",
         r->a >= 0 ? '+' : '-', fabs(r->a),
         r->b >= 0 ? '+' : '-', fabs(r->b));
  printf("(brent)   %16s : ", "Status");
  switch(r->error_key) {
    case brent_continue:
      break;
    case brent_success:
      printf("Success\n");
      break;
    case brent_error_root_not_bracketed:
      printf("Failure\n");
      printf("(brent)   %16s : ", "Error message");
      printf("Initial interval does not bracket the root.\n");
      break;
    case brent_error_max_iter:
      printf("Failure\n");
      printf("(brent)   %16s : ", "Error message");
      printf("Maximum number of iterations (%d) exceeded.\n", r->max_iters);
      break;
  }

  // Step 2: If succeeded, print detailed success message
  if( !r->error_key ) {
    printf("(brent)   %16s : %d\n", "Iterations", r->n_iters);
    printf("(brent)   %16s : %.15e\n", "Root", r->root);
    printf("(brent)   %16s : %.15e\n", "Residual", r->residual);
  }
}

/*
 * Function   : swap
 * Author     : Leo Werneck
 *
 * Swaps the values of two doubles a and b.
 *
 * Parameters : a        - First number.
 *            : b        - Second number.
 *
 * Returns    : Nothing.
 */
static void
swap(
    double *restrict a,
    double *restrict b ) {

  const double c = *a;
  *a = *b;
  *b = c;
}

/*
 * Function   : cicle
 * Author     : Leo Werneck
 *
 * From inputs a, b, c, cicle a <- b <- c <- a.
 *
 * Parameters : a        - First number.
 *            : b        - Second number.
 *            : c        - Third number.
 *
 * Returns    : Nothing.
 */
static void
cicle(
    double *restrict a,
    double *restrict b,
    double *restrict c ) {

  const double d = *a;
  *a = *b;
  *b = *c;
  *c = d;
}

/*
 * Function   : ensure_b_is_closest_to_root
 * Author     : Leo Werneck
 *
 * Given a, b, f(a), f(b), ensures |f(b)| < |f(a)|.
 *
 * Parameters : a        - First point where f(x) is evaluated at.
 *            : b        - Second point where f(x) is evaluated at.
 * Parameters : fa       - f(a)
 *            : fb       - f(b).
 *
 * Returns    : Nothing.
 */
static void
ensure_b_is_closest_to_root(
    double *restrict a,
    double *restrict b,
    double *restrict fa,
    double *restrict fb ) {

  if( fabs(*fa) < fabs(*fb) ) {
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
 *            : r        - Pointer to Brent parameters.
 *
 * Returns    : One the following error keys:
 *                 - brent_success if the root is found
 *                 - brent_continue if the root is not found but no errors
 *                   occurred
 *                 - brent_error_root_not_bracketed if the interval [a,b]
 *                   does not bracket a root of f(x)
 */
static brent_error_t
check_a_b_compute_fa_fb(
    double f(double const, void *restrict),
    void   *restrict fparams,
    double *restrict a,
    double *restrict b,
    double *restrict fa,
    double *restrict fb,
    brent_params *restrict r ) {

  // Step 1: Compute fa; check if a is the root.
  *fa = f(*a, fparams);
  if( *fa == 0.0 ) {
    r->root     = *a;
    r->residual = *fa;
    return (r->error_key = brent_success);
  }

  // Step 2: Compute fb; check if b is the root.
  *fb = f(*b, fparams);
  if( *fb == 0.0 ) {
    r->root     = *b;
    r->residual = *fb;
    return (r->error_key = brent_success);
  }

  // Step 3: Ensure the root is in [a,b]
  if( (*fa)*(*fb) > 0 )
    return (r->error_key = brent_error_root_not_bracketed);

  // Step 4: Ensure b contains the best approximation to the root
  ensure_b_is_closest_to_root(a, b, fa, fb);

  // Step 5: If [a,b] is too small, return b
  if( fabs(*a - *b) < r->tol ) {
    r->root     = *b;
    r->residual = *fb;
    return (r->error_key = brent_success);
  }

  // Step 6: Root not found.
  return brent_continue;
}

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
 *            : r        - Pointer to brent library parameters (see brent.h).
 *                         The root is stored in r->root.
 *
 * Returns    : One the following error keys:
 *                 - brent_success if the root is found
 *                 - brent_error_root_not_bracketed if the interval [a,b]
 *                   does not bracket a root of f(x)
 *                 - brent_error_max_iter if the maximum allowed number of
 *                   iterations is exceeded
 *
 * References : Press et al., Numerical Recipes, Ch. 9.3
 *              Freely available at: http://numerical.recipes/book/book.html
 *            : Brent, Algorithms for Minimization Without Derivatives (1973)
 */
static brent_error_t
brent(
    double f(double const, void *restrict),
    void *restrict fparams,
    double a,
    double b,
    brent_params *restrict r ) {

  // Step 0: Set basic info to the brent_params struct
  r->a = a;
  r->b = b;

  // Step 1: Check whether a or b is the root; compute fa and fb
  double fa, fb;
  if( check_a_b_compute_fa_fb(f, fparams, &a, &b, &fa, &fb, r) >= brent_success )
    return r->error_key;

  // Step 2: Declare/initialize auxiliary variables
  double c  = b;
  double fc = fb;
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
    if( fabs(b-a) < tol || fb == 0.0 ) {
      r->root     = b;
      r->residual = fb;
      return (r->error_key = brent_success);
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
  return (r->error_key = brent_error_max_iter);
}

#endif // BRENT_H_
