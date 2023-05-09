#ifndef ROOTS_H_
#define ROOTS_H_

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.220446049250313e-16
#endif

/*
 * Error keys for root-finding routines
 */
typedef enum {
  roots_continue=-1,
  roots_success,
  roots_error_root_not_bracketed,
  roots_error_max_iter,
  roots_error_not_finite
} roots_error_t;

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
typedef struct roots_params {
  roots_error_t error_key;
  char routine_name[256];
  unsigned n_iters, max_iters;
  double a, b, root, residual, tol;
} roots_params;

/*
 * Function   : roots_info
 * Author     : Leo Werneck
 *
 * Prints information about the root-finding process.
 *
 * Parameters : r        - Pointer to roots parameters.
 *
 * Returns    : Nothing.
 */
static inline void
roots_info(const roots_params *restrict r) {

  // Step 1: Print basic message to the user
  fprintf(stderr, "(roots) Root-finding information:\n");
  fprintf(stderr, "(roots)   %16s : %s\n", "Routine", r->routine_name);
  fprintf(stderr, "(roots)   %16s : [%c%21.15e, %c%21.15e]\n",
         "Initial interval",
         r->a >= 0 ? '+' : '-', fabs(r->a),
         r->b >= 0 ? '+' : '-', fabs(r->b));
  fprintf(stderr, "(roots)   %16s : ", "Status");
  switch(r->error_key) {
    case roots_continue:
      break;
    case roots_success:
      fprintf(stderr, "Success\n");
      break;
    case roots_error_root_not_bracketed:
      fprintf(stderr, "Failure\n");
      fprintf(stderr, "(roots)   %16s : ", "Error message");
      fprintf(stderr, "Initial interval does not bracket the root.\n");
      break;
    case roots_error_max_iter:
      fprintf(stderr, "Failure\n");
      fprintf(stderr, "(roots)   %16s : ", "Error message");
      fprintf(stderr, "Maximum number of iterations (%d) exceeded.\n", r->max_iters);
      break;
    case roots_error_not_finite:
      fprintf(stderr, "Failure\n");
      fprintf(stderr, "(roots)   %16s : ", "Error message");
      fprintf(stderr, "Found NAN or INF during root-finding procedure.\n");
      break;
  }

  // Step 2: If succeeded, print detailed success message
  if( !r->error_key ) {
    fprintf(stderr, "(roots)   %16s : %d\n", "Iterations", r->n_iters);
    fprintf(stderr, "(roots)   %16s : %.15e\n", "Root", r->root);
    fprintf(stderr, "(roots)   %16s : %.15e\n", "Residual", r->residual);
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
static inline void
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
 * From inputs a, b, c, cicle a, b, c = b, c, b.
 *
 * Parameters : a        - First number.
 *            : b        - Second number.
 *            : c        - Third number.
 *
 * Returns    : Nothing.
 */
// static inline void
// cicle(
//     double *restrict a,
//     double *restrict b,
//     double *restrict c ) {

//   *a = *b;
//   *b = *c;
//   *c = *a;
// }

/*
 * Function   : sign
 * Author     : Leo Werneck
 *
 * Returns the sign of a number x.
 *
 * Parameters : x        - Number
 *            : b        - Second number.
 *
 * Returns    : +1 if x>=0, -1 otherwise.
 */
static inline int
sign( const double x ) {

  return (x>0) - (x<0);
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
 * Returns    : One the following error keys:
 *                 - roots_success if the root is found
 *                 - roots_continue if the root is not found but no errors
 *                   occurred
 *                 - roots_error_root_not_bracketed if the interval [a,b]
 *                   does not bracket a root of f(x)
 */
// static inline roots_error_t
// check_a_b_compute_fa_fb(
//     double f(double const, void *restrict),
//     void   *restrict fparams,
//     double *restrict a,
//     double *restrict b,
//     double *restrict fa,
//     double *restrict fb,
//     roots_params *restrict r ) {

//   // Step 1: Compute fa; check if a is the root.
//   *fa = f(*a, fparams);
//   if( *fa == 0.0 ) {
//     r->root     = *a;
//     r->residual = *fa;
//     return (r->error_key = roots_success);
//   }

//   // Step 2: Compute fb; check if b is the root.
//   *fb = f(*b, fparams);
//   if( *fb == 0.0 ) {
//     r->root     = *b;
//     r->residual = *fb;
//     return (r->error_key = roots_success);
//   }

//   // Step 3: Ensure the root is in [a,b]
//   if( (*fa)*(*fb) > 0 )
//     return (r->error_key = roots_error_root_not_bracketed);

//   // Step 4: Ensure b contains the best approximation to the root
//   ensure_b_is_closest_to_root(a, b, fa, fb);

//   // Step 5: If [a,b] is too small, return b
//   if( fabs(*a - *b) < r->tol ) {
//     r->root     = *b;
//     r->residual = *fb;
//     return (r->error_key = roots_success);
//   }

//   // Step 6: Root not found.
//   return roots_continue;
// }


/*
 * Function prototypes
 */
// roots_error_t
// brent(
//     double f(const double, void *restrict),
//     void *restrict fparams,
//     double a,
//     double b,
//     roots_params *restrict r );

roots_error_t
toms748(
    double f(const double, void *restrict),
    void *restrict fparams,
    double a,
    double b,
    roots_params *restrict r );

#endif // ROOTS_H_
