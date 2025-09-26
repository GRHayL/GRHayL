#ifndef ROOTS_H_
#define ROOTS_H_

#include "ghl.h"
#include <stdio.h>
#include <float.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.220446049250313e-16
#endif

/*
 * typedef    : fparams_struct
 * Author     : Leo Werneck
 *
 * Structure containing parameters used by f(x), the function for which we
 * wish to find the root when using the Palenzuela et al. con2prim.
 *
 * Parameters : evolve_T              - Whether or not to evolve the temperature.
 *            : temp_guess            - Temperature guess.
 *            : Y_e                   - Electron fraction.
 *            : q                     - Auxiliary quantity: tau/D.
 *            : r                     - Auxiliary quantity: S^{2}/D^{2}.
 *            : s                     - Auxiliary quantity: B^{2}/D.
 *            : t                     - Auxiliary quantity: B.S/D^{3/2}.
 *            : eos                   - Pointer to const ghl_eos_parameters struct.
 *            : cons_undens           - Pointer to const ghl_conservative_quantities struct.
 *            : compute_rho_P_eps_T_W - Function pointer provided by the user.
 */
typedef struct fparams_struct {
  bool evolve_T;
  double q, r, s, t;
  void (*compute_rho_P_eps_T_W)(
      const double x,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      struct fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims,
      double *restrict W_ptr);
  void (*compute_rho_P_eps_W)(
      const double x,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      struct fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims,
      double *restrict W_ptr);
} fparams_struct;

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
GHL_DEVICE static inline void
roots_info(const roots_params *restrict r) {

  // Step 1: Print basic message to the user
  fprintf(stderr, "(roots) Root-finding information:\n");
  fprintf(stderr, "(roots)   %16s : %s\n", "Routine", r->routine_name);
  fprintf(stderr, "(roots)   %16s : [%c%21.15e, %c%21.15e]\n",
         "Initial interval",
         r->a >= 0 ? '+' : '-', fabs(r->a),
         r->b >= 0 ? '+' : '-', fabs(r->b));
  fprintf(stderr, "(roots)   %16s : ", "Status");
  fprintf(stderr, "(roots)   %16s : %d\n", "Iterations", r->n_iters);
  fprintf(stderr, "(roots)   %16s : %.15e\n", "Root", r->root);
  fprintf(stderr, "(roots)   %16s : %.15e\n", "Residual", r->residual);
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
GHL_DEVICE static inline void
swap(
    double *restrict a,
    double *restrict b) {

  const double c = *a;
  *a = *b;
  *b = c;
}



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
GHL_DEVICE static inline int
sign(const double x) {

  return (x>0) - (x<0);
}

/*
 * Function prototypes
 */
GHL_DEVICE
ghl_error_codes_t ghl_brent(
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
      roots_params *restrict r);

#endif // ROOTS_H_
