#ifndef ROOTS_H_
#define ROOTS_H_

#include "ghl.h"
#include <stdio.h>
#include <float.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.220446049250313e-16
#endif

/**
 * @ingroup c2p_internal
 * @brief Structure containing parameters used by f(x)
 *
 * @details
 * This struct contains several values needed by the function for which we
 * wish to find the root when using the Palenzuela et al. con2prim.
 *
 * @todo
 * I think most functions should have ghl_parameters, so I would remove evolve_T from
 * here. Also, I would like to see us move away from function pointers as much as
 * possible to move towards GPU support. Reworking so we don't need to carry around
 * these fpointers would be great.
 */
typedef struct fparams_struct {
  /** Whether or not to evolve the temperature. */
  bool evolve_T;
  /** Auxiliary quantity \f$ \frac{\tau}{D} \f$ */
  double q;
  /** Auxiliary quantity \f$ \frac{S^2}{D^2} \f$ */
  double r;
  /** Auxiliary quantity \f$ \frac{B^2}{D} \f$ */
  double s;
  /** Auxiliary quantity \f$ \frac{B \cdot S}{D^{3/2}} \f$ */
  double t;
  /** Function pointer provided by the surrounding functions */
  void (*compute_rho_P_eps_T_W)(
      const double x,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      struct fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims,
      double *restrict W_ptr);
  /** Function pointer provided by surrounding functions */
  void (*compute_rho_P_eps_W)(
      const double x,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      struct fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims,
      double *restrict W_ptr);
} fparams_struct;

/**
 * @ingroup c2p_internal
 * @brief Parameters for Brent root-finding routine
 *
 * @todo
 * Several of these are included in ghl_parameters and ghl_con2prim_diagnostics.
 * I think these duplications should be reduced as much as possible.
 */
typedef struct roots_params {
  unsigned n_iters;   /**< Number of iterations used to find the root */
  unsigned max_iters; /**< Maximum allowed iterations to find the root */
  double a;           /**< Initial floor of interval */
  double b;           /**< Initial ceiling of interval */
  double root;        /**< The root */
  double residual;    /**< f(root) */
  double tol;         /**< Tolerance to trigger success: \f$ |b_n - a_n| < tol \f$ */
} roots_params;

/**
 * @ingroup c2p_internal
 * @brief Prints information about the root-finding process.
 *
 * @param[in] r: pointer to roots_params
 *
 * @param[in] calling_function: Calling function (just pass \_\_func\_\_)
 *
 * @returns void
 */
static inline void
roots_info(
      const roots_params *restrict r,
      const char *calling_function) {

  // Step 1: Print basic message to the user
  fprintf(stderr, "(roots) Root-finding information:\n");
  fprintf(stderr, "(roots)   %16s : %s\n", "Routine", calling_function);
  fprintf(stderr, "(roots)   %16s : [%c%21.15e, %c%21.15e]\n",
         "Initial interval",
         r->a >= 0 ? '+' : '-', fabs(r->a),
         r->b >= 0 ? '+' : '-', fabs(r->b));
  fprintf(stderr, "(roots)   %16s : ", "Status");
  fprintf(stderr, "(roots)   %16s : %d\n", "Iterations", r->n_iters);
  fprintf(stderr, "(roots)   %16s : %.15e\n", "Root", r->root);
  fprintf(stderr, "(roots)   %16s : %.15e\n", "Residual", r->residual);
}

/**
 * @ingroup c2p_internal
 * @brief Swaps the value of two doubles
 *
 * @param[in,out] a: first number
 * @param[in,out] b: second number
 *
 * @returns void
 */
static inline void
swap(
    double *restrict a,
    double *restrict b) {

  const double c = *a;
  *a = *b;
  *b = c;
}

/**
 * @ingroup c2p_internal
 * @brief Returns the sign of a double
 *
 * @param[in] x: value
 *
 * @returns +1 if x>=0, -1 otherwise
 */
static inline int
sign(const double x) {

  return (x>0) - (x<0);
}

/*
 * Function prototypes
 */
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

ghl_error_codes_t ghl_toms748(
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
