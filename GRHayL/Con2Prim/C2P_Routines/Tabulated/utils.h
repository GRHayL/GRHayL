#ifndef PALENZUELA_H_
#define PALENZUELA_H_

#include "con2prim.h"
#include "roots.h"

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
 *            : eos                   - Pointer to const eos_parameters struct.
 *            : cons_undens           - Pointer to const conservative_quantities struct.
 *            : compute_rho_P_eps_T_W - Function pointer provided by the user.
 */
typedef struct fparams_struct {
  bool evolve_T;
  double temp_guess, Y_e, q, r, s, t;
  const eos_parameters *eos;
  const conservative_quantities *cons_undens;
  void (*compute_rho_P_eps_T_W)(
      const double x,
      struct fparams_struct *restrict fparams,
      double *restrict rho_ptr,
      double *restrict P_ptr,
      double *restrict eps_ptr,
      double *restrict T_ptr,
      double *restrict W_ptr );
} fparams_struct;

/*
 * Function   : compute_S_squared
 * Author     : Leo Werneck
 *
 * Computes S^{2} = gamma^{ij}S_{i}S_{j}.
 *
 * Parameters : metric   - Pointer to GRHayL metric_quantities struct.
 *            : SD       - Array containing S_{i}.
 *
 * Returns    : S^{2}
 */
static inline double
compute_S_squared(
      const metric_quantities *restrict metric,
      const double *restrict SD ) {

  return metric->adm_gupxx * SD[0] * SD[0] +
         metric->adm_gupyy * SD[1] * SD[1] +
         metric->adm_gupzz * SD[2] * SD[2] +
   2.0 * metric->adm_gupxy * SD[0] * SD[1] +
   2.0 * metric->adm_gupxz * SD[0] * SD[2] +
   2.0 * metric->adm_gupyz * SD[1] * SD[2];
}

/*
 * Function   : raise_vector_3d
 * Author     : Leo Werneck
 *
 * Raises a 3-vector v^{i} = gamma^{ij} v_{j}.
 *
 * Parameters : metric   - Pointer to GRHayL metric_quantities struct.
 *            : vD       - Array containing v_{i}.
 *            : vU       - Array containing v^{i}, the output.
 *
 * Returns    : Nothing.
 */
static inline void
raise_vector_3d(
      const metric_quantities *restrict metric,
      const double *restrict vD,
      double *restrict vU ) {

  // v^{x} = gamma^{xj}v_{j} = gamma^{jx}v_{j}
  vU[0] = metric->adm_gupxx * vD[0]
        + metric->adm_gupxy * vD[1]
        + metric->adm_gupxz * vD[2];

  // v^{y} = gamma^{yj}v_{j} = gamma^{jy}v_{j}
  vU[1] = metric->adm_gupxy * vD[0]
        + metric->adm_gupyy * vD[1]
        + metric->adm_gupyz * vD[2];

  // v^{z} = gamma^{zj}v_{j} = gamma^{jz}v_{j}
  vU[2] = metric->adm_gupxz * vD[0]
        + metric->adm_gupyz * vD[1]
        + metric->adm_gupzz * vD[2];
}

/*
 * Function   : compute_rho_W_from_x_and_conservatives
 * Author     : Leo Werneck
 *
 * Computes rho and W from x and the conservative variables using Eq. (42)
 * of https://arxiv.org/pdf/1712.07538.pdf and rho = D/W.
 *
 * Parameters : x        - Pointer to fparams_struct.
 *            : fparams  - Array containing v_{i}.
 *            : rho_ptr  - Stores the output of rho.
 *            : W_ptr    - Stores the output of W.
 *
 * Returns    : Nothing.
 */
static inline void
compute_rho_W_from_x_and_conservatives(
  const double x,
  const fparams_struct *restrict fparams,
  double *restrict rho_ptr,
  double *restrict W_ptr ) {

  // Step 1: Unpack the fparams struct
  const double r = fparams->r;
  const double s = fparams->s;
  const double t = fparams->t;

  // Step 2: Compute W
  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t  )/ (x*x*(x+s)*(x+s));
  Wminus2        = MIN(MAX(Wminus2, fparams->eos->inv_W_max_squared ), 1.0);
  const double W = pow(Wminus2, -0.5);

  // Step 3: Compute rho
  const double rho = fparams->cons_undens->rho/W;

  // Step 4: Set the output
  *rho_ptr = rho;
  *W_ptr   = W;
}

/*
 * Function : froot
 * Author   : Leo Werneck
 *
 * Computes Eq. (33) of Siegel et al., 2018 (arXiv: 1712.07538). The function
 * arguments follow the standards set by the roots.h file.
 *
 * Parameters : x        - The point at which f(x) is evaluated.
 *            : fparams  - Pointer to parameter structed containing auxiliary
 *                         variables needed by this function (see definition
 *                         above).
 *
 * Returns    : Nothing.
 */
static inline double
froot(
      const double x,
      void *restrict fparams ) {

  double rho, P, eps, T, W;
  ((fparams_struct *)fparams)->compute_rho_P_eps_T_W(
    x, fparams, &rho, &P, &eps, &T, &W);

  // Eq: (33) of https://arxiv.org/pdf/1712.07538.pdf
  return x - (1.0 + eps + P/rho)*W;
}

/*
 * Function : compute_BU_SU_Bsq_Ssq_BdotS
 * Author   : Leo Werneck
 *
 * Computes B^{i}, S^{i}, B^2, S^2, B.S = B^{i}S_{i}.
 *
 * Parameters : metric       - Metric quantities
 *            : cons_undens  - Undensitized conservatives
 *            : prims        - Input primitives (for B^{i})
 *            : BU           - Stores B^{i}.
 *            : SU           - Stores S^{i}.
 *            : Bsq          - Stores B^2.
 *            : Ssq          - Stores S^2.
 *            : BdotS        - Stores B.S = B^{i}S_{i}.
 *
 * Returns    : Nothing.
 */
static inline void
compute_BU_SU_Bsq_Ssq_BdotS(
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      const primitive_quantities *restrict prims,
      double *restrict BU,
      double *restrict SU,
      double *restrict Bsq,
      double *restrict Ssq,
      double *restrict BdotS ) {

  // Step 1: Compute S^{2} = gamma^{ij}S_{i}S_{j}
  double SD[3] = {cons_undens->S_x, cons_undens->S_y, cons_undens->S_z};
  double S_squared = compute_S_squared(metric, SD);

  // Step 2: Enforce ceiling on S^{2} (Eq. A5 of [1])
  // Step 2.1: Compute maximum allowed value for S^{2}
  const double S_squared_max = SQR(cons_undens->tau + cons_undens->rho);
  if( S_squared > S_squared_max ) {
    // Step 2.2: Rescale S_{i}
    const double rescale_factor = sqrt(0.9999*S_squared_max/S_squared);
    for(int i=0;i<3;i++)
      SD[i] *= rescale_factor;

    // Step 2.3: Recompute S^{2}
    S_squared = compute_S_squared(metric, SD);
  }
  *Ssq = S_squared;

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  BU[0] = prims->Bx * ONE_OVER_SQRT_4PI;
  BU[1] = prims->By * ONE_OVER_SQRT_4PI;
  BU[2] = prims->Bz * ONE_OVER_SQRT_4PI;
  *Bsq = compute_Bsq_from_Bup(metric, BU);

  // Step 4: Compute B.S = B^{i}S_{i}
  *BdotS = 0.0;
  for(int i=0;i<3;i++) *BdotS += BU[i]*SD[i];

  // Step 5: Compute S^{i}
  raise_vector_3d(metric, SD, SU);
}

/*
 * Function prototypes
 */
int Tabulated_Palenzuela1D(
      void compute_rho_P_eps_T_W(
            const double x,
            fparams_struct *restrict fparams,
            double *restrict rho_ptr,
            double *restrict P_ptr,
            double *restrict eps_ptr,
            double *restrict T_ptr,
            double *restrict W_ptr ),
      const GRHayL_parameters *restrict grhayl_params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics );

int Tabulated_Newman1D(
      void compute_rho_P_eps_T_W(
            const double x,
            fparams_struct *restrict fparams,
            double *restrict rho_ptr,
            double *restrict P_ptr,
            double *restrict eps_ptr,
            double *restrict T_ptr,
            double *restrict W_ptr ),
      const GRHayL_parameters *restrict grhayl_params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics );

#endif // PALENZUELA_H_