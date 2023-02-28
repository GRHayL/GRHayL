#include "con2prim.h"
#include "roots.h"

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

typedef struct fparams_struct {
  bool evolve_T;
  double temp_guess, q, r, s, t;
  const eos_parameters *eos;
  const conservative_quantities *cons_undens;
} fparams_struct;

static void
compute_rho_Ye_T_P_eps_W_from_x_and_conservatives(
  const double x,
  const fparams_struct *restrict fparams,
  double *restrict rho_ptr,
  double *restrict Y_e_ptr,
  double *restrict T_ptr,
  double *restrict P_ptr,
  double *restrict eps_ptr,
  double *restrict W_ptr ) {

  // Step 1: Unpack the fparams struct
  const double q = fparams->q;
  const double r = fparams->r;
  const double s = fparams->s;
  const double t = fparams->t;
  const eos_parameters *eos = fparams->eos;
  const conservative_quantities *cons_undens = fparams->cons_undens;

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t  )/ (x*x*(x+s)*(x+s));
  Wminus2        = fmin(fmax(Wminus2, eos->inv_W_max_squared ), 1.0);
  const double W = pow(Wminus2, -0.5);

  double rho     = cons_undens->rho/W;
  double ye      = cons_undens->Y_e/cons_undens->rho;
  double temp    = fparams->temp_guess;
  double P       = 0.0;
  double eps     = 0.0;

  if( fparams->evolve_T ) {
    double ent = cons_undens->entropy/W;
    rho = MIN(MAX(rho, eos->table_rho_min), eos->table_rho_max);
    ent = MIN(MAX(ent, eos->table_ent_min), eos->table_ent_max);
    ye  = MIN(MAX(ye , eos->table_Ye_min ), eos->table_Ye_max );
    eos->tabulated_compute_P_eps_T_from_S( eos, rho, ye, ent, &P, &eps, &temp );
  }
  else
    eos->tabulated_compute_P_eps_from_T( eos, rho, ye, temp, &P, &eps );

  // Set the output
  *rho_ptr = rho;
  *Y_e_ptr = ye;
  *T_ptr   = temp;
  *P_ptr   = P;
  *eps_ptr = eps;
  *W_ptr   = W;
}

/*
 * Function : froot
 * Author   : Leo Werneck (adapted from Siegel et al)
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
static double
froot( const double x, void *restrict fparams ) {

  double rho, Y_e, T, P, eps, W;
  compute_rho_Ye_T_P_eps_W_from_x_and_conservatives(
    x, fparams, &rho, &Y_e, &T, &P, &eps, &W );

  //printf("ans = %.15e\n", x - (1.0 + eps + P/rho)*W);

  return x - (1.0 + eps + P/rho)*W;
}

/* Function    : initialize_palenzuela_quantities
 * Description : Initializes the palenzuela_quantities struct.
 *
 * Inputs      : metric         - metric_quantities struct containing local
 *                                metric data (gamma_{ij}, g_{mu nu}, etc).
 *             : prims          - primitive_quantities struct containing
 *                                local hydrodynamic quantities.
 *             : cons           - conservative_quantities struct containing
 *                                local conservative quantities.
 *
 * Outputs     : pal            - palenzuela_quantities struct that will be
 *                                initialized.
 * References  : [1] Palenzuela et al. (2015) arXiv:1505.01607
 *             : [2] Siegel et al. (2018) arXiv:1712.07538
 */
int Tabulated_Palenzuela1D_entropy(
      const GRHayL_parameters *restrict grhayl_params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

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

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  const double Bup[3] = {prims->Bx * ONE_OVER_SQRT_4PI,
                         prims->By * ONE_OVER_SQRT_4PI,
                         prims->Bz * ONE_OVER_SQRT_4PI};
  const double B_squared = compute_Bsq_from_Bup(metric, Bup);

  // Step 4: Compute B.S = B^{i}S_{i}
  double BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += Bup[i]*SD[i];

  // Step 5: Set specific quantities for this routine (Eq. A7 of [1])
  fparams_struct fparams;
  const double invD = 1.0/cons_undens->rho;
  fparams.q = cons_undens->tau * invD;
  fparams.r = S_squared * invD * invD;
  fparams.s = B_squared * invD;
  fparams.t = BdotS/pow(cons_undens->rho, 1.5);
  fparams.eos = eos;
  fparams.cons_undens = cons_undens;

  // printf("q = %.15e\n", fparams.q);
  // printf("r = %.15e\n", fparams.r);
  // printf("s = %.15e\n", fparams.s);
  // printf("t = %.15e\n", fparams.t);

  // Step 6: Bracket x (Eq. A8 of [1])
  double xlow = 1 + fparams.q - fparams.s;
  double xup  = 2*(1 + fparams.q) - fparams.s;

  printf("[%.15e, %.15e]\n", xlow, xup);
  printf("[%.15e, %.15e]\n", froot(xlow, &fparams), froot(xup, &fparams));

  // Step 7: Set initial guess for temperature
  fparams.temp_guess = prims->temperature;

  // Step 8: Call the main function and perform the con2prim
  roots_params rparams;
  rparams.tol = 1e-15;
  rparams.max_iters = 300;
  brent(froot, &fparams, xlow, xup, &rparams);
  if( rparams.error_key != roots_success ) {
    fparams.temp_guess = eos->T_max;
    brent(froot, &fparams, xlow, xup, &rparams);
    if( rparams.error_key != roots_success ) {
      grhayl_warn("Brent's method did not find a root (see error report below)\n");
      roots_info(&rparams);
      return rparams.error_key;
    }
  }

  double x = rparams.root;

  // printf(" x  = %.15e\n", x);
  // printf("res = %.15e\n", bparams.residual);

  // Step 9: Set core primitives using the EOS and the root
  double W;
  compute_rho_Ye_T_P_eps_W_from_x_and_conservatives(
    x, &fparams, &prims->rho, &prims->Y_e, &prims->temperature,
    &prims->press, &prims->eps, &W );

  // Step 10: Compute the velocities using Eq. (24) in [2]. Note, however, that
  //          GRHayL expects the velocity tilde(u)^{i} := W v^{i},
  // Step 10.a: Compute S^{i}
  double SU[3];
  raise_vector_3d(metric, SD, SU);

  // Step 10.b: Set Z
  const double Z = x*prims->rho*W;

  // Step 10.c: Compute tilde(u)^{i}
  prims->vx = W*(SU[0] + BdotS*prims->Bx/Z)/(Z+B_squared);
  prims->vy = W*(SU[1] + BdotS*prims->By/Z)/(Z+B_squared);
  prims->vz = W*(SU[2] + BdotS*prims->Bz/Z)/(Z+B_squared);

  return roots_success;
}
