#include "con2prim.h"

double
zbrent(
      double func(const eos_parameters *restrict,
                  const conservative_quantities *restrict,
                  const double *restrict,
                  const bool,
                  const double,
                  primitive_quantities *restrict,
                  double *restrict),
      const eos_parameters *restrict eos,
      const conservative_quantities *restrict cons_undens,
      const double *restrict params,
      const bool evolve_T,
      primitive_quantities *restrict prims,
      double *restrict temp_guess,
      double x1,
      double x2,
      double tol_x,
      con2prim_diagnostics *restrict stats );

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

static inline void
raise_vector_3d(
      const metric_quantities *restrict metric,
      const double *restrict SD,
      double *restrict SU ) {

  // S^{x} = gamma^{xj}S_{j} = gamma^{jx}S_{j}
  SU[0] = metric->adm_gupxx * SD[0]
        + metric->adm_gupxy * SD[1]
        + metric->adm_gupxz * SD[2];

  // S^{y} = gamma^{yj}S_{j} = gamma^{jy}S_{j}
  SU[1] = metric->adm_gupxy * SD[0]
        + metric->adm_gupyy * SD[1]
        + metric->adm_gupyz * SD[2];

  // S^{z} = gamma^{zj}S_{j} = gamma^{jz}S_{j}
  SU[2] = metric->adm_gupxz * SD[0]
        + metric->adm_gupyz * SD[1]
        + metric->adm_gupzz * SD[2];
}

void
compute_rho_Ye_T_P_eps_W_from_x_and_conservatives(
  const eos_parameters *restrict eos,
  const conservative_quantities *restrict cons_undens,
  const double *restrict params,
  const bool evolve_T,
  const double x,
  const double *temp_guess,
  double *restrict rho_ptr,
  double *restrict Y_e_ptr,
  double *restrict T_ptr,
  double *restrict P_ptr,
  double *restrict eps_ptr,
  double *restrict W_ptr ) {

  // computes f(x) from x and q,r,s,t

  const double q = params[0];
  const double r = params[1];
  const double s = params[2];
  const double t = params[3];

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t  )/ (x*x*(x+s)*(x+s));
  Wminus2        = fmin(fmax(Wminus2, eos->inv_W_max_squared ), 1.0);
  const double W = pow(Wminus2, -0.5);

  double rho     = cons_undens->rho/W;
  double ye      = cons_undens->Y_e/cons_undens->rho;
  double temp    = *temp_guess;
  double P       = 0.0;
  double eps     = 0.0;

  if( evolve_T ) {
    eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W) );
    eos->tabulated_compute_P_T_from_eps( eos, rho, ye, eps, &P, &temp );
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

double
func_root(
      const eos_parameters *restrict eos,
      const conservative_quantities *restrict cons_undens,
      const double *restrict params,
      const bool evolve_T,
      const double x,
      primitive_quantities *restrict prims,
      double *restrict temp_guess ) {

  double rho, Y_e, T, P, eps, W;
  compute_rho_Ye_T_P_eps_W_from_x_and_conservatives(
    eos, cons_undens, params, evolve_T, x, temp_guess,
    &rho, &Y_e, &T, &P, &eps, &W );

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
int Tabulated_Palenzuela1D(
      const GRHayL_parameters *restrict grhayl_params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims_guess,
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
  const double Bup[3] = {prims_guess->Bx * ONE_OVER_SQRT_4PI,
                         prims_guess->By * ONE_OVER_SQRT_4PI,
                         prims_guess->Bz * ONE_OVER_SQRT_4PI};
  const double B_squared = compute_Bsq_from_Bup(metric, Bup);

  // Step 4: Compute B.S = B^{i}S_{i}
  double BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += Bup[i]*SD[i];

  // Step 5: Set specific quantities for this routine (Eq. A7 of [1])
  const double invD  = 1.0/cons_undens->rho;
  const double q     = cons_undens->tau * invD;
  const double r     = S_squared * invD * invD;
  const double s     = B_squared * invD;
  const double t     = BdotS/pow(cons_undens->rho, 1.5);

  // Step 6: Bracket x (Eq. A8 of [1])
  double xlow = 1 + q - s;
  double xup  = 2*(1 + q) - s;

  // Step 7: Set initial guess for temperature
  double temp_guess = prims_guess->temperature;

  // Step 8: Call the main function and perform the con2prim
  const double params[4] = {q, r, s, t};
  const double tol_x     = 1e-10;

  double x = zbrent(func_root, eos, cons_undens, params,
                    grhayl_params->evolve_temp, prims_guess,
                    &temp_guess, xlow, xup, tol_x, diagnostics);

  // Step 9: Set core primitives using the EOS and the root
  double W;
  compute_rho_Ye_T_P_eps_W_from_x_and_conservatives(
    eos, cons_undens, params, grhayl_params->evolve_temp, x, &temp_guess,
    &prims_guess->rho, &prims_guess->Y_e, &prims_guess->temperature,
    &prims_guess->press, &prims_guess->eps, &W );

  // Step 10: Compute the velocities using Eq. (24) in [2]. Note, however, that
  //          GRHayL expects the velocity tilde(u)^{i} := W v^{i},
  // Step 10.a: Compute S^{i}
  double SU[3];
  raise_vector_3d(metric, SD, SU);

  // Step 10.b: Set Z
  const double Z = x*prims_guess->rho*W;

  // Step 10.c: Compute tilde(u)^{i}
  prims_guess->vx = W*(SU[0] + BdotS*prims_guess->Bx/Z)/(Z+B_squared);
  prims_guess->vy = W*(SU[1] + BdotS*prims_guess->By/Z)/(Z+B_squared);
  prims_guess->vz = W*(SU[2] + BdotS*prims_guess->Bz/Z)/(Z+B_squared);

  return diagnostics->c2p_failed;
}
