#include "nrpyeos_tabulated.h"
#include "../../Con2Prim/roots.h"

#define munu_index(ir, it, iy) \
  (NRPyEOS_munu_key            \
   + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

#define mu_p_index(ir, it, iy) \
  (NRPyEOS_mu_p_key            \
   + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))
   
#define mu_n_index(ir, it, iy) \
  (NRPyEOS_mu_n_key            \
  + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

#define entropy_index(ir, it, iy) \
  (NRPyEOS_entropy_key         \
    + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

typedef struct {
  double target_entropy;
  int ir;
  double T;   // <--- new field for temperature
  ghl_eos_parameters *eos;
} brent_params;     

static double find_Ye_st_munu_is_zero(
      const int n,
      const double *restrict Ye,
      const double *restrict munu) {

  int i0 = -1;
  for(int i = 0; i < n - 1; i++) {
    if(munu[i] * munu[i + 1] < 0) {
      i0 = i;
      break;
    }
  }

  // If munu does not cross zero, return the minimum Ye
  if(i0 == -1) {
    return Ye[0];
  }

  const int i1 = i0 + 1;
  const double x0 = Ye[i0];
  const double x1 = Ye[i1];
  const double y0 = munu[i0];
  const double y1 = munu[i1];

  return (y0 * x1 - y1 * x0) / (y0 - y1);
}

static int find_left_index_uniform_array(
      const int nx,
      const double *restrict x_arr,
      const double x) {

  return (x - x_arr[0]) / (x_arr[1] - x_arr[0]) + 0.5;
}

static int
find_left_index_bisection(const int nx, const double *restrict x_arr, const double x) {

  int ia = 0;
  double a = x_arr[ia] - x;
  if(a == 0) {
    return ia;
  }

  int ib = nx - 1;
  double b = x_arr[ib] - x;
  if(b == 0) {
    return ib;
  }

  if(a * b >= 0) {
    ghl_error("Interval [%g, %g] does not bracket the root %g\n", a, b, x);
  }

  do {
    int ic = (ia + ib) / 2;
    double c = x_arr[ic] - x;

    if(c == 0) {
      return ic;
    }
    else if(b * c < 0) {
      ia = ic;
      a = c;
    }
    else {
      ib = ic;
      b = c;
    }
  } while(ib - ia > 1);

  if(a > 0 || b < 0) {
    ghl_error("Bisection failed: %g not in [%g, %g]\n", x, a, b);
  }
  return ia;
}

// This is a simple linear interpolation (linterp) function.
static double linterp(
      const int nx,
      const double *restrict x_arr,
      const double *restrict y_arr,
      const double x,
      int (*find_left_index)(const int, const double *restrict, const double)) {

  if(x < x_arr[0] || x > x_arr[nx - 1]) {
    ghl_error("Point (%e) out of array bounds [%e, %e]\n", x, x_arr[0], x_arr[nx - 1]);
  }

  // Set up basic quantities for the interpolation
  const int i0 = find_left_index(nx, x_arr, x);
  const int i1 = i0 + 1;
  const double x0 = x_arr[i0];
  const double x1 = x_arr[i1];
  const double y0 = y_arr[i0];
  const double y1 = y_arr[i1];

  // Perform the linear interpolation
  return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
}

// Compute Ye(logrho) using linear interpolation
double NRPyEOS_tabulated_compute_Ye_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho) {

  return linterp(
        eos->N_rho, eos->table_logrho, eos->Ye_of_lr, log(rho),
        find_left_index_uniform_array);
}

// Interpolate logP(logrho) using linear interpolation
double NRPyEOS_tabulated_compute_P_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho) {

  return exp(linterp(
        eos->N_rho, eos->table_logrho, eos->lp_of_lr, log(rho),
        find_left_index_uniform_array));
}

// Interpolate logrho(logP) using linear interpolation
double NRPyEOS_tabulated_compute_rho_from_P(
      const ghl_eos_parameters *restrict eos,
      const double P) {

  return exp(linterp(
        eos->N_rho, eos->lp_of_lr, eos->table_logrho, log(P),
        find_left_index_bisection));
}

// Interpolate logeps(logrho) using linear interpolation
double NRPyEOS_tabulated_compute_eps_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho) {

  return exp(linterp(
               eos->N_rho, eos->table_logrho, eos->le_of_lr, log(rho),
               find_left_index_bisection))
         - eos->energy_shift;
}

void NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos) {

  const int it = ghl_tabulated_get_index_T(eos, T);
  const int nr = eos->N_rho;
  const int ny = eos->N_Ye;

  if(eos->Ye_of_lr == NULL) {
    eos->Ye_of_lr = (double *)malloc(sizeof(double) * nr);
  }
  double *munu_of_Ye = (double *)malloc(sizeof(double) * ny);

  for(int ir = 0; ir < nr; ir++) {
    for(int iy = 0; iy < ny; iy++) {
      munu_of_Ye[iy] = eos->table_all[munu_index(ir, it, iy)];
    }
    eos->Ye_of_lr[ir] = find_Ye_st_munu_is_zero(ny, eos->table_Y_e, munu_of_Ye);
  }
  free(munu_of_Ye);
}

void NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos) {

  // Start by obtaining Ye(logrho)
  NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(T, eos);

  // Now allocate memory for logP(logrho) and logeps(logrho)
  if(eos->lp_of_lr == NULL) {
    eos->lp_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  }
  if(eos->le_of_lr == NULL) {
    eos->le_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  }

  // Compute logP(logrho) and logeps(logrho)
  for(int ir = 0; ir < eos->N_rho; ir++) {
    const double rho = exp(eos->table_logrho[ir]);
    const double Y_e = eos->Ye_of_lr[ir];
    double P, eps;
    ghl_tabulated_compute_P_eps_from_T(eos, rho, Y_e, T, &P, &eps);
    eos->lp_of_lr[ir] = log(P);
    eos->le_of_lr[ir] = log(eps + eos->energy_shift);
  }
}

static double entropy_at_Ye_T(
    double Ye,
    double T,
    const int ir,
    ghl_eos_parameters *restrict eos) {

  Ye = MAX(Ye, eos->Y_e_min);
  Ye = MIN(Ye, eos->Y_e_max);

  T = MAX(T, eos->T_min);
  T = MIN(T, eos->T_max);

  const int nt = eos->N_T;
  const int ny = eos->N_Ye;

  const int it = find_left_index_bisection(nt, eos->table_logT, log(T));

  if (it < 0 || it >= nt - 1) {
    ghl_error("Error: Temperature index out of bounds.\n");
    return NAN; // Indicate an error
  }

  double *entropy_of_Ye = (double *)malloc(sizeof(double) * ny);
  if (entropy_of_Ye == NULL) {
    ghl_error("Memory allocation failed.\n");
    return NAN; // Indicate an error
  }

  for (int iy = 0; iy < ny; iy++) {
    entropy_of_Ye[iy] = eos->table_all[entropy_index(ir, it, iy)];
  }
  
  double entropy_val = linterp(ny, eos->table_Y_e, entropy_of_Ye, Ye, find_left_index_uniform_array);
  free(entropy_of_Ye);
  return entropy_val;
}

static double mu_p__minus__mu_n__at_Ye_T(
    double Ye,
    double T,
    const int ir,
    ghl_eos_parameters *restrict eos) {

  Ye = MAX(Ye, eos->Y_e_min);
  Ye = MIN(Ye, eos->Y_e_max);

  T = MAX(T, eos->T_min);
  T = MIN(T, eos->T_max);

  const int nt = eos->N_T;
  const int ny = eos->N_Ye;

  const int it = find_left_index_bisection(nt, eos->table_logT, log(T));

  if (it < 0 || it >= nt - 1) {
    ghl_error("Error: Temperature index out of bounds.\n");
    return NAN; // Indicate an error
  }

  double *mu_p__minus__mu_n__of_Ye = (double *)malloc(sizeof(double) * ny);
  if (mu_p__minus__mu_n__of_Ye == NULL) {
    ghl_error("Memory allocation failed.\n");
    return NAN; // Indicate an error
  }

  for (int iy = 0; iy < ny; iy++) {
    // munu_of_Ye[iy] = eos->table_all[munu_index(ir, it, iy)];
    // The above is commented out because this is for beta-equilibrium for a cold, neutrino-less star.
    // Given constant entropy, with a temperature gradient, we aim to model a newly formed, hot proto-neutron star.
    // Thus, we need beta-equilibrium, with mu_n + mu_nu = mu_p + mu_e.
    // --> mu_p - mu_n = 0.
    mu_p__minus__mu_n__of_Ye[iy] = eos->table_all[mu_p_index(ir, it, iy)] - eos->table_all[mu_n_index(ir, it, iy)];
  }

  double mu_p__minus__mu_n_val = linterp(ny, eos->table_Y_e, mu_p__minus__mu_n__of_Ye, Ye, find_left_index_uniform_array);
  free(mu_p__minus__mu_n__of_Ye);
  return mu_p__minus__mu_n_val;
}

// Add a wrapper for mu_p__minus__mu_n__at_Ye_T so that it matches the prototype expected by ghl_brent.
static double mu_p__minus__mu_n_wrapper(const double Ye,
                           const ghl_parameters *restrict params,
                           const ghl_eos_parameters *restrict eos,
                           const ghl_conservative_quantities *restrict cons_undens,
                           fparams_struct *restrict fparams,
                           ghl_primitive_quantities *restrict prims) {
  // Cast back fparams to our brent_params struct.
  brent_params *bparams = (brent_params *)fparams;
  return mu_p__minus__mu_n__at_Ye_T(Ye, bparams->T, bparams->ir, bparams->eos);
}

static double brent_objective(
  const double T,
  const ghl_parameters *restrict params,
  const ghl_eos_parameters *restrict eos,
  const ghl_conservative_quantities *restrict cons_undens,
  fparams_struct *restrict fparams,
  ghl_primitive_quantities *restrict prims)
{
  // Cast the fparams pointer to retrieve our brent_params structure
  brent_params *bparams = (brent_params *)fparams;
  bparams->T = T;

  double ye_beq = eos->Y_e_min;
  roots_params root_params = {0};
  root_params.tol = 1e-10; // Adjust tolerance as needed
  root_params.max_iters = 100;

  // Use the wrapper function for mu_p__minus__mu_n_wrapper solve
  ghl_brent(mu_p__minus__mu_n_wrapper, params, eos, cons_undens, fparams, prims,
            bparams->eos->table_Y_e[1],
            bparams->eos->table_Y_e[bparams->eos->N_Ye - 1],
            &root_params);
  ye_beq = root_params.root;

  double entropy_beq = entropy_at_Ye_T(ye_beq, T, bparams->ir, bparams->eos);
  return entropy_beq - bparams->target_entropy;
}

void NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_entropy(
    const double entropy,
    ghl_eos_parameters *restrict eos) {

  const int nr = eos->N_rho;

  if (eos->Ye_of_lr == NULL) {
    eos->Ye_of_lr = (double *)malloc(sizeof(double) * nr);
  }
  if (eos->lp_of_lr == NULL) {
    eos->lp_of_lr = (double *)malloc(sizeof(double) * nr);
  }
  if (eos->le_of_lr == NULL) {
    eos->le_of_lr = (double *)malloc(sizeof(double) * nr);
  }

  roots_params root_params = {0};
  root_params.tol = 1e-10; // Adjust tolerance as needed
  root_params.max_iters = 100;

  brent_params bparams = {0};
  bparams.target_entropy = entropy;
  bparams.eos = eos;

  double (*entropy_func)(const double, const ghl_parameters *, const ghl_eos_parameters *, const ghl_conservative_quantities *, fparams_struct *, ghl_primitive_quantities *);

  entropy_func = (double (*)(const double, const ghl_parameters *, const ghl_eos_parameters *, const ghl_conservative_quantities *, fparams_struct *, ghl_primitive_quantities *))brent_objective;

  for (int ir = 0; ir < nr; ir++) {
    bparams.ir = ir;
    ghl_brent(entropy_func, NULL, eos, NULL, (fparams_struct*)&bparams, NULL, eos->table_logT[0], eos->table_logT[eos->N_T - 1], &root_params);
    double temp_beq = root_params.root;

    double ye_beq = eos->Y_e_min;
    // double mu_p__minus__mu_n_wrapper_beq = 0.0;
    roots_params ye_params = {0};
    ye_params.tol = 1e-8;
    ye_params.max_iters = 100;

    double (*ye_func)(const double, const ghl_parameters *, const ghl_eos_parameters *,
      const ghl_conservative_quantities *, fparams_struct *,
      ghl_primitive_quantities *);

    ye_func = mu_p__minus__mu_n_wrapper;
    ghl_brent(ye_func, NULL, eos, NULL, (fparams_struct*)&bparams, NULL,
              eos->table_Y_e[1], eos->table_Y_e[eos->N_Ye-1], &ye_params);

    ye_beq = ye_params.root;

    eos->Ye_of_lr[ir] = ye_beq;

    double P, eps;
    ghl_tabulated_compute_P_eps_from_T(eos, exp(eos->table_logrho[ir]), ye_beq, temp_beq, &P, &eps);
    eos->lp_of_lr[ir] = log(P);
    eos->le_of_lr[ir] = log(eps + eos->energy_shift);
  }
}

void NRPyEOS_tabulated_free_beq_quantities(ghl_eos_parameters *restrict eos) {
  if (eos->Ye_of_lr) {
    free(eos->Ye_of_lr);
    eos->Ye_of_lr = NULL;
  }
  if (eos->lp_of_lr) {
    free(eos->lp_of_lr);
    eos->lp_of_lr = NULL;
  }
  if (eos->le_of_lr) {
    free(eos->le_of_lr);
    eos->le_of_lr = NULL;
  }
}

