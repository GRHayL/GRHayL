#include "ghl_nrpyeos_tabulated.h"

#define munu_index(ir, it, iy) \
  (NRPyEOS_munu_key            \
   + NRPyEOS_ntablekeys * ((ir) + eos->N_rho * ((it) + eos->N_T * (iy))))

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
      const double x,
      int *restrict ix) {

  const double dx = x_arr[1] - x_arr[0];
  int idx = (int)((x - x_arr[0]) / dx);

  if(idx < 0) {
    idx = 0;
  }
  else if(idx > nx - 1) {
    idx = nx - 1;
  }

  *ix = idx;
  return ghl_success;
}

static int find_left_index_bisection(
      const int nx,
      const double *restrict x_arr,
      const double x,
      int *restrict ix) {

  *ix = -1;

  int ia = 0;
  double a = x_arr[ia] - x;
  if(a == 0) {
    *ix = ia;
    return ghl_success;
  }

  int ib = nx - 1;
  double b = x_arr[ib] - x;
  if(b == 0) {
    *ix = ib;
    return ghl_success;
  }

  if(a * b >= 0) {
    return ghl_error_root_not_bracketed;
  }

  do {
    int ic = (ia + ib) / 2;
    double c = x_arr[ic] - x;

    if(c == 0) {
      *ix = ic;
      return ghl_success;
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
    return ghl_error_table_bisection;
  }

  *ix = ia;
  return ghl_success;
}

// This is a simple linear interpolation (linterp) function.
ghl_error_codes_t linterp(
      const int nx,
      const double *restrict x_arr,
      const double *restrict y_arr,
      const double x,
      double *restrict y,
      int (*find_left_index)(const int, const double *restrict, const double, int *restrict)) {

  if(x < x_arr[0] || x > x_arr[nx - 1]) {
    return ghl_error_root_not_bracketed;
  }

  // Set up basic quantities for the interpolation
  int i0 = -1;
  ghl_error_codes_t err = find_left_index(nx, x_arr, x, &i0);
  if(err != ghl_success) {
    return err;
  }
  if(i0 == nx - 1) {
    i0--;
  }
  const int i1 = i0 + 1;
  const double x0 = x_arr[i0];
  const double x1 = x_arr[i1];
  const double y0 = y_arr[i0];
  const double y1 = y_arr[i1];

  // Perform the linear interpolation
  *y = (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
  return ghl_success;
}

// This is a simple discrete derivative.
ghl_error_codes_t discrete_derivative(
      const int nx,
      const double *restrict x_arr,
      const double *restrict y_arr,
      const double x,
      double *restrict dydx,
      int (*find_left_index)(const int, const double *restrict, const double, int *restrict)) {

  if(x < x_arr[0] || x > x_arr[nx - 1]) {
    return ghl_error_root_not_bracketed;
  }

  // Set up basic quantities for the interpolation
  int i0 = -1;
  ghl_error_codes_t err = find_left_index(nx, x_arr, x, &i0);
  if(err != ghl_success) {
    return err;
  }
  if(i0 == nx - 1) {
    i0--;
  }
  const int i1 = i0 + 1;
  const double x0 = x_arr[i0];
  const double x1 = x_arr[i1];
  const double y0 = y_arr[i0];
  const double y1 = y_arr[i1];

  // Perform the derivative
  *dydx = (y1 - y0) / (x1 - x0);
  return ghl_success;
}

// Compute Ye(logrho) using linear interpolation
ghl_error_codes_t NRPyEOS_tabulated_compute_Ye_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      double *restrict Ye) {

  *Ye = NAN;
  ghl_error_codes_t err = linterp(
                          eos->N_rho, eos->table_logrho, eos->Ye_of_lr, log(rho),
                          Ye, find_left_index_uniform_array);

  return err;
}

// Interpolate logP(logrho) using linear interpolation
ghl_error_codes_t NRPyEOS_tabulated_compute_P_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      double *restrict P) {

  *P = NAN;
  double logP = NAN;
  ghl_error_codes_t err = linterp(
                          eos->N_rho, eos->table_logrho, eos->lp_of_lr, log(rho),
                          &logP, find_left_index_uniform_array);
  *P = exp(logP);
  return err;
}

// Interpolate logrho(logP) using linear interpolation
ghl_error_codes_t NRPyEOS_tabulated_compute_rho_from_P(
      const ghl_eos_parameters *restrict eos,
      const double P,
      double *restrict rho) {

  *rho = NAN;
  double logrho = NAN;
  ghl_error_codes_t err = linterp(
                          eos->N_rho, eos->lp_of_lr, eos->table_logrho, log(P),
                          &logrho, find_left_index_bisection);
  *rho = exp(logrho);
  return err;
}

// Interpolate logeps(logrho) using linear interpolation
ghl_error_codes_t NRPyEOS_tabulated_compute_eps_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      double *restrict eps) {

  *eps = NAN;
  double logeps_shifted = NAN;
  ghl_error_codes_t err = linterp(
                          eos->N_rho, eos->table_logrho, eos->le_of_lr, log(rho),
                          &logeps_shifted, find_left_index_bisection);
  *eps = exp(logeps_shifted) - eos->energy_shift;
  return err;
}

ghl_error_codes_t NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos) {

  const int it = ghl_tabulated_get_index_T(eos, T);
  const int nr = eos->N_rho;
  const int ny = eos->N_Ye;

  if(eos->Ye_of_lr == NULL) {
    eos->Ye_of_lr = (double *)malloc(sizeof(double) * nr);
    if(eos->Ye_of_lr == NULL) {
      return ghl_error_out_of_memory;
    }
  }
  double *munu_of_Ye = (double *)malloc(sizeof(double) * ny);
  if(munu_of_Ye == NULL) {
    return ghl_error_out_of_memory;
  }

  for(int ir = 0; ir < nr; ir++) {
    for(int iy = 0; iy < ny; iy++) {
      munu_of_Ye[iy] = eos->table_all[munu_index(ir, it, iy)];
    }
    eos->Ye_of_lr[ir] = find_Ye_st_munu_is_zero(ny, eos->table_Y_e, munu_of_Ye);
  }
  free(munu_of_Ye);

  return ghl_success;
}

ghl_error_codes_t NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos) {

  // Start by obtaining Ye(logrho)
  ghl_error_codes_t err = NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(T, eos);
  if(err != ghl_success) {
    return err;
  }

  // Now allocate memory for logP(logrho) and logeps(logrho)
  double *lp_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  double *le_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  double *lh_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  if(lp_of_lr == NULL || le_of_lr == NULL || lh_of_lr == NULL) {
    free(lp_of_lr);
    free(le_of_lr);
    free(lh_of_lr);
    return ghl_error_out_of_memory;
  }

  // Compute logP(logrho) and logeps(logrho) and logh(logrho)
  for(int ir = 0; ir < eos->N_rho; ir++) {
    const double rho = exp(eos->table_logrho[ir]);
    const double Y_e = eos->Ye_of_lr[ir];
    double P, eps;
    err = ghl_tabulated_compute_P_eps_from_T(eos, rho, Y_e, T, &P, &eps);
    if(err != ghl_success) {
      free(lp_of_lr);
      free(le_of_lr);
      free(lh_of_lr);
      return err;
    }
    lp_of_lr[ir] = log(P);
    le_of_lr[ir] = log(eps + eos->energy_shift);
    lh_of_lr[ir] = log(1.0 + eps + P / rho);
  }

  free(eos->lp_of_lr);
  free(eos->le_of_lr);
  free(eos->lh_of_lr);
  eos->lp_of_lr = lp_of_lr;
  eos->le_of_lr = le_of_lr;
  eos->lh_of_lr = lh_of_lr;

  return ghl_success;
}

void NRPyEOS_tabulated_free_beq_quantities(ghl_eos_parameters *restrict eos) {
  free(eos->Ye_of_lr);
  free(eos->lp_of_lr);
  free(eos->le_of_lr);
  free(eos->lh_of_lr);

  eos->Ye_of_lr = NULL;
  eos->lp_of_lr = NULL;
  eos->le_of_lr = NULL;
  eos->lh_of_lr = NULL;
}

ghl_error_codes_t NRPyEOS_tabulated_compute_dP_drho_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      double *restrict dP_drho) {

  double P = NAN;
  double dlogP_dlogrho = NAN;

  ghl_error_codes_t err = ghl_success;

  err = NRPyEOS_tabulated_compute_P_from_rho(eos, rho, &P);
  if(err != ghl_success) {
    return err;
  }
  err = discrete_derivative(eos->N_rho, eos->table_logrho, eos->lp_of_lr, log(rho),
                            &dlogP_dlogrho, find_left_index_uniform_array);
  if(err != ghl_success) {
    return err;
  }

  *dP_drho = dlogP_dlogrho * P / rho;
  return ghl_success;
}

ghl_error_codes_t NRPyEOS_tabulated_compute_deps_dP_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      double *restrict deps_dP) {

  double P = NAN;
  double eps = NAN;
  double dP_drho = NAN;
  ghl_error_codes_t err = ghl_success;

  err = NRPyEOS_tabulated_compute_P_from_rho(eos, rho, &P);
  if(err != ghl_success) {
    return err;
  }
  err = NRPyEOS_tabulated_compute_eps_from_rho(eos, rho, &eps);
  if(err != ghl_success) {
    return err;
  }
  err = NRPyEOS_tabulated_compute_dP_drho_from_rho(eos, rho, &dP_drho);
  if(err != ghl_success) {
    return err;
  }

  const double energy = rho*(1.+eps);
  const double rhoh = energy + P;

  *deps_dP = rhoh / (dP_drho*rho);
  return ghl_success;
}
