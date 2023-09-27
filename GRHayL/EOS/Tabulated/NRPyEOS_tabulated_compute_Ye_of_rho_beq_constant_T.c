#include "nrpyeos_tabulated.h"

#define munu_index(ir, it, iy) ( NRPyEOS_munu_key + NRPyEOS_ntablekeys*((ir) + eos->N_rho*((it) + eos->N_T*(iy))) )

static double
find_Ye_st_munu_is_zero(
      const int n,
      const double *restrict Ye,
      const double *restrict munu ) {

  int i0 = -1;
  for(int i=0;i<n-1;i++) {
    if(munu[i]*munu[i+1] < 0) {
      i0 = i;
      break;
    }
  }

  // If munu does not cross zero, return the minimum Ye
  if(i0 == -1)
    return Ye[0];

  const int    i1 = i0+1;
  const double x0 = Ye[i0];
  const double x1 = Ye[i1];
  const double y0 = munu[i0];
  const double y1 = munu[i1];

  return (y0*x1 - y1*x0)/(y0-y1);
}

// This is a simple linear interpolation (linterp) function.
static double
linterp(
      const int nx,
      const double *restrict x_arr,
      const double *restrict y_arr,
      const double x ) {

  if( x < x_arr[0] || x > x_arr[nx-1] )
    ghl_error("Point (%e) out of array bounds [%e, %e]\n", x, x_arr[0], x_arr[nx-1]);


  // Set up basic quantities for the interpolation
  const double dx = x_arr[1] - x_arr[0];
  const int    i0 = (x - x_arr[0])/dx + 0.5;
  const int    i1 = i0 + 1;
  const double x0 = x_arr[i0];
  const double x1 = x_arr[i1];
  const double y0 = y_arr[i0];
  const double y1 = y_arr[i1];

  // Perform the linear interpolation
  return ( y0*( x1 - x ) + y1*( x - x0 ) )/dx;
}

// Interpolate Ye(logrho) using linear interpolation
double
NRPyEOS_tabulated_compute_Ye_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho ) {

  return linterp(eos->N_rho, eos->table_logrho, eos->Ye_of_lr, log(rho));
}

// Interpolate logP(logrho) using linear interpolation
double
NRPyEOS_tabulated_compute_P_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho ) {

  return exp(linterp(eos->N_rho, eos->table_logrho, eos->lp_of_lr, log(rho)));
}

// Interpolate logrho(logP) using linear interpolation
double
NRPyEOS_tabulated_compute_rho_from_P(
      const ghl_eos_parameters *restrict eos,
      const double P ) {

  return exp(linterp(eos->N_rho, eos->lp_of_lr, eos->table_logrho, log(P)));
}

// Interpolate logeps(logrho) using linear interpolation
double
NRPyEOS_tabulated_compute_eps_from_rho(
      const ghl_eos_parameters *restrict eos,
      const double rho ) {

  return exp(linterp(eos->N_rho, eos->table_logrho, eos->le_of_lr, log(rho))) - eos->energy_shift;
}

void
NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos ) {

  const int it = ghl_tabulated_get_index_T(eos, T);
  const int nr = eos->N_rho;
  const int ny = eos->N_Ye;

  eos->Ye_of_lr      = (double *)malloc(sizeof(double)*nr);
  double *munu_of_Ye = (double *)malloc(sizeof(double)*ny);

  for(int ir=0;ir<nr;ir++) {
    for(int iy=0;iy<ny;iy++)
      munu_of_Ye[iy] = eos->table_all[munu_index(ir, it, iy)];
    eos->Ye_of_lr[ir] = find_Ye_st_munu_is_zero(ny, eos->table_Y_e, munu_of_Ye);
  }

  free(munu_of_Ye);
}

void
NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(
      const double T,
      ghl_eos_parameters *restrict eos ) {

  // Start by obtaining Ye(logrho)
  NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(T, eos);

  // Now allocate memory for logP(logrho) and logeps(logrho)
  eos->lp_of_lr = (double *)malloc(sizeof(double)*eos->N_rho);
  eos->le_of_lr = (double *)malloc(sizeof(double)*eos->N_rho);

  // Compute logP(logrho) and logeps(logrho)
  for(int ir=0;ir<eos->N_rho;ir++) {
    const double rho = exp(eos->table_logrho[ir]);
    const double Y_e = eos->Ye_of_lr[ir];
    double P, eps;
    ghl_tabulated_compute_P_eps_from_T(eos, rho, Y_e, T, &P, &eps);
    eos->lp_of_lr[ir] = log(P);
    eos->le_of_lr[ir] = log(eps + eos->energy_shift);
  }
}

void
NRPyEOS_tabulated_free_beq_quantities( ghl_eos_parameters *restrict eos ) {
  if( eos->Ye_of_lr ) { free(eos->Ye_of_lr); eos->Ye_of_lr = NULL; }
  if( eos->lp_of_lr ) { free(eos->lp_of_lr); eos->lp_of_lr = NULL; }
  if( eos->le_of_lr ) { free(eos->le_of_lr); eos->le_of_lr = NULL; }
}
