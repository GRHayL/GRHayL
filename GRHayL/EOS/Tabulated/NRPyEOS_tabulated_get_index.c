#include "nrpyeos_tabulated.h"

static int bisect_left(
      const int n,
      const double *restrict x_arr,
      const double x) {

  int ia = 0;
  int ib = n-1;
  double a = x_arr[ia];
  double b = x_arr[ib];
  if(x < a || x > b)
    ghl_error("Value %g is out of array bounds [%g, %g]\n", x, a, b);

  if(x == a) return ia;
  if(x == b) return ib;

  while(true) {
    if(ib - ia == 1)
      break;
    int ic = (ia + ib)/2;
    double c = x_arr[ic];
    if(c == x) {
      printf("Here : ");
      a = c;
      ia = ic;
      break;
    }
    else if(x > c) {
      ia = ic;
      a = c;
    }
    else {
      ib = ic;
      b = c;
    }
  }
  if(fabs(1-x/a) < 1e-14 || fabs(1-x/b) > 1e-14)
    return ia;
  return ib;
}

int NRPyEOS_tabulated_get_index_rho(const ghl_eos_parameters *restrict eos, const double rho) {
  return bisect_left(eos->N_rho, eos->table_logrho, log(rho));
}

int NRPyEOS_tabulated_get_index_T(const ghl_eos_parameters *restrict eos, const double T) {
  return bisect_left(eos->N_T, eos->table_logT, log(T));
}

int NRPyEOS_tabulated_get_index_Ye(const ghl_eos_parameters *restrict eos, const double Y_e) {
  return bisect_left(eos->N_Ye, eos->table_Y_e, Y_e);
}
