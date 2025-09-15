#include "ghl_nrpyeos_tabulated.h"

// Set up the criterion for whether two elements are the same
static bool comp(const double a, const double b) {
  if(a == 0 || b == 0) {
    return (fabs(a - b) < 1e-15);
  }
  return (fabs(1 - a / b) < 1e-15);
}

// For a given array x_arr and an input value x, this function returns the
// index i such that x in (x_arr[i],  x_arr[i+1]), unless x is equal to one
// of the elements in x_arr, in which case the element index is returned.
static int bisect_left(
    const int n,
    const double *restrict x_arr,
    const double x,
    bool (*comp)(const double, const double)) {

  int ia = 0;
  int ib = n - 1;
  double a = x_arr[ia];
  double b = x_arr[ib];
  if(x < a || x > b) {
    ghl_error("Value %g is out of array bounds [%g, %g]\n", x, a, b);
  }

  if(comp(a, x)) {
    return ia;
  }
  if(comp(b, x)) {
    return ib;
  }

  while(ib - ia > 1) {
    int ic = (ia + ib) / 2;
    double c = x_arr[ic];
    if(comp(c, x)) {
      return ic;
    }
    else if(c < x) {
      ia = ic;
      a = c;
    }
    else {
      ib = ic;
      b = c;
    }
  }
  return ia;
}

int NRPyEOS_tabulated_get_index_T(
    const ghl_eos_parameters *restrict eos,
    const double T) {
  return bisect_left(eos->N_T, eos->table_logT, log(T), comp);
}

// Uncomment these functions if we ever need them
// int NRPyEOS_tabulated_get_index_rho(
//     const ghl_eos_parameters *restrict eos,
//     const double rho) {
//   return bisect_left(eos->N_rho, eos->table_logrho, log(rho), comp);
// }
//
// int NRPyEOS_tabulated_get_index_Ye(
//     const ghl_eos_parameters *restrict eos,
//     const double Y_e) {
//   return bisect_left(eos->N_Ye, eos->table_Y_e, Y_e, comp);
// }
