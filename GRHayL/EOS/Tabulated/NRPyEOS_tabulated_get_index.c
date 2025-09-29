#include "ghl_nrpyeos_tabulated.h"

// Set up the criterion for whether two elements are the same
GHL_DEVICE
static bool comp(const double a, const double b) {
  if(a == 0 || b == 0) {
    return (fabs(a - b) < 1e-15);
  }
  return (fabs(1 - a / b) < 1e-15);
}

// For a given array x_arr and an input value x, this function returns the
// index i such that x in (x_arr[i],  x_arr[i+1]), unless x is equal to one
// of the elements in x_arr, in which case the element index is returned.
GHL_DEVICE
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
    return ghl_error_root_not_bracketed;
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

GHL_DEVICE
int NRPyEOS_tabulated_get_index_T(
    const ghl_eos_parameters *restrict eos,
    const double T) {
  return bisect_left(eos->N_T, eos->table_logT, log(T), comp);
}
