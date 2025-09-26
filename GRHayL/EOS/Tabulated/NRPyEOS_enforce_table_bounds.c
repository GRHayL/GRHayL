#include "ghl.h"
GHL_DEVICE
void NRPyEOS_enforce_table_bounds_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict T) {

  // Enforce bounds on rho
  *rho = fmin(fmax(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = fmin(fmax(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on T
  *T = fmin(fmax(*T, eos->T_min), eos->T_max);
}

void NRPyEOS_enforce_table_bounds_rho_Ye_eps(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict eps) {

  // Enforce bounds on rho
  *rho = fmin(fmax(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = fmin(fmax(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on eps
  *eps = fmin(fmax(*eps, eos->eps_min), eos->eps_max);
}

void NRPyEOS_enforce_table_bounds_rho_Ye_S(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict S) {

  // Enforce bounds on rho
  *rho = fmin(fmax(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = fmin(fmax(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on S
  *S = fmin(fmax(*S, eos->entropy_min), eos->entropy_max);
}

void NRPyEOS_enforce_table_bounds_rho_Ye_P(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict P) {

  // Enforce bounds on rho
  *rho = fmin(fmax(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = fmin(fmax(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on P
  *P = fmin(fmax(*P, eos->press_min), eos->press_max);
}
