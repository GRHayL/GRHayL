#include "ghl.h"

/// @brief Enforce table bounds on density, electron fraction, and temperature
/// @param eos Pointer to GRHayL EOS parameters struct
/// @param rho Pointer to density
/// @param Y_e Pointer to electron fraction
/// @param T Pointer to temperature
void NRPyEOS_enforce_table_bounds_rho_Ye_T(
    const ghl_eos_parameters *restrict eos,
    double *restrict rho,
    double *restrict Y_e,
    double *restrict T) {

  // Enforce bounds on rho
  *rho = MIN(MAX(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = MIN(MAX(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on T
  *T = MIN(MAX(*T, eos->T_min), eos->T_max);
}

/// @brief Enforce table bounds on density, electron fraction, and specific internal energy
/// @param eos Pointer to GRHayL EOS parameters struct
/// @param rho Pointer to density
/// @param Y_e Pointer to electron fraction
/// @param eps Pointer to specific internal energy
void NRPyEOS_enforce_table_bounds_rho_Ye_eps(
    const ghl_eos_parameters *restrict eos,
    double *restrict rho,
    double *restrict Y_e,
    double *restrict eps) {

  // Enforce bounds on rho
  *rho = MIN(MAX(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = MIN(MAX(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on eps
  *eps = MIN(MAX(*eps, eos->eps_min), eos->eps_max);
}

/// @brief Enforce table bounds on density, electron fraction, and specific entropy
/// @param eos Pointer to GRHayL EOS parameters struct
/// @param rho Pointer to density
/// @param Y_e Pointer to electron fraction
/// @param S Pointer to specific entropy
void NRPyEOS_enforce_table_bounds_rho_Ye_S(
    const ghl_eos_parameters *restrict eos,
    double *restrict rho,
    double *restrict Y_e,
    double *restrict S) {

  // Enforce bounds on rho
  *rho = MIN(MAX(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = MIN(MAX(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on S
  *S = MIN(MAX(*S, eos->entropy_min), eos->entropy_max);
}

/// @brief Enforce table bounds on density, electron fraction, and pressure
/// @param eos Pointer to GRHayL EOS parameters struct
/// @param rho Pointer to density
/// @param Y_e Pointer to electron fraction
/// @param P Pointer to pressure
void NRPyEOS_enforce_table_bounds_rho_Ye_P(
    const ghl_eos_parameters *restrict eos,
    double *restrict rho,
    double *restrict Y_e,
    double *restrict P) {

  // Enforce bounds on rho
  *rho = MIN(MAX(*rho, eos->rho_min), eos->rho_max);

  // Enforce bounds on Ye
  *Y_e = MIN(MAX(*Y_e, eos->Y_e_min), eos->Y_e_max);

  // Enforce bounds on P
  *P = MIN(MAX(*P, eos->press_min), eos->press_max);
}
