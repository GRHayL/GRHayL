#ifndef GHL_WENOZ_RECONSTRUCTION_HELPERS_H_
#define GHL_WENOZ_RECONSTRUCTION_HELPERS_H_

#include "ghl_reconstruction.h"

static inline double ghl_wenoz_compute_beta(
      const double *restrict q,
      const int substencil_size,
      const double *restrict beta_matrix) {
  double beta = 0.0;

  for(int i = 0; i < substencil_size; ++i) {
    for(int j = 0; j < substencil_size; ++j) {
      beta += q[i] * beta_matrix[i * substencil_size + j] * q[j];
    }
  }

  return beta;
}

static inline double ghl_wenoz_reconstruct_right_face(
      const int substencil_size,
      const double *restrict U,
      const double *restrict alpha_coeffs,
      const double *restrict gamma,
      const double *restrict beta_matrices,
      const double *restrict tau_coeffs) {
  const double eps = 1e-100;

  double beta[5];
  double w[5];
  double tau = 0.0;
  double q_face = 0.0;
  double wsum = 0.0;

  for(int k = 0; k < substencil_size; ++k) {
    beta[k] =
        ghl_wenoz_compute_beta(
            &U[k],
            substencil_size,
            &beta_matrices[k * substencil_size * substencil_size]) +
        eps;
    tau += tau_coeffs[k] * beta[k];
  }

  tau = fabs(tau);

  for(int k = 0; k < substencil_size; ++k) {
    double q_candidate = 0.0;

    w[k] = gamma[k] * (1.0 + tau / beta[k]) + eps;
    wsum += w[k];

    for(int j = 0; j < substencil_size; ++j) {
      q_candidate += alpha_coeffs[k * substencil_size + j] * U[k + j];
    }

    q_face += w[k] * q_candidate;
  }

  wsum = 1.0 / wsum;
  q_face *= wsum;

  return q_face;
}

#endif // GHL_WENOZ_RECONSTRUCTION_HELPERS_H_
