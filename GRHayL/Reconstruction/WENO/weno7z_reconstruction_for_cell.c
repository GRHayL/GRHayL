#include "wenoz_reconstruction_helpers.h"

/*
 * These tables are for finite-volume WENO-Z reconstruction from cell averages
 * on a uniform grid, with the center cell located at x = 0 and the right face
 * at x = 1/2.
 *
 * For each substencil S_k = {k-3, ..., k}, let p_k(x) be the cubic polynomial
 * whose cell averages match the four cell averages on S_k. Then:
 *   - w7alpha[k][j] is p_k(1/2) written as a linear combination of the cell
 *     averages on substencil k.
 *   - w7gamma[k] are the optimal linear weights such that the weighted sum of
 *     the four candidate reconstructions matches the degree-6 polynomial built
 *     from the full 7-cell stencil.
 *   - w7beta[k] stores the symmetric Jiang-Shu quadratic form B_k in
 *       beta_k = q_k^T B_k q_k,
 *     where q_k is the vector of cell averages on substencil k and
 *       beta_k = sum_{m=1}^{3} int_{-1/2}^{1/2} (d^m p_k / dx^m)^2 dx.
 *   - tau7 = |beta0 + 3 beta1 - 3 beta2 - beta3| is the standard 7th-order
 *     WENO-Z combination from Castro, Costa, and Don.
 *
 * The tables below were generated offline with exact rational arithmetic in
 * SymPy. The full script is derive_wenoz_coefficients.py in this directory.
 * The core calls used for WENO7Z were:
 *
 *   from derive_wenoz_coefficients import candidate_matrix
 *   from derive_wenoz_coefficients import optimal_weights
 *   from derive_wenoz_coefficients import smoothness_matrices
 *
 *   r = 4
 *   alpha = candidate_matrix(r)
 *   gamma = optimal_weights(r, alpha)
 *   beta = smoothness_matrices(r)
 */

void ghl_weno7z_reconstruction_for_cell(
      const double U[7],
      double *restrict Ur,
      double *restrict Ul) {

  const double w7alpha[4][4] = {
      {-1.0 / 4.0, 13.0 / 12.0, -23.0 / 12.0, 25.0 / 12.0},
      {1.0 / 12.0, -5.0 / 12.0, 13.0 / 12.0, 1.0 / 4.0},
      {-1.0 / 12.0, 7.0 / 12.0, 7.0 / 12.0, -1.0 / 12.0},
      {1.0 / 4.0, 13.0 / 12.0, -5.0 / 12.0, 1.0 / 12.0}};

  const double w7gamma[4] = {
      1.0 / 35.0, 12.0 / 35.0, 18.0 / 35.0, 4.0 / 35.0};

  const double w7beta[4][4][4] = {
      {{547.0 / 240.0, -647.0 / 80.0, 2321.0 / 240.0, -309.0 / 80.0},
       {-647.0 / 80.0, 7043.0 / 240.0, -8623.0 / 240.0, 3521.0 / 240.0},
       {2321.0 / 240.0, -8623.0 / 240.0, 11003.0 / 240.0, -1567.0 / 80.0},
       {-309.0 / 80.0, 3521.0 / 240.0, -1567.0 / 80.0, 2107.0 / 240.0}},
      {{89.0 / 80.0, -821.0 / 240.0, 267.0 / 80.0, -247.0 / 240.0},
       {-821.0 / 240.0, 2843.0 / 240.0, -2983.0 / 240.0, 961.0 / 240.0},
       {267.0 / 80.0, -2983.0 / 240.0, 3443.0 / 240.0, -1261.0 / 240.0},
       {-247.0 / 240.0, 961.0 / 240.0, -1261.0 / 240.0, 547.0 / 240.0}},
      {{547.0 / 240.0, -1261.0 / 240.0, 961.0 / 240.0, -247.0 / 240.0},
       {-1261.0 / 240.0, 3443.0 / 240.0, -2983.0 / 240.0, 267.0 / 80.0},
       {961.0 / 240.0, -2983.0 / 240.0, 2843.0 / 240.0, -821.0 / 240.0},
       {-247.0 / 240.0, 267.0 / 80.0, -821.0 / 240.0, 89.0 / 80.0}},
      {{2107.0 / 240.0, -1567.0 / 80.0, 3521.0 / 240.0, -309.0 / 80.0},
       {-1567.0 / 80.0, 11003.0 / 240.0, -8623.0 / 240.0, 2321.0 / 240.0},
       {3521.0 / 240.0, -8623.0 / 240.0, 7043.0 / 240.0, -647.0 / 80.0},
       {-309.0 / 80.0, 2321.0 / 240.0, -647.0 / 80.0, 547.0 / 240.0}}};

  const double tau7[4] = {1.0, 3.0, -3.0, -1.0};

  double Urev[7];

  for(int i = 0; i < 7; ++i) {
    Urev[i] = U[6 - i];
  }

  const double qr =
      ghl_wenoz_reconstruct_right_face(
          4,
          U,
          &w7alpha[0][0],
          w7gamma,
          &w7beta[0][0][0],
          tau7);

  const double ql =
      ghl_wenoz_reconstruct_right_face(
          4,
          Urev,
          &w7alpha[0][0],
          w7gamma,
          &w7beta[0][0][0],
          tau7);

  *Ur = qr;
  *Ul = ql;
}
