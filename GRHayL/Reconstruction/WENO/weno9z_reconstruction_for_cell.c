#include "wenoz_reconstruction_helpers.h"

/*
 * These tables are for finite-volume WENO-Z reconstruction from cell averages
 * on a uniform grid, with the center cell located at x = 0 and the right face
 * at x = 1/2.
 *
 * For each substencil S_k = {k-4, ..., k}, let p_k(x) be the quartic
 * polynomial whose cell averages match the five cell averages on S_k. Then:
 *   - w9alpha[k][j] is p_k(1/2) written as a linear combination of the cell
 *     averages on substencil k.
 *   - w9gamma[k] are the optimal linear weights such that the weighted sum of
 *     the five candidate reconstructions matches the degree-8 polynomial built
 *     from the full 9-cell stencil.
 *   - w9beta[k] stores the symmetric Jiang-Shu quadratic form B_k in
 *       beta_k = q_k^T B_k q_k,
 *     where q_k is the vector of cell averages on substencil k and
 *       beta_k = sum_{m=1}^{4} int_{-1/2}^{1/2} (d^m p_k / dx^m)^2 dx.
 *   - tau9 = |beta0 + 2 beta1 - 6 beta2 + 2 beta3 + beta4| is the standard
 *     9th-order WENO-Z combination from Castro, Costa, and Don.
 *
 * The tables below were generated offline with exact rational arithmetic in
 * SymPy. The full script is derive_wenoz_coefficients.py in this directory.
 * The core calls used for WENO9Z were:
 *
 *   from derive_wenoz_coefficients import candidate_matrix
 *   from derive_wenoz_coefficients import optimal_weights
 *   from derive_wenoz_coefficients import smoothness_matrices
 *
 *   r = 5
 *   alpha = candidate_matrix(r)
 *   gamma = optimal_weights(r, alpha)
 *   beta = smoothness_matrices(r)
 */

void ghl_weno9z_reconstruction_for_cell(
      const double U[9],
      double *restrict Ur,
      double *restrict Ul) {

  const double w9alpha[5][5] = {
      {1.0 / 5.0, -21.0 / 20.0, 137.0 / 60.0, -163.0 / 60.0, 137.0 / 60.0},
      {-1.0 / 20.0, 17.0 / 60.0, -43.0 / 60.0, 77.0 / 60.0, 1.0 / 5.0},
      {1.0 / 30.0, -13.0 / 60.0, 47.0 / 60.0, 9.0 / 20.0, -1.0 / 20.0},
      {-1.0 / 20.0, 9.0 / 20.0, 47.0 / 60.0, -13.0 / 60.0, 1.0 / 30.0},
      {1.0 / 5.0, 77.0 / 60.0, -43.0 / 60.0, 17.0 / 60.0, -1.0 / 20.0}};

  const double w9gamma[5] = {
      1.0 / 126.0, 10.0 / 63.0, 10.0 / 21.0, 20.0 / 63.0, 5.0 / 126.0};

  const double w9beta[5][5][5] = {
      {{11329.0 / 2520.0, -208501.0 / 10080.0, 121621.0 / 3360.0,
        -288007.0 / 10080.0, 86329.0 / 10080.0},
       {-208501.0 / 10080.0, 482963.0 / 5040.0, -142033.0 / 840.0,
        679229.0 / 5040.0, -411487.0 / 10080.0},
       {121621.0 / 3360.0, -142033.0 / 840.0, 507131.0 / 1680.0,
        -68391.0 / 280.0, 252941.0 / 3360.0},
       {-288007.0 / 10080.0, 679229.0 / 5040.0, -68391.0 / 280.0,
        1020563.0 / 5040.0, -649501.0 / 10080.0},
       {86329.0 / 10080.0, -411487.0 / 10080.0, 252941.0 / 3360.0,
        -649501.0 / 10080.0, 53959.0 / 2520.0}},
      {{1727.0 / 1260.0, -60871.0 / 10080.0, 33071.0 / 3360.0,
        -70237.0 / 10080.0, 18079.0 / 10080.0},
       {-60871.0 / 10080.0, 138563.0 / 5040.0, -3229.0 / 70.0,
        168509.0 / 5040.0, -88297.0 / 10080.0},
       {33071.0 / 3360.0, -3229.0 / 70.0, 135431.0 / 1680.0,
        -25499.0 / 420.0, 55051.0 / 3360.0},
       {-70237.0 / 10080.0, 168509.0 / 5040.0, -25499.0 / 420.0,
        242723.0 / 5040.0, -140251.0 / 10080.0},
       {18079.0 / 10080.0, -88297.0 / 10080.0, 55051.0 / 3360.0,
        -140251.0 / 10080.0, 11329.0 / 2520.0}},
      {{1727.0 / 1260.0, -51001.0 / 10080.0, 7547.0 / 1120.0,
        -38947.0 / 10080.0, 8209.0 / 10080.0},
       {-51001.0 / 10080.0, 104963.0 / 5040.0, -24923.0 / 840.0,
        89549.0 / 5040.0, -38947.0 / 10080.0},
       {7547.0 / 1120.0, -24923.0 / 840.0, 77051.0 / 1680.0,
        -24923.0 / 840.0, 7547.0 / 1120.0},
       {-38947.0 / 10080.0, 89549.0 / 5040.0, -24923.0 / 840.0,
        104963.0 / 5040.0, -51001.0 / 10080.0},
       {8209.0 / 10080.0, -38947.0 / 10080.0, 7547.0 / 1120.0,
        -51001.0 / 10080.0, 1727.0 / 1260.0}},
      {{11329.0 / 2520.0, -140251.0 / 10080.0, 55051.0 / 3360.0,
        -88297.0 / 10080.0, 18079.0 / 10080.0},
       {-140251.0 / 10080.0, 242723.0 / 5040.0, -25499.0 / 420.0,
        168509.0 / 5040.0, -70237.0 / 10080.0},
       {55051.0 / 3360.0, -25499.0 / 420.0, 135431.0 / 1680.0,
        -3229.0 / 70.0, 33071.0 / 3360.0},
       {-88297.0 / 10080.0, 168509.0 / 5040.0, -3229.0 / 70.0,
        138563.0 / 5040.0, -60871.0 / 10080.0},
       {18079.0 / 10080.0, -70237.0 / 10080.0, 33071.0 / 3360.0,
        -60871.0 / 10080.0, 1727.0 / 1260.0}},
      {{53959.0 / 2520.0, -649501.0 / 10080.0, 252941.0 / 3360.0,
        -411487.0 / 10080.0, 86329.0 / 10080.0},
       {-649501.0 / 10080.0, 1020563.0 / 5040.0, -68391.0 / 280.0,
        679229.0 / 5040.0, -288007.0 / 10080.0},
       {252941.0 / 3360.0, -68391.0 / 280.0, 507131.0 / 1680.0,
        -142033.0 / 840.0, 121621.0 / 3360.0},
       {-411487.0 / 10080.0, 679229.0 / 5040.0, -142033.0 / 840.0,
        482963.0 / 5040.0, -208501.0 / 10080.0},
       {86329.0 / 10080.0, -288007.0 / 10080.0, 121621.0 / 3360.0,
        -208501.0 / 10080.0, 11329.0 / 2520.0}}};

  const double tau9[5] = {1.0, 2.0, -6.0, 2.0, 1.0};

  double Urev[9];

  for(int i = 0; i < 9; ++i) {
    Urev[i] = U[8 - i];
  }

  const double qr =
      ghl_wenoz_reconstruct_right_face(
          5,
          U,
          &w9alpha[0][0],
          w9gamma,
          &w9beta[0][0][0],
          tau9);

  const double ql =
      ghl_wenoz_reconstruct_right_face(
          5,
          Urev,
          &w9alpha[0][0],
          w9gamma,
          &w9beta[0][0][0],
          tau9);

  *Ur = qr;
  *Ul = ql;
}
