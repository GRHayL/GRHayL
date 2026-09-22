#include "ghl_reconstruction.h"

/*
 * Function     : ghl_wenoz_reconstruction_right_left_faces()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i+1/2-epsilon)
 *                    Ul(i) = U(i-1/2+epsilon)
 *                using the WENO-z reconstruction algorithm,
 *                i.e. it reconstructs at x-1/2*delta x and
 *                x+1/2*delta x
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/wenoz

The following code has been adapted from the phoebus code:

//========================================================================================
// (C) (or copyright) 2021. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
// for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================

BSD 3-Clause License

Copyright (c) 2021, Los Alamos National Laboratory
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// WENO interpolation. See Tchekhovskoy et al. 2007 (T07), Shu 2011 (S11)
// Implemented by Monika Moscibrodzka

static double mc(const double dm, const double dp, const double alpha) {
  const double dc = ((dm > 0.0 && dp > 0.0) || (dm < 0.0 && dp < 0.0)) * 0.5 * (dm + dp);
  return copysign(
      fmin(fabs(dc), alpha * fmin(fabs(dm), fabs(dp))), dc);
}

static void wenoz_weights(
      const double beta[3],
      const double tau5,
      const double gamma[3],
      double weights[3],
      double *restrict alpha) {

  const double eps = 1e-100;
  const double beta_min = fmin(beta[0], fmin(beta[1], beta[2]));
  const double scale = fmax(beta_min, tau5);
  const double scaled_beta_min = beta_min / scale;
  const double scaled_tau5 = tau5 / scale;
  const double scaled_sum = scaled_beta_min + scaled_tau5;
  const double common_beta = scaled_beta_min / scaled_sum;
  const double common_tau5 = scaled_tau5 / scaled_sum;

  for(int i = 0; i < 3; i++) {
    weights[i] = gamma[i] * (common_beta + common_tau5 * (beta_min / beta[i]))
                 + eps * common_beta;
  }

  const double weight_sum = weights[0] + weights[1] + weights[2];
  for(int i = 0; i < 3; i++) {
    weights[i] /= weight_sum;
  }

  const double weight_min = fmin(weights[0], fmin(weights[1], weights[2]));
  if(weight_min == 0.0) {
    *alpha = eps;
  }
  else {
    *alpha = 3.0 * weight_min
                   / (gamma[2] * (weight_min / weights[2])
                      + gamma[1] * (weight_min / weights[1])
                      + gamma[0] * (weight_min / weights[0]))
             + eps;
  }
}

static __attribute__((noinline, cold)) void wenoz_reconstruct_fallback(
      const double U[5],
      const double beta0,
      const double beta1,
      const double beta2,
      const double tau5,
      double *restrict Ur,
      double *restrict Ul) {

  const double w5alpha[3][3] = { { 1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0 },
                                 { -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0 },
                                 { 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0 } };
  const double w5gamma[3] = { 0.1, 0.6, 0.3 };

  double weights_l[3], weights_r[3], alpha_l, alpha_r;
  const double beta[3] = { beta0, beta1, beta2 };
  wenoz_weights(beta, tau5, w5gamma, weights_l, &alpha_l);
  const double beta_reversed[3] = { beta2, beta1, beta0 };
  wenoz_weights(beta_reversed, tau5, w5gamma, weights_r, &alpha_r);

  double qr = weights_l[0]
              * (w5alpha[0][0] * U[0] + w5alpha[0][1] * U[1] + w5alpha[0][2] * U[2]);
  qr += weights_l[1]
        * (w5alpha[1][0] * U[1] + w5alpha[1][1] * U[2] + w5alpha[1][2] * U[3]);
  qr += weights_l[2]
        * (w5alpha[2][0] * U[2] + w5alpha[2][1] * U[3] + w5alpha[2][2] * U[4]);

  double ql = weights_r[0]
              * (w5alpha[0][0] * U[4] + w5alpha[0][1] * U[3] + w5alpha[0][2] * U[2]);
  ql += weights_r[1]
        * (w5alpha[1][0] * U[3] + w5alpha[1][1] * U[2] + w5alpha[1][2] * U[1]);
  ql += weights_r[2]
        * (w5alpha[2][0] * U[2] + w5alpha[2][1] * U[1] + w5alpha[2][2] * U[0]);

  double dq = U[3] - U[2];
  dq = mc(U[2] - U[1], dq, 2.0);

  const double alpha_lin = 2.0 * alpha_l * alpha_r / (alpha_l + alpha_r);
  *Ur = alpha_lin * qr + (1.0 - alpha_lin) * (U[2] + 0.5 * dq);
  *Ul = alpha_lin * ql + (1.0 - alpha_lin) * (U[2] - 0.5 * dq);
}

void ghl_wenoz_reconstruction_right_left_faces(
      const double U[5],
      double *restrict Ur,
      double *restrict Ul) {

  const double q0 = U[0];
  const double q1 = U[1];
  const double q2 = U[2];
  const double q3 = U[3];
  const double q4 = U[4];

  const double w5alpha[3][3] = {{ 1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0},
                                {-1.0 / 6.0,  5.0 / 6.0,  1.0 / 3.0},
                                { 1.0 / 3.0,  5.0 / 6.0, -1.0 / 6.0}};

  const double w5gamma[3] = {0.1, 0.6, 0.3};
  const double eps = 1e-100;
  const double thirteen_thirds = 13.0 / 3.0;

  double a = q0 - 2 * q1 + q2;
  double b = q0 - 4.0 * q1 + 3.0 * q2;
  double beta0 = thirteen_thirds * a * a + b * b + eps;
  a = q1 - 2.0 * q2 + q3;
  b = q3 - q1;
  double beta1 = thirteen_thirds * a * a + b * b + eps;
  a = q2 - 2.0 * q3 + q4;
  b = q4 - 4.0 * q3 + 3.0 * q2;
  double beta2 = thirteen_thirds * a * a + b * b + eps;
  const double tau5 = fabs(beta2 - beta0);

  // With beta >= 1e-100, this bounds the two possibly large raw weights so
  // their product remains finite.
  const bool use_fast_weights = tau5 <= 1e50;
  if(use_fast_weights) {
    const double beta_scaled[3]
          = { 1.0 + tau5 / beta0, 1.0 + tau5 / beta1, 1.0 + tau5 / beta2 };

    double w0 = w5gamma[0] * beta_scaled[0] + eps;
    double w1 = w5gamma[1] * beta_scaled[1] + eps;
    double w2 = w5gamma[2] * beta_scaled[2] + eps;
    double wsum = 1.0 / (w0 + w1 + w2);
    double qr = w0 * (w5alpha[0][0] * q0 + w5alpha[0][1] * q1 + w5alpha[0][2] * q2);
    qr += w1 * (w5alpha[1][0] * q1 + w5alpha[1][1] * q2 + w5alpha[1][2] * q3);
    qr += w2 * (w5alpha[2][0] * q2 + w5alpha[2][1] * q3 + w5alpha[2][2] * q4);
    qr *= wsum;
    const double alpha_l
          = 3.0 * wsum * w0 * w1 * w2
                  / (w5gamma[2] * w0 * w1 + w5gamma[1] * w0 * w2 + w5gamma[0] * w1 * w2)
            + eps;

    w0 = w5gamma[0] * beta_scaled[2] + eps;
    w1 = w5gamma[1] * beta_scaled[1] + eps;
    w2 = w5gamma[2] * beta_scaled[0] + eps;
    wsum = 1.0 / (w0 + w1 + w2);
    double ql = w0 * (w5alpha[0][0] * q4 + w5alpha[0][1] * q3 + w5alpha[0][2] * q2);
    ql += w1 * (w5alpha[1][0] * q3 + w5alpha[1][1] * q2 + w5alpha[1][2] * q1);
    ql += w2 * (w5alpha[2][0] * q2 + w5alpha[2][1] * q1 + w5alpha[2][2] * q0);
    ql *= wsum;
    const double alpha_r
          = 3.0 * wsum * w0 * w1 * w2
                  / (w5gamma[2] * w0 * w1 + w5gamma[1] * w0 * w2 + w5gamma[0] * w1 * w2)
            + eps;
    double dq = q3 - q2;
    dq = mc(q2 - q1, dq, 2.0);

    const double alpha_lin = 2.0 * alpha_l * alpha_r / (alpha_l + alpha_r);
    *Ur = alpha_lin * qr + (1.0 - alpha_lin) * (q2 + 0.5 * dq);
    *Ul = alpha_lin * ql + (1.0 - alpha_lin) * (q2 - 0.5 * dq);
    return;
  }

  wenoz_reconstruct_fallback(U, beta0, beta1, beta2, tau5, Ur, Ul);
}
