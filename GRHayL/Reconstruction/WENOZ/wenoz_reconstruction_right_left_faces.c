#include "reconstruction.h"

/*
 * Function     : ghl_wenoz_reconstruction_right_left_faces()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i+1/2+epsilon)
 *                    Ul(i) = U(i-1/2-epsilon)
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
  const double dc = (dm * dp > 0.0) * 0.5 * (dm + dp);
  return copysign(
      MIN(fabs(dc), alpha * MIN(fabs(dm), fabs(dp))), dc);
}

void ghl_wenoz_reconstruction_right_left_faces(
      const double U[5],
      double *restrict Ur,
      double *restrict Ul) {

  double qr, ql;

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

  beta0 = (beta0 + tau5) / beta0;
  beta1 = (beta1 + tau5) / beta1;
  beta2 = (beta2 + tau5) / beta2;

  double w0 = w5gamma[0] * beta0 + eps;
  double w1 = w5gamma[1] * beta1 + eps;
  double w2 = w5gamma[2] * beta2 + eps;
  double wsum = 1.0 / (w0 + w1 + w2);
  qr =  w0 * (w5alpha[0][0] * q0 + w5alpha[0][1] * q1 + w5alpha[0][2] * q2);
  qr += w1 * (w5alpha[1][0] * q1 + w5alpha[1][1] * q2 + w5alpha[1][2] * q3);
  qr += w2 * (w5alpha[2][0] * q2 + w5alpha[2][1] * q3 + w5alpha[2][2] * q4);
  qr *= wsum;
  const double alpha_l =
      3.0 * wsum * w0 * w1 * w2 /
          (w5gamma[2] * w0 * w1 + w5gamma[1] * w0 * w2 + w5gamma[0] * w1 * w2) +
      eps;

  w0 = w5gamma[0] * beta2 + eps;
  w1 = w5gamma[1] * beta1 + eps;
  w2 = w5gamma[2] * beta0 + eps;
  wsum = 1.0 / (w0 + w1 + w2);
  ql =  w0 * (w5alpha[0][0] * q4 + w5alpha[0][1] * q3 + w5alpha[0][2] * q2);
  ql += w1 * (w5alpha[1][0] * q3 + w5alpha[1][1] * q2 + w5alpha[1][2] * q1);
  ql += w2 * (w5alpha[2][0] * q2 + w5alpha[2][1] * q1 + w5alpha[2][2] * q0);
  ql *= wsum;
  const double alpha_r =
      3.0 * wsum * w0 * w1 * w2 /
          (w5gamma[2] * w0 * w1 + w5gamma[1] * w0 * w2 + w5gamma[0] * w1 * w2) +
      eps;

  double dq = q3 - q2;
  dq = mc(q2 - q1, dq, 2.0);

  const double alpha_lin = 2.0 * alpha_l * alpha_r / (alpha_l + alpha_r);
  *Ur = alpha_lin * qr + (1.0 - alpha_lin) * (q2 + 0.5 * dq);
  *Ul = alpha_lin * ql + (1.0 - alpha_lin) * (q2 - 0.5 * dq);
 }