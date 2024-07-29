#include "reconstruction.h"

/*
 * Function     : ghl_weno5_reconstruction_right_left_faces()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i+1/2+epsilon)
 *                    Ul(i) = U(i-1/2-epsilon)
 *                using the WENO-5 reconstruction algorithm,
 *                i.e. it reconstructs at x-1/2*delta x and 
 *                x+1/2*delta x
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/weno5

The following code has been adapted from the nubhlight code:

Miller, Jonah, Ryan, Benjamin, Dolence, Joshua, and USDOE. lanl/nubhlight. C
omputer software. https://www.osti.gov//servlets/purl/1605099. USDOE. 3 Feb. 2020. Web. doi:10.11578/dc.20200317.1.

copyright 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
operated by Triad National Security, LLC for the U.S.  Department of
Energy/National Nuclear Security Administration. All rights in the
program are reserved by Triad National Security, LLC, and the
U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting
on its behalf a nonexclusive, paid-up, irrevocable worldwide license
in this material to reproduce, prepare derivative works, distribute
copies to the public, perform publicly and display publicly, and to
permit others to do so.

This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
 
2.Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
 
3.Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.
*/

// WENO interpolation. See Tchekhovskoy et al. 2007 (T07), Shu 2011 (S11)
// Implemented by Monika Moscibrodzka
void ghl_weno5_reconstruction_right_left_faces(
      const double U[5],
      double *restrict Ur,
      double *restrict Ul) {

  const double x1 = U[0];
  const double x2 = U[1];
  const double x3 = U[2];
  const double x4 = U[3];
  const double x5 = U[4];

  // S11 1, 2, 3
  double vr[3], vl[3];
  vr[0] = (3. / 8.) * x1 - (5. / 4.) * x2 + (15. / 8.) * x3;
  vr[1] = (-1. / 8.) * x2 + (3. / 4.) * x3 + (3. / 8.) * x4;
  vr[2] = (3. / 8.) * x3 + (3. / 4.) * x4 - (1. / 8.) * x5;

  vl[0] = (3. / 8.) * x5 - (5. / 4.) * x4 + (15. / 8.) * x3;
  vl[1] = (-1. / 8.) * x4 + (3. / 4.) * x3 + (3. / 8.) * x2;
  vl[2] = (3. / 8.) * x3 + (3. / 4.) * x2 - (1. / 8.) * x1;

  // Smoothness indicators, T07 A18 or S11 8
  double beta[3];
  beta[0] = (13. / 12.) * pow(x1 - 2. * x2 + x3, 2) +
            (1. / 4.) * pow(x1 - 4. * x2 + 3. * x3, 2);
  beta[1] =
      (13. / 12.) * pow(x2 - 2. * x3 + x4, 2) + (1. / 4.) * pow(x4 - x2, 2);
  beta[2] = (13. / 12.) * pow(x3 - 2. * x4 + x5, 2) +
            (1. / 4.) * pow(x5 - 4. * x4 + 3. * x3, 2);

  // Nonlinear weights S11 9
  double den, wtr[3], Wr, wr[3], wtl[3], Wl, wl[3], eps;
  eps = 1.e-26;

  den = eps + beta[0];
  den *= den;
  wtr[0] = (1. / 16.) / den;
  den    = eps + beta[1];
  den *= den;
  wtr[1] = (5. / 8.) / den;
  den    = eps + beta[2];
  den *= den;
  wtr[2] = (5. / 16.) / den;
  Wr     = wtr[0] + wtr[1] + wtr[2];
  wr[0]  = wtr[0] / Wr;
  wr[1]  = wtr[1] / Wr;
  wr[2]  = wtr[2] / Wr;

  den = eps + beta[2];
  den *= den;
  wtl[0] = (1. / 16.) / den;
  den    = eps + beta[1];
  den *= den;
  wtl[1] = (5. / 8.) / den;
  den    = eps + beta[0];
  den *= den;
  wtl[2] = (5. / 16.) / den;
  Wl     = wtl[0] + wtl[1] + wtl[2];
  wl[0]  = wtl[0] / Wl;
  wl[1]  = wtl[1] / Wl;
  wl[2]  = wtl[2] / Wl;

  *Ul = vl[0] * wl[0] + vl[1] * wl[1] + vl[2] * wl[2];
  *Ur = vr[0] * wr[0] + vr[1] * wr[1] + vr[2] * wr[2];
 }