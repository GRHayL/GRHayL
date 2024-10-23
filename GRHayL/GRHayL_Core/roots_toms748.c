/*
 * This is a modified version of the file toms748_solve.hpp from the Boost
 * library (https://www.boost.org). The copyright notice from the original
 * file is maintained below, as well as information of how to obtain a copy
 * of the Boost Software License, Version 1.0. The specific file on which
 * this file is based on can be found at:
 * https://www.boost.org/doc/libs/1_81_0/boost/math/tools/toms748_solve.hpp
 */

//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "ghl_roots.h"

static void
bracket(
    double f(const double, void *restrict),
    void *restrict fparams,
    double *restrict a,
    double *restrict b,
    double c,
    double *restrict fa,
    double *restrict fb,
    double *restrict d,
    double *restrict fd) {

  //
  // Given a point c inside the existing enclosing interval
  // [a, b] sets a = c if f(c) == 0, otherwise finds the new
  // enclosing interval: either [a, c] or [c, b] and sets
  // d and fd to the point that has just been removed from
  // the interval.  In other words d is the third best guess
  // to the root.
  //
  const double tol = 2.0 * DBL_EPSILON;
  //
  // If the interval [a,b] is very small, or if c is too close
  // to one end of the interval then we need to adjust the
  // location of c accordingly:
  //
  if((*b - *a) < 2 * tol * (*a)) {
    c = *a + (*b - *a) / 2;
  }
  else if(c <= *a + fabs(*a) * tol) {
    c = *a + fabs(*a) * tol;
  }
  else if(c >= *b - fabs(*b) * tol) {
    c = *b - fabs(*b) * tol;
  }
  //
  // OK, lets invoke f(c):
  //
  double fc = f(c, fparams);
  //
  // if we have a zero then we have an exact solution to the root:
  //
  if(fc == 0) {
    *a = c;
    *fa = 0;
    *d = 0;
    *fd = 0;
    return;
  }
  //
  // Non-zero fc, update the interval:
  //
  if(sign(*fa) * sign(fc) < 0) {
    *d = *b;
    *fd = *fb;
    *b = c;
    *fb = fc;
  }
  else {
    *d = *a;
    *fd = *fa;
    *a = c;
    *fa = fc;
  }
}

static inline double safe_div(double num, double denom, double r) {

  //
  // return num / denom without overflow,
  // return r if overflow would occur.
  //
  if(fabs(denom) < 1) {
    if(fabs(denom * DBL_MAX) <= fabs(num)) {
      return r;
    }
  }
  return num / denom;
}

static inline double
secant_interpolate(const double a, const double b, const double fa, const double fb) {

  //
  // Performs standard secant interpolation of [a,b] given
  // function evaluations f(a) and f(b).  Performs a bisection
  // if secant interpolation would leave us very close to either
  // a or b.  Rationale: we only call this function when at least
  // one other form of interpolation has already failed, so we know
  // that the function is unlikely to be smooth with a root very
  // close to a or b.
  //
  double tol = 5 * DBL_EPSILON;
  double c = a - (fa / (fb - fa)) * (b - a);
  if((c <= a + fabs(a) * tol) || (c >= b - fabs(b) * tol)) {
    return (a + b) / 2;
  }
  return c;
}

static double quadratic_interpolate(
      const double a,
      const double b,
      const double d,
      const double fa,
      const double fb,
      const double fd,
      unsigned count) {
  //
  // Performs quadratic interpolation to determine the next point,
  // takes count Newton steps to find the location of the
  // quadratic polynomial.
  //
  // Point d must lie outside of the interval [a,b], it is the third
  // best approximation to the root, after a and b.
  //
  // Note: this does not guarantee to find a root
  // inside [a, b], so we fall back to a secant step should
  // the result be out of range.
  //
  // Start by obtaining the coefficients of the quadratic polynomial:
  //
  double B = safe_div(fb - fa, b - a, DBL_MAX);
  double A = safe_div(fd - fb, d - b, DBL_MAX);
  A = safe_div((A - B), (d - a), 0.0);

  if(A == 0) {
    // failure to determine coefficients, try a secant step:
    return secant_interpolate(a, b, fa, fb);
  }
  //
  // Determine the starting point of the Newton steps:
  //
  double c;
  if(sign(A) * sign(fa) > 0) {
    c = a;
  }
  else {
    c = b;
  }
  //
  // Take the Newton steps:
  //
  for(unsigned i = 1; i <= count; ++i) {
    // c -= safe_div(B * c, (B + A * (2 * c - a - b)), 1 + c - a);
    c -= safe_div(fa + (B + A * (c - b)) * (c - a), B + A * (2 * c - a - b), 1 + c - a);
  }
  if((c <= a) || (c >= b)) {
    // Oops, failure, try a secant step:
    c = secant_interpolate(a, b, fa, fb);
  }
  return c;
}

static double cubic_interpolate(
      const double a,
      const double b,
      const double d,
      const double e,
      const double fa,
      const double fb,
      const double fd,
      const double fe) {

  //
  // Uses inverse cubic interpolation of f(x) at points
  // [a,b,d,e] to obtain an approximate root of f(x).
  // Points d and e lie outside the interval [a,b]
  // and are the third and forth best approximations
  // to the root that we have found so far.
  //
  // Note: this does not guarantee to find a root
  // inside [a, b], so we fall back to quadratic
  // interpolation in case of an erroneous result.
  //
  double q11 = (d - e) * fd / (fe - fd);
  double q21 = (b - d) * fb / (fd - fb);
  double q31 = (a - b) * fa / (fb - fa);
  double d21 = (b - d) * fd / (fd - fb);
  double d31 = (a - b) * fb / (fb - fa);
  double q22 = (d21 - q11) * fb / (fe - fb);
  double q32 = (d31 - q21) * fa / (fd - fa);
  double d32 = (d31 - q21) * fd / (fd - fa);
  double q33 = (d32 - q22) * fa / (fe - fa);
  double c = q31 + q32 + q33 + a;

  if((c <= a) || (c >= b)) {
    // Out of bounds step, fall back to quadratic interpolation:
    c = quadratic_interpolate(a, b, d, fa, fb, fd, 3);
  }
  return c;
}

ghl_error_codes_t ghl_toms748(
      double f(const double, void *restrict),
      void *restrict fparams,
      double a,
      double b,
      roots_params *restrict r) {

  // Step 0: Set basic info to the roots_params struct
  sprintf(r->routine_name, __func__);
  r->a = a;
  r->b = b;

  int count = r->max_iters;
  double c, u, fu, d, fd, e, fe;
  static const double mu = 0.5;

  // We proceed by assuming a < b
  if(a >= b) {
    swap(&a, &b);
  }

  // Now compute fa and fb
  double fa = f(a, fparams);
  double fb = f(b, fparams);

  if(sign(fa) * sign(fb) > 0)
    return ghl_error_root_not_bracketed;

  // Check if we already have a root
  if(fabs(b - a) < r->tol || (fa == 0) || (fb == 0)) {
    r->n_iters = 0;
    if(fabs(fa) < fabs(fb)) {
      r->root = a;
      r->residual = fa;
      return ghl_success;
    }
    r->root = b;
    r->residual = fb;
    return ghl_success;
  }

  // dummy value for fd, e and fe:
  fe = e = fd = 1e5;

  if(fa != 0) {
    //
    // On the first step we take a secant step:
    //
    c = secant_interpolate(a, b, fa, fb);
    bracket(f, fparams, &a, &b, c, &fa, &fb, &d, &fd);
    --count;

    if(count && (fa != 0) && fabs(b - a) > r->tol) {
      //
      // On the second step we take a quadratic interpolation:
      //
      c = quadratic_interpolate(a, b, d, fa, fb, fd, 2);
      e = d;
      fe = fd;
      bracket(f, fparams, &a, &b, c, &fa, &fb, &d, &fd);
      --count;
    }
  }

  while(count && (fa != 0) && fabs(b - a) > r->tol) {
    // save our brackets:
    double a0 = a;
    double b0 = b;
    //
    // Starting with the third step taken
    // we can use either quadratic or cubic interpolation.
    // Cubic interpolation requires that all four function values
    // fa, fb, fd, and fe are distinct, should that not be the case
    // then variable prof will get set to true, and we'll end up
    // taking a quadratic step instead.
    //
    double min_diff = 32 * DBL_MIN;
    bool prof = (fabs(fa - fb) < min_diff) || (fabs(fa - fd) < min_diff)
                || (fabs(fa - fe) < min_diff) || (fabs(fb - fd) < min_diff)
                || (fabs(fb - fe) < min_diff) || (fabs(fd - fe) < min_diff);
    if(prof) {
      c = quadratic_interpolate(a, b, d, fa, fb, fd, 2);
    }
    else {
      c = cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
    }
    //
    // re-bracket, and check for termination:
    //
    e = d;
    fe = fd;
    bracket(f, fparams, &a, &b, c, &fa, &fb, &d, &fd);
    if((0 == --count) || (fa == 0) || fabs(b - a) < r->tol) {
      break;
    }

    //
    // Now another interpolated step:
    //
    prof = (fabs(fa - fb) < min_diff) || (fabs(fa - fd) < min_diff)
           || (fabs(fa - fe) < min_diff) || (fabs(fb - fd) < min_diff)
           || (fabs(fb - fe) < min_diff) || (fabs(fd - fe) < min_diff);
    if(prof) {
      c = quadratic_interpolate(a, b, d, fa, fb, fd, 3);
    }
    else {
      c = cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
    }

    //
    // Bracket again, and check termination condition, update e:
    //
    bracket(f, fparams, &a, &b, c, &fa, &fb, &d, &fd);
    if((0 == --count) || (fa == 0) || fabs(b - a) < r->tol) {
      break;
    }

    //
    // Now we take a double-length secant step:
    //
    if(fabs(fa) < fabs(fb)) {
      u = a;
      fu = fa;
    }
    else {
      u = b;
      fu = fb;
    }
    c = u - 2 * (fu / (fb - fa)) * (b - a);
    if(fabs(c - u) > (b - a) / 2) {
      c = a + (b - a) / 2;
    }

    //
    // Bracket again, and check termination condition:
    //
    e = d;
    fe = fd;
    bracket(f, fparams, &a, &b, c, &fa, &fb, &d, &fd);
    if((0 == --count) || (fa == 0) || fabs(b - a) < r->tol) {
      break;
    }

    //
    // And finally... check to see if an additional bisection step is
    // to be taken, we do this if we're not converging fast enough:
    //
    if((b - a) < mu * (b0 - a0)) {
      continue;
    }

    //
    // bracket again on a bisection:
    //
    e = d;
    fe = fd;
    bracket(f, fparams, &a, &b, a + (b - a) / 2, &fa, &fb, &d, &fd);
    --count;
  } // while loop

  r->n_iters = r->max_iters - count;
  if(fabs(fa) < fabs(fb)) {
    r->root = a;
    r->residual = fa;
    return ghl_success;
  }
  r->root = b;
  r->residual = fb;
  return ghl_success;
}
