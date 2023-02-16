#include "con2prim.h"

#define ITMAX 300   //Maximum allowed number of iterations.
#define PREC 3.0e-15
#define EXTRA_ITER 0
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double
zbrent(
      double func(const eos_parameters *restrict,
                  const conservative_quantities *restrict,
                  const double *restrict,
                  const bool,
                  const double,
                  primitive_quantities *restrict,
                  double *restrict),
      const eos_parameters *restrict eos,
      const conservative_quantities *restrict cons_undens,
      const double *restrict params,
      const bool evolve_T,
      primitive_quantities *restrict prims,
      double *restrict temp_guess,
      double x1,
      double x2,
      double tol_x,
      con2prim_diagnostics *restrict stats ) {

  // Implementation of Brent’s method to find the root of a function func to accuracy tol_x. The root
  // must be bracketed by x1 and x2.
  // Based on W. H. Press et al. 2007, Numerical Recipes in C: The Art of Scientific Computing, 3rd edition,
  // Cambridge University Press, ISBN 978-0-521-88068-8

  int iter = 0;

  //termination critera definitions
  double ans = 0.0;

  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1=0.0,min2=0.0;
  double temp_guess_a = *temp_guess;
  double temp_guess_b = *temp_guess;
  double temp_guess_c = *temp_guess;
  double fa=func(eos, cons_undens, params, evolve_T, a, prims, &temp_guess_a);
  double fb=func(eos, cons_undens, params, evolve_T, b, prims, &temp_guess_b);
  double fc,p,q,r,s,tol1;

  // root must be bracketed
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    stats->c2p_failed = true;
    return b;
  }
  fc=fb;
  temp_guess_c = temp_guess_b;

  int i_extra = 0;
  int doing_extra = 0;

  if (EXTRA_ITER==0) {
    i_extra=-1;
  }

  bool keep_iterating = 1;
  double maxerror = 0.0;

  while (keep_iterating) {

    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      // reset c, adjust bracket interval d
      c=a;
      fc=fa;
      e=d=b-a;
      temp_guess_c = temp_guess_a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
      temp_guess_a = temp_guess_b;
      temp_guess_b = temp_guess_c;
      temp_guess_c = temp_guess_a;
    }

    tol1=2.0*PREC*fabs(b)+0.5*tol_x;
    maxerror=0.5*(c-b);
    // check convergence
    if (fabs(maxerror) <= tol1 || fb == 0.0){
      // we are done
      stats->c2p_failed = false;
      keep_iterating = 0;
      ans = b;
    }
    else if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      // inverse quadratic interpolation
      s=fb/fa;

      if (a == c) {
        // fake quad interpolation
        p=2.0*maxerror*s;
        q=1.0-s;
      } else {
        // actual quad interpolation
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*maxerror*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }

      if (p > 0.0) q = -q;  // check whether in bounds
      p=fabs(p);
      min1=3.0*maxerror*q-fabs(tol1*q);
      min2=fabs(e*q);

      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        // interpolation successful
        e=d;
        d=p/q;
      } else {
        // interpolation failed, use bisection
        d=maxerror;
        e=d;
      }

      // update a to be last best guess
      a=b;
      fa=fb;
      temp_guess_a = temp_guess_b;

      // update trial root
      if (fabs(d) > tol1)  b += d;
      else                 b += SIGN(tol1,maxerror);
      func(eos, cons_undens, params, evolve_T, b, prims, &temp_guess_b);


    } else {

      // Bounds decreasing too slowly, use bisection

      d=maxerror;
      e=d;

      // update a to be last best guess
      a=b;
      fa=fb;
      temp_guess_a = temp_guess_b;

      // update trial root
      if (fabs(d) > tol1)  b += d;
      else                 b += SIGN(tol1, maxerror);
      func(eos, cons_undens, params, evolve_T, b, prims, &temp_guess_b);
    }

    ++iter;

    // termination criterion
    if( (fabs(maxerror) <= tol1 ) && (doing_extra == 0) && (EXTRA_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(maxerror) <= tol1)&&(doing_extra == 0))
        || (i_extra > EXTRA_ITER) || (iter >= (ITMAX)) ) {
      keep_iterating = 0;
    }

  } //end for loop

  if( (!isfinite(fb)) ) {
    stats->c2p_failed = true;
    return b;
  } if( (fabs(maxerror) <= tol1 || fb == 0.0)){
    stats->c2p_failed = false;
    return ans;
  } else if( (fabs(maxerror) <= tol1) && (fabs(maxerror) > tol1) ){
    stats->c2p_failed = false;
    return ans;
  } else {
    stats->c2p_failed = true;
    return b;
  }

}
