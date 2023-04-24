#include "reconstruction.h"

/* Function    : grhayl_minmod()
 * Description : evaluates the minmod function, defined as
 *
 *
 *          { a if |a| < |b| and ab > 0,
 * result = { b if |b| < |a| and ab > 0
 *          { 0 if ab <= 0
 *
 * Inputs      : a
 *             : b
 *
 * Outputs     : result - minmod evaluation
 */

double grhayl_minmod(const double a,
                     const double b) {

  const double ab = a*b;

  if(      (fabs(a) < fabs(b)) && (ab > 0)) return a;
  else if( (fabs(b) < fabs(a)) && (ab > 0)) return b;
  else                                      return 0.;

  /* code from ChatGPT, but may not be faster  
    // calculate result
    int diff = abs(a) - abs(b);
    int cond1 = (abs(a) < abs(b)) && (ab > 0);
    int cond2 = (abs(b) < abs(a)) && (ab > 0);
    *result = cond1 ? a : (cond2 ? b : 0);
  */
}
