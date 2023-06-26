#include "reconstruction.h"

/* Function    : ghl_maxmod()
 * Description : evaluates the minmod function, defined as
 *
 *
 *          { a if |a| > |b| and ab > 0,
 * result = { b if |b| > |a| and ab > 0
 *          { 0 if ab <= 0
 *
 * Inputs      : a
 *             : b
 *
 * Outputs     : result - minmod evaluation
 */

double ghl_maxmod(
      const double a,
      const double b) {

  const double ab = a*b;

  if(      (fabs(a) > fabs(b)) && (ab > 0)) return a;
  else if( (fabs(b) > fabs(a)) && (ab > 0)) return b;
  else                                      return 0.;
}
