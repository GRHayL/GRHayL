#include "reconstruction.h"

/* Function     : ghl_minmod()
 * Description  : evaluates the minmod function, defined as
 *                         { a if |a| < |b| and ab > 0,
 *                result = { b if |b| < |a| and ab > 0
 *                         { 0 if ab <= 0
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_minmod
*/

double ghl_minmod(
      const double a,
      const double b) {

  const double ab = a*b;

  if(ab > 0) {
    if(fabs(a) < fabs(b)) {
      return a;
    } else {
      return b;
    }
  } else {
    return 0.0;
  }
}
