#include "reconstruction.h"

/*
 * Function     : ghl_maxmod()
 * Description  : evaluates the maxmod function, defined as
 *                         { a if |a| > |b| and ab > 0
 *                result = { b if |b| > |a| and ab > 0
 *                         { 0 if ab <= 0
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_maxmod
*/

double ghl_maxmod(
      const double a,
      const double b) {

  const double ab = a*b;

  if(ab > 0) {
    if(fabs(a) > fabs(b)) {
      return a;
    } else {
      return b;
    }
  } else {
    return 0.0;
  }
}
