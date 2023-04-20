#include "induction.h"

double avg(const int size, const double f[size][size][size], const int imin, const int imax, const int jmin, const int jmax, const int kmin, const int kmax) {
  double retval=0.0,num_in_sum=0.0;
  for(int kk=kmin; kk<=kmax; kk++)
    for(int jj=jmin; jj<=jmax; jj++)
      for(int ii=imin; ii<=imax; ii++) {
        retval+=f[kk][jj][ii]; num_in_sum++;
      }
  return retval/num_in_sum;
}
