#include "nrpyeos_tabulated.h"

#define munu_index(ir, it, iy) ( NRPyEOS_munu_key + NRPyEOS_ntablekeys*((ir) + eos->N_rho*((it) + eos->N_T*(iy))) )

static double find_Ye_st_munu_is_zero(
      const int n,
      const double *restrict Ye,
      const double *restrict munu ) {

  int i0 = -1;
  for(int i=0;i<n-1;i++) {
    if( munu[i]*munu[i+1] < 0 ) {
      i0 = i;
      break;
    }
  }

  /* If munu does not cross zero, return the minimum Ye */
  if( i0 == -1 )
    return Ye[0];

  const int i1  = i0+1;
  const double x0 = Ye[i0];
  const double x1 = Ye[i1];
  const double y0 = munu[i0];
  const double y1 = munu[i1];

  /*
   * The linear interpolation formula reads
   *
   * y(x) = ( y0*(x1-x) + y1*(x-x0) )/(x1-x0) .
   *
   * Let y = munu and x = Ye. We thus want to find
   * x such that y(x) = 0, i.e.,
   *
   *     0 = ( y0*(x1-x) + y1*(x-x0) )/(x1-x0)
   * =>  0 = y0*(x1-x) + y1*(x-x0)
   * =>  0 = y0*x1 - y1*x0 - x*(y0-y1)
   * =>  x = (y0*x1 - y1*x0)/(y0-y1)
   */
  return (y0*x1 - y1*x0)/(y0-y1);
}

void NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T(
      const ghl_eos_parameters *restrict eos,
      const double T,
      double **Ye_of_rho_out ) {

  const int it = ghl_tabulated_get_index_T(eos, T);
  const int nr = eos->N_rho;
  const int ny = eos->N_Ye;

  double *Ye_of_rho  = (double *)malloc(sizeof(double)*nr);
  double *munu_of_Ye = (double *)malloc(sizeof(double)*ny);

  for(int ir=0;ir<nr;ir++) {
    for(int iy=0;iy<ny;iy++)
      munu_of_Ye[iy] = eos->table_all[munu_index(ir, it, iy)];
    Ye_of_rho[ir] = find_Ye_st_munu_is_zero(ny, eos->table_Y_e, munu_of_Ye);
  }

  free(munu_of_Ye);
  *Ye_of_rho_out = Ye_of_rho;
}
