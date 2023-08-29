#include "../../harm_u2p_util.h"

void ghl_func_1D(
      const ghl_eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double dummy,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df,
      int *restrict n_iter) {

  const double Z = x[0];
  double vsq = ghl_vsq_calc(harm_aux, Z);
  const double Zsq = Z*Z;

  // Make sure that v^2 is physically reasonable, and if not make it so:
  const double dv = 1.0e-15;
  vsq = ( vsq < -dv ) ?  0.      : fabs(vsq);
  vsq = ( vsq > 1. )  ?  (1.-dv) : vsq;

  // Compute P from Z and v^2
  const double X  = harm_aux->Bsq + Z;
  const double X3 = X*X*X;
  const double QBsq_Z2 = harm_aux->QdotBsq/Zsq;

  const double dvsq = -2.0*(harm_aux->Qtsq + (QBsq_Z2/Z) * (3*Z*X + harm_aux->Bsq*harm_aux->Bsq)) / X3;

  /*** For hybrid EOS ***/
  double p_tmp, dPdvsq, dPdZ;
  ghl_compute_func_auxiliaries(eos, Z, vsq , harm_aux->D, &p_tmp, &dPdvsq, &dPdZ);
  /*** For hybrid EOS ***/


  // Compute the residual and the needed Jacobian component
  resid[0]  = Z + 0.5 * harm_aux->Bsq * ( 1. + vsq ) - 0.5*QBsq_Z2 + harm_aux->Qdotn - p_tmp;
  jac[0][0] = 1. - dPdZ + QBsq_Z2/Z + 0.5*harm_aux->Bsq*dvsq;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];
  *df   = - resid[0]*resid[0];
  *f    = -0.5 * ( *df );
}
