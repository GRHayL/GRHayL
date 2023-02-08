#include "../harm_u2p_util.h"

double dpdW_calc_vsq(const eos_parameters *restrict eos, const double W, const double vsq);
double dvsq_dW(const harm_aux_vars_struct *restrict harm_aux, double W);

void func_1d_orig(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *f,
      double *df,
      int n,
      double dummy) {

  const double W = x[0];
  double vsq = vsq_calc(harm_aux, W);
  const double Wsq = W*W;

  // Make sure that v^2 is physically reasonable, and if not make it so:
  const double dv = 1.0e-15;
  vsq = ( vsq < -dv ) ?  0.      : fabs(vsq);
  vsq = ( vsq > 1. )  ?  (1.-dv) : vsq;


  // Compute P from W and v^2
  const double dvsq = dvsq_dW(harm_aux, W);

  /*** For hybrid EOS ***/
  const double p_tmp  = pressure_W_vsq(eos, W, vsq , harm_aux->D);
  const double dPdvsq = dpdvsq_calc(eos, W, vsq, harm_aux->D);
  const double dPdW   = dpdW_calc_vsq(eos, W, vsq ) + dPdvsq*dvsq;
  /*** For hybrid EOS ***/


  // Compute the residual and the needed Jacobian component
  resid[0]  = W + 0.5 * harm_aux->Bsq * ( 1. + vsq ) - 0.5*harm_aux->QdotBsq/Wsq + harm_aux->Qdotn - p_tmp;
  jac[0][0] = 1. - dPdW + harm_aux->QdotBsq/(Wsq*W) + 0.5*harm_aux->Bsq*dvsq;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];

  *df = - resid[0]*resid[0];
  *f  = -0.5 * ( *df );

}

double dvsq_dW(const harm_aux_vars_struct *restrict harm_aux, const double W)
{

  double X  = harm_aux->Bsq + W;
  double W3 = W*W*W;
  double X3 = X*X*X;

  return( -2.*( harm_aux->Qtsq/X3 + harm_aux->QdotBsq * (3*W*X + harm_aux->Bsq*harm_aux->Bsq) / ( W3 * X3 ) )  );
}
