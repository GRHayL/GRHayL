#include "../harm_u2p_util.h"

void grhayl_func_1d_orig(
      const eos_parameters *restrict eos,
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

  const double W = x[0];
  double vsq = grhayl_vsq_calc(harm_aux, W);
  const double Wsq = W*W;

  //TODO: consider adding fabs to vsq calc and replacing below with
  //  vsq = ( ( vtmp > 1. ) ? (1.0 - 1.e-15) : vtmp );
  // Make sure that v^2 is physically reasonable, and if not make it so:
  const double dv = 1.0e-15;
  vsq = ( vsq < -dv ) ?  0.      : fabs(vsq);
  vsq = ( vsq > 1. )  ?  (1.-dv) : vsq;

  // Compute P from W and v^2
  const double dvsq = grhayl_dvsq_dW(harm_aux, W);

  /*** For hybrid EOS ***/
  const double p_tmp  = grhayl_pressure_W_vsq(eos, W, vsq , harm_aux->D);
  const double dPdvsq = grhayl_dpdvsq_calc(eos, W, vsq, harm_aux->D);
  const double dPdW   = grhayl_dpdW_calc_vsq(eos, W, vsq ) + dPdvsq*dvsq;
  /*** For hybrid EOS ***/


  // Compute the residual and the needed Jacobian component
  resid[0]  = W + 0.5 * harm_aux->Bsq * ( 1. + vsq ) - 0.5*harm_aux->QdotBsq/Wsq + harm_aux->Qdotn - p_tmp;
  jac[0][0] = 1. - dPdW + harm_aux->QdotBsq/(Wsq*W) + 0.5*harm_aux->Bsq*dvsq;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];
  *df   = - resid[0]*resid[0];
  *f    = -0.5 * ( *df );
}
