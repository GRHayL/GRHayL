#include "../harm_u2p_util.h"

/**********************************************************************/
/*********************************************************************************
   func_Z()

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=Z here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
*********************************************************************************/
void ghl_func_Z(
      const ghl_eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double rho_in,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df,
      int *restrict n_iter) {

  const double W = x[0];

  // Set input rho for NR
  double rho      =  rho_in;
  double x_rho[1] = {rho_in};
  // Find rho from W
  int retval, ntries = 0;
  while (  (retval = ghl_newton_raphson_1d(eos, harm_aux, ndim, W, n_iter, &rho, func_rho)) &&  ( ntries++ < 10 )  ) {
    rho     *= 10.;
    x_rho[0] = rho;
  }
  // Set rho to the NR output
  rho = x_rho[0];

  // Auxiliary variable declarations
  // Gamma
  const double Gamma = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos,rho)];
  // Gamma-1
  const double Gm1 = Gamma - 1.0;
  // B^2 + B^2
  const double two_Bsq = harm_aux->Bsq + harm_aux->Bsq;
  // D^2
  const double t2 = harm_aux->D*harm_aux->D;
  // (Q.B)^2 * D^2
  const double t4 = harm_aux->QdotBsq*t2;
  // rho^2
  const double t6 = rho*rho; // Leo says: added, since will be used twice
  // B^4
  const double t7 = harm_aux->Bsq*harm_aux->Bsq;
  // rho^2 - D^2
  const double t15 = t6 - t2; // Leo says: changed from -(D-rho)*(D+rho)
  // D^(-2)
  const double t24 = 1/t2;
  // W + 2B^2
  const double t200 = W + two_Bsq;
  // (Q.B)^2 * B^2 * D^2
  const double t300 = harm_aux->Bsq*t4; // Leo says: changed from QdotBsq*Bsq*t2
  // \tilde{Q}^2 * D^2
  const double t400 = harm_aux->Qtsq*t2;
  // D * Gamma * gamma * S
  const double s200 = harm_aux->D * Gamma * harm_aux->W_times_S;
  // rho^(Gamma - 1)
  const double rho_Gm1 = pow(rho,Gm1);
  // -rho^2 / ( W * rho - rho^(Gamma-1) * D * Gamma * gamma * S )
  const double drho_dW = -t6/( -rho_Gm1*s200 + W*rho); // Leo says: changed rho*rho -> t6
  // rho^(Gamma-1) * drho_dW (see variable above)
  const double t1000 = rho*drho_dW;

  // Compute the residual and the needed Jacobian component
  resid[0]  = (t300+(t4+t4+(t400+t15*(t7+(t200)*W))*W)*W)*t24;
  jac[0][0] = 2*(t4+(t400+t15*t7+(3.0*t15*harm_aux->Bsq+t7*t1000+(t15+t15+t1000*(t200))*W)*W)*W)*t24;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];
  *df   = - resid[0]*resid[0];
  *f    = - 0.5*(*df);
}
