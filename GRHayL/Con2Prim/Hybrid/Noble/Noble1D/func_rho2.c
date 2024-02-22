#include "../../../utils_Noble.h"

/**********************************************************************/
/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho Qtsq equation with
            the definition of Z
        Z  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho
              substituted in.

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
*********************************************************************************/
void ghl_func_rho2(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double dummy,
      const double x[],
      double dx[],
      double *restrict f,
      double *restrict df) {

  // Set rho
  const double rho = x[0];

  // Auxiliary variable declarations
  // Gamma
  const double Gamma = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos,rho)];
  // Gamma-1
  double Gm1 = Gamma - 1.0;
  // rho^2
  const double rhosq = rho*rho;
  // Gamma/(Gamma-1) * W * S * rho^(Gamma-1)
  const double t100 = (Gamma/Gm1) * harm_aux->W_times_S * pow(rho,Gm1);
  // D^(-2)
  const double t200 = 1.0/(harm_aux->D*harm_aux->D);
  // Z
  const double Z = harm_aux->D * ( harm_aux->D + t100 ) / rho;
  // dZdrho
  const double dZdrho = harm_aux->D * ( -harm_aux->D + t100*(Gamma-2.0) ) / rhosq;
  // Z^2
  const double t1 = Z*Z;
  // B^2 + Z
  const double t2 = harm_aux->Bsq+Z;
  // (B^2 + Z)^2
  const double t3 = t2*t2;
  // v^2
  const double rel_err = (harm_aux->D != 0.0) ? fabs((harm_aux->D-rho)/harm_aux->D) : ( ( rho != 0.0 ) ? fabs((harm_aux->D-rho)/rho) : 0.0 );
  const double vsq = ( rel_err > 1e-15 ) ? ((harm_aux->D-rho)*(harm_aux->D+rho)*t200) : 0.0;
  // d(v^2)/drho
  const double dvsqdrho = -2*rho*t200;
  // B^4
  const double t12 = harm_aux->Bsq*harm_aux->Bsq;
  // v^2 * dwrho
  const double t17 = dZdrho*vsq;

  // Compute the residual and the needed Jacobian component
  const double resid  = t1*(harm_aux->Qtsq-vsq*t3)+harm_aux->QdotBsq*(t2+Z);
  const double jac = 2*harm_aux->QdotBsq*dZdrho
    +((harm_aux->Qtsq-vsq*t12)*2*dZdrho+(-6*t17*harm_aux->Bsq-dvsqdrho*t12
                                        +(-2*dvsqdrho*harm_aux->Bsq-4*t17-dvsqdrho*Z)*Z)*Z)*Z;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid/jac;
  *df   = - resid*resid;
  *f    = - 0.5*(*df);
}
