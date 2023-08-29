#include "../../../utils_Noble.h"

void ghl_func_2D(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double dummy,
      const double x[],
      double dx[],
      double resid[],
      double jac[][2],
      double *restrict f,
      double *restrict df) {

  const double Z   = x[0];
  const double vsq = x[1];
  const double Zsq = Z*Z;

  /*** For hybrid EOS ***/
  double p_tmp, dPdvsq, dPdZ;
  ghl_compute_func_auxiliaries(eos, Z, vsq , harm_aux->D, &p_tmp, &dPdvsq, &dPdZ);
  /*** For hybrid EOS ***/

  // These expressions were calculated using Mathematica, but made into efficient
  // code using Maple.  Since we know the analytic form of the equations, we can
  // explicitly calculate the Newton-Raphson step:

  const double t2  = -0.5*harm_aux->Bsq+dPdvsq;
  const double t3  = harm_aux->Bsq+Z;
  const double t4  = t3*t3;
  const double t9  = 1/Zsq;
  const double t11 = harm_aux->Qtsq-vsq*t4+harm_aux->QdotBsq*(harm_aux->Bsq+2.0*Z)*t9;
  const double t16 = harm_aux->QdotBsq*t9;
  const double t18 = -harm_aux->Qdotn-0.5*harm_aux->Bsq*(1.0+vsq)+0.5*t16-Z+p_tmp;
  const double t21 = 1/t3;
  const double t23 = 1/Z;
  const double t24 = t16*t23;
  const double t25 = -1.0+dPdZ-t24;
  const double t35 = t25*t3+(harm_aux->Bsq-2.0*dPdvsq)*(harm_aux->QdotBsq+vsq*Zsq*Z)*t9*t23;
  const double t36 = 1/t35;
  const double t40 = (vsq+t24)*t3;
  dx[0] = -(t2*t11+t4*t18)*t21*t36;
  dx[1] = -(-t25*t11-2.0*t40*t18)*t21*t36;
  //detJ = t3*t35; // <- set but not used...
  jac[0][0] = -2.0*t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];
  *f  = -0.5 * ( *df );

}
