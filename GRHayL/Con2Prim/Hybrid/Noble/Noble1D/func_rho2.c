#include "../../../utils_Noble.h"

/*********************************************************************************
   func_rho2():

        -- residual/Jacobian routine for the momentum equation after
           substituting Z(rho). The enthalpy uses the full hybrid EOS,
           including its thermal Gamma and piecewise-polytropic integration
           constants.

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
*********************************************************************************/
void ghl_func_rho2(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double dummy,
      const double x[],
      double dx[],
      double *restrict f,
      double *restrict df) {

  const double rho = x[0];
  const int index = ghl_hybrid_find_polytropic_index(eos, rho);
  const double Gamma = eos->Gamma_ppoly[index];
  const double Gamma_th = eos->Gamma_th;
  const double P_cold = eos->K_ppoly[index] * pow(rho, Gamma);
  const double eps_cold = P_cold / (rho * (Gamma - 1.0)) + eos->eps_integ_const[index];
  const double press = harm_aux->W_times_S * pow(rho, Gamma) / harm_aux->D;

  // Enthalpy density and its derivative within the active polytropic piece.
  const double w
        = rho * (1.0 + eps_cold) + (Gamma_th * press - P_cold) / (Gamma_th - 1.0);
  const double dwdrho = 1.0 + eps_cold + P_cold / rho
                        + Gamma * (Gamma_th * press - P_cold) / ((Gamma_th - 1.0) * rho);

  const double rhosq = rho * rho;
  const double Dsq = harm_aux->D * harm_aux->D;
  const double Z = Dsq * w / rhosq;
  const double dZdrho = Dsq * (dwdrho - 2.0 * w / rho) / rhosq;
  const double vsq = (harm_aux->D - rho) * (harm_aux->D + rho) / Dsq;
  const double dvsqdrho = -2.0 * rho / Dsq;
  const double Bsq_plus_Z = harm_aux->Bsq + Z;
  const double momentum_resid = harm_aux->Qtsq - vsq * Bsq_plus_Z * Bsq_plus_Z;

  // Compute the residual and the needed Jacobian component
  const double resid
        = Z * Z * momentum_resid + harm_aux->QdotBsq * (harm_aux->Bsq + 2.0 * Z);
  const double jac
        = 2.0 * dZdrho
                * (Z * momentum_resid - Z * Z * vsq * Bsq_plus_Z + harm_aux->QdotBsq)
          - Z * Z * dvsqdrho * Bsq_plus_Z * Bsq_plus_Z;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = -resid / jac;
  *df = -resid * resid;
  *f = -0.5 * (*df);
}
