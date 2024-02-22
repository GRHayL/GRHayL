#include "../../../utils_Noble.h"

/**********************************************************************/
/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho from W via the polytrope:

        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
      IGM update:
        eos   = struct containing polytropic eos quantities
    W_times_S = HARM's Sc variable, set in the Utoprim_1d_ee() function
        W_in  = HARM's W_for_gnr2 variable, set in the Utoprim_new_body() function
          D   = W * rho, set in the Utoprim_new_body() functoin
*********************************************************************************/
// for the isentropic version: eq. (27)
void ghl_func_rho(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double W_in,
      const double x[],
      double dx[],
      double *restrict f,
      double *restrict df) {

  // Set rho and W
  const double rho = x[0];
  const double W   = W_in;

  // Auxiliary variable declarations
  // Gamma
  const double Gamma = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos,rho)];
  // Gamma-1
  const double Gm1 = Gamma - 1.0;
  // rho^(Gamma-1)
  const double t40 = pow(rho,Gm1);
  // rho^(Gamma-2)
  const double t14 = t40/rho;
  // Gamma/(Gamma-1) * W * S
  const double s100 = Gamma/Gm1 * harm_aux->W_times_S;
  // D * Gamma * W * S
  const double s200 = harm_aux->D * Gamma * harm_aux->W_times_S;

  // Compute the residual and the needed Jacobian component
  const double resid  = (rho*W+(-t40*s100-harm_aux->D)*harm_aux->D);
  const double jac = -t14*s200 + W;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid/jac;
  *df   = - resid*resid;
  *f    = - 0.5*(*df);
}
