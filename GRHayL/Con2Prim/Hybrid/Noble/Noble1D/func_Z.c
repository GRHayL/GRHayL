#include "../../../utils_Noble.h"

/**
 * @ingroup c2p_internal
 * @brief Evaluates the entropy momentum residual and Newton step in Z.
 * @details First solves for density with the constant-Gamma entropy closure.
 * Every density attempt has a fresh Newton budget and contributes to n_iter.
 * Exhausted retries set harm_aux->error before any unconverged density is used.
 */
void ghl_func_Z(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double rho_in,
      const double x[],
      double dx[],
      double *restrict f,
      double *restrict df) {

  const double Z = x[0];

  // Set input rho for NR
  double rho = rho_in;
  double x_rho[1] = { rho_in };
  // Find rho from Z
  ghl_error_codes_t error;
  for(int ntries = 0; ntries <= 10; ntries++) {
    error = ghl_general_newton_raphson(
          eos, harm_aux, 1, Z, x_rho, ghl_validate_1D_entropy, ghl_func_rho);
    if(error == ghl_success) {
      break;
    }
    rho *= 10.0;
    x_rho[0] = rho;
  }
  if(error != ghl_success) {
    harm_aux->error = error;
    return;
  }
  // The nested solve validates with fabs and rejects rho == 0 through its
  // nonfinite Jacobian, so successful output is finite and strictly positive.
  rho = x_rho[0];

  // Auxiliary variable declarations
  // Gamma
  const double Gamma = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos, rho)];
  // Gamma-1
  const double Gm1 = Gamma - 1.0;
  // B^2 + B^2
  const double two_Bsq = harm_aux->Bsq + harm_aux->Bsq;
  // D^2
  const double t2 = harm_aux->D * harm_aux->D;
  // (Q.B)^2 * D^2
  const double t4 = harm_aux->QdotBsq * t2;
  // rho^2
  const double t6 = rho * rho; // Leo says: added, since will be used twice
  // B^4
  const double t7 = harm_aux->Bsq * harm_aux->Bsq;
  // rho^2 - D^2
  const double t15 = t6 - t2; // Leo says: changed from -(D-rho)*(D+rho)
  // D^(-2)
  const double t24 = 1 / t2;
  // Z + 2B^2
  const double t200 = Z + two_Bsq;
  // (Q.B)^2 * B^2 * D^2
  const double t300 = harm_aux->Bsq * t4; // Leo says: changed from QdotBsq*Bsq*t2
  // \tilde{Q}^2 * D^2
  const double t400 = harm_aux->Qtsq * t2;
  // D * Gamma * gamma * S
  const double s200 = harm_aux->D * Gamma * harm_aux->W_times_S;
  // rho^(Gamma - 1)
  const double rho_Gm1 = pow(rho, Gm1);
  // -rho^2 / ( Z * rho - rho^(Gamma-1) * D * Gamma * gamma * S )
  const double drho_dZ
        = -t6 / (-rho_Gm1 * s200 + Z * rho); // Leo says: changed rho*rho -> t6
  // rho^(Gamma-1) * drho_dZ (see variable above)
  const double t1000 = rho * drho_dZ;

  // Compute the residual and the needed Jacobian component
  const double resid = (t300 + (t4 + t4 + (t400 + t15 * (t7 + (t200)*Z)) * Z) * Z) * t24;
  const double jac = 2
                     * (t4
                        + (t400 + t15 * t7
                           + (3.0 * t15 * harm_aux->Bsq + t7 * t1000
                              + (t15 + t15 + t1000 * (t200)) * Z)
                                   * Z)
                                * Z)
                     * t24;
  if(!isfinite(jac) || jac == 0.0) {
    harm_aux->error = ghl_error_c2p_singular;
    return;
  }
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = -resid / jac;
  *df = -resid * resid;
  *f = -0.5 * (*df);
}
