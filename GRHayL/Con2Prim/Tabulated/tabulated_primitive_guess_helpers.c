#include "ghl_con2prim.h"

void ghl_tabulated_compute_primitive_guess_auxiliaries(
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      ghl_tabulated_primitive_guess_aux *restrict aux) {

  ghl_compute_SU_Bsq_Ssq_BdotS(metric_adm, cons_undens, prims,
                               aux->SU, &aux->B_squared,
                               &aux->S_squared, &aux->BdotS);

  const double invD = 1.0 / cons_undens->rho;
  aux->q = cons_undens->tau * invD;
  aux->r = aux->S_squared * invD * invD;
  aux->s = aux->B_squared * invD;
  aux->t = aux->BdotS / pow(cons_undens->rho, 1.5);
}

void ghl_tabulated_primitive_guess_from_x_and_W(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_tabulated_primitive_guess_aux *restrict aux,
      double x,
      double W,
      ghl_primitive_quantities *restrict prims) {

  const double q = aux->q;
  const double r = aux->r;
  const double s = aux->s;
  const double t = aux->t;

  if(!(W > 0.0) || !isfinite(W)) {
    double Wminus2 = 1.0 - (x * x * r + (2.0 * x + s) * t * t)
                   / (x * x * (x + s) * (x + s));
    Wminus2 = ghl_clamp(Wminus2, params->inv_sq_max_Lorentz_factor, 1.0);
    W = pow(Wminus2, -0.5);
  }
  else {
    W = ghl_clamp(W, 1.0, params->max_Lorentz_factor);
  }

  prims->rho = cons_undens->rho / W;
  prims->Y_e = cons_undens->Y_e / cons_undens->rho;
  prims->u0  = W * metric_adm->lapseinv;

  prims->eps = - 1.0 + (1.0 - W * W) * x / W
             + W * (1.0 + q - s + t * t / (2.0 * x * x) + s / (2.0 * W * W));

  prims->temperature = eos->T_max;
  ghl_tabulated_enforce_bounds_rho_Ye_eps(eos, &prims->rho, &prims->Y_e, &prims->eps);
  ghl_tabulated_compute_P_S_T_from_eps(eos, prims->rho, prims->Y_e, prims->eps,
                                       &prims->press, &prims->entropy, &prims->temperature);

  const double z = x * prims->rho * W;
  double utildeU[3] = {
    W * (aux->SU[0] + aux->BdotS * prims->BU[0] / z) / (z + aux->B_squared),
    W * (aux->SU[1] + aux->BdotS * prims->BU[1] / z) / (z + aux->B_squared),
    W * (aux->SU[2] + aux->BdotS * prims->BU[2] / z) / (z + aux->B_squared),
  };
  ghl_limit_utilde_and_compute_v(params, metric_adm, utildeU, prims);
}
