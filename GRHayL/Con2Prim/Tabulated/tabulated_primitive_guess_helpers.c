#include "ghl_con2prim.h"

static void set_tabulated_atmosphere_guess(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      ghl_primitive_quantities *restrict prims) {

  const double BU[3] = { prims->BU[0], prims->BU[1], prims->BU[2] };
  ghl_initialize_primitives(
        eos->rho_atm, eos->press_atm, eos->eps_atm,
        -metric_adm->betaU[0], -metric_adm->betaU[1], -metric_adm->betaU[2],
        BU[0], BU[1], BU[2],
        eos->entropy_atm, eos->Y_e_atm, eos->T_atm, prims);
  prims->u0 = metric_adm->lapseinv;
}

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

void ghl_tabulated_primitive_guess_from_x(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_tabulated_primitive_guess_aux *restrict aux,
      double x,
      ghl_primitive_quantities *restrict prims) {

  if(!isfinite(cons_undens->rho) || cons_undens->rho <= 0.0
        || !isfinite(x)
        || !isfinite(aux->q) || !isfinite(aux->r)
        || !isfinite(aux->s) || !isfinite(aux->t)
        || !isfinite(aux->B_squared) || !isfinite(aux->BdotS)
        || !isfinite(aux->SU[0]) || !isfinite(aux->SU[1])
        || !isfinite(aux->SU[2])) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
    return;
  }

  const double q = aux->q;
  const double r = aux->r;
  const double s = aux->s;
  const double t = aux->t;

  const double x_squared = x * x;
  const double x_plus_s = x + s;
  const double numerator = x_squared * r + (2.0 * x + s) * t * t;
  const double denominator = x_squared * x_plus_s * x_plus_s;
  if(!isfinite(numerator) || !isfinite(denominator) || denominator == 0.0) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
    return;
  }
  double Wminus2 = 1.0 - numerator / denominator;
  if(!isfinite(Wminus2)) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
    return;
  }
  Wminus2 = ghl_clamp(Wminus2, params->inv_sq_max_Lorentz_factor, 1.0);
  const double W = pow(Wminus2, -0.5);

  prims->rho = cons_undens->rho / W;
  prims->Y_e = cons_undens->Y_e / cons_undens->rho;
  prims->u0  = W * metric_adm->lapseinv;

  prims->eps = - 1.0 + (1.0 - W * W) * x / W
             + W * (1.0 + q - s + t * t / (2.0 * x * x) + s / (2.0 * W * W));

  if(!isfinite(W) || !isfinite(prims->rho) || prims->rho <= 0.0
        || !isfinite(prims->Y_e) || !isfinite(prims->u0)
        || !isfinite(prims->eps)) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
    return;
  }

  if(ghl_tabulated_enforce_bounds_rho_Ye_eps == NULL
        || ghl_tabulated_compute_P_S_T_from_eps == NULL) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
    return;
  }

  prims->temperature = eos->T_max;
  ghl_tabulated_enforce_bounds_rho_Ye_eps(eos, &prims->rho, &prims->Y_e, &prims->eps);
  const ghl_error_codes_t eos_error = ghl_tabulated_compute_P_S_T_from_eps(
        eos, prims->rho, prims->Y_e, prims->eps,
        &prims->press, &prims->entropy, &prims->temperature);
  if(eos_error != ghl_success) {
    // A failed inversion only rejects this algebraic initial guess; it does not
    // make the conserved state invalid. Preserve rho, Y_e, and the established
    // T_max root-finding seed instead of replacing the composition with
    // atmosphere values. Pressure and entropy are not used by the recovery
    // seed, but keep them finite and deterministic.
    prims->press = eos->press_atm;
    prims->entropy = eos->entropy_atm;
    prims->temperature = eos->T_max;
  }

  const double z = x * prims->rho * W;
  const double velocity_denominator = z + aux->B_squared;
  if(!isfinite(prims->rho) || !isfinite(prims->Y_e)
        || !isfinite(prims->eps) || !isfinite(prims->press)
        || !isfinite(prims->entropy) || !isfinite(prims->temperature)
        || !isfinite(z) || z == 0.0 || !isfinite(velocity_denominator)
        || velocity_denominator == 0.0) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
    return;
  }
  double utildeU[3] = {
    W * (aux->SU[0] + aux->BdotS * prims->BU[0] / z) / velocity_denominator,
    W * (aux->SU[1] + aux->BdotS * prims->BU[1] / z) / velocity_denominator,
    W * (aux->SU[2] + aux->BdotS * prims->BU[2] / z) / velocity_denominator,
  };
  if(!isfinite(utildeU[0]) || !isfinite(utildeU[1]) || !isfinite(utildeU[2])) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
    return;
  }
  ghl_limit_utilde_and_compute_v(params, metric_adm, utildeU, prims);
  if(!isfinite(prims->u0) || !isfinite(prims->vU[0])
        || !isfinite(prims->vU[1]) || !isfinite(prims->vU[2])) {
    set_tabulated_atmosphere_guess(eos, metric_adm, prims);
  }
}
