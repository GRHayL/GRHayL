#include "../../utils_Noble.h"
#include "ghl_nrpyeos_tabulated.h"

#include <float.h>

/**
 * @ingroup c2p_internal
 * @brief Finds enthalpy bounds over temperature at fixed table coordinates.
 */
static ghl_error_codes_t get_min_max(
      const int nx,
      const int ny,
      const int nz,
      const double *restrict x_arr,
      const double *restrict y_arr,
      const double *restrict f_arr,
      const double x,
      const double y,
      double *restrict f_min,
      double *restrict f_max) {

  if(!isfinite(x) || !isfinite(y) || x < x_arr[0] || x > x_arr[nx - 1] || y < y_arr[0]
     || y > y_arr[ny - 1]) {
    return ghl_error_table_bisection;
  }

  const double dx = x_arr[1] - x_arr[0];
  const double dy = y_arr[1] - y_arr[0];

  const int ix = ghl_iclamp((x - x_arr[0]) / dx, 0, nx - 2);
  const int iy = ghl_iclamp((y - y_arr[0]) / dy, 0, ny - 2);

  const double tx = (x - x_arr[ix]) / dx;
  const double ty = (y - y_arr[iy]) / dy;

  const double w00 = (1.0 - tx) * (1.0 - ty);
  const double w10 = tx * (1.0 - ty);
  const double w01 = (1.0 - tx) * ty;
  const double w11 = tx * ty;

  double min = +DBL_MAX;
  double max = -DBL_MAX;

  for(int iz = 0; iz < nz; iz++) {
#define GHL_IDX3D(ix, iy, iz) ((ix) + nx * ((iz) + nz * (iy)))
    const int i00 = GHL_IDX3D(ix + 0, iy + 0, iz);
    const int i01 = GHL_IDX3D(ix + 0, iy + 1, iz);
    const int i10 = GHL_IDX3D(ix + 1, iy + 0, iz);
    const int i11 = GHL_IDX3D(ix + 1, iy + 1, iz);
#undef GHL_IDX3D

    const double f
          = w00 * f_arr[i00] + w10 * f_arr[i10] + w01 * f_arr[i01] + w11 * f_arr[i11];

    if(!isfinite(f)) {
      return ghl_error_table_bisection;
    }
    if(f < min) {
      min = f;
    }
    if(f > max) {
      max = f;
    }
  }

  *f_min = min;
  *f_max = max;

  return ghl_success;
}

/**
 * @ingroup c2p_internal
 * @brief Bounds a finite state for tabulated enthalpy inversion.
 * @return ghl_success or a recoverable table-inversion error.
 */
static ghl_error_codes_t enforce_bounds_rho_Ye_h(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict Y_e,
      double *restrict h) {

  if(!isfinite(*rho) || !isfinite(*Y_e) || !isfinite(*h)) {
    return ghl_error_table_bisection;
  }

  *rho = ghl_clamp(*rho, eos->rho_min, eos->rho_max);
  *Y_e = ghl_clamp(*Y_e, eos->Y_e_min, eos->Y_e_max);

  // Runtime density bounds lie inside the table. Correct log(exp(endpoint))
  // roundoff in table coordinates before the strict interpolation range check.
  const double logrho = log(*rho);
  if(!isfinite(logrho)) {
    return ghl_error_table_bisection;
  }
  const double bounded_logrho
        = ghl_clamp(logrho, eos->table_logrho[0], eos->table_logrho[eos->N_rho - 1]);
  double logh_min = 0.0;
  double logh_max = 0.0;
  const ghl_error_codes_t error = get_min_max(
        eos->N_rho, eos->N_Ye, eos->N_T, eos->table_logrho, eos->table_Y_e,
        eos->table_logh, bounded_logrho, *Y_e, &logh_min, &logh_max);
  if(error != ghl_success) {
    return error;
  }

  const double h_min = exp(logh_min);
  const double h_max = exp(logh_max);
  if(!isfinite(h_min) || !isfinite(h_max) || h_min <= 0.0) {
    return ghl_error_table_bisection;
  }
  *h = ghl_clamp(*h, h_min, h_max);
  return ghl_success;
}

/**
 * @ingroup c2p_internal
 * @brief Initializes the tabulated Noble enthalpy and squared-velocity guesses.
 */
static ghl_error_codes_t tabulated_initialize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      harm_aux_vars_struct *restrict harm_aux,
      double *restrict Z_ptr,
      double *restrict vsq_ptr) {

  if(!isfinite(cons_undens->rho) || cons_undens->rho <= 0.0) {
    return ghl_error_neg_rho;
  }

  ghl_error_codes_t error = ghl_initialize_Noble(
        params, eos, metric_adm, metric_aux, cons_undens, prims, harm_aux, Z_ptr);
  if(error) {
    return error;
  }

  const double v_valenciaU[3] = {
    metric_adm->lapseinv * (prims->vU[0] + metric_adm->betaU[0]),
    metric_adm->lapseinv * (prims->vU[1] + metric_adm->betaU[1]),
    metric_adm->lapseinv * (prims->vU[2] + metric_adm->betaU[2]),
  };
  double vsq = ghl_compute_vec2_from_vec3D(metric_adm->gammaDD, v_valenciaU);
  if(!isfinite(vsq)) {
    return ghl_error_invalid_utsq;
  }
  vsq = ghl_clamp(vsq, 0.0, 1.0 - 1e-15);
  const double Wsq = 1.0 / (1.0 - vsq);
  const double W = sqrt(Wsq);

  double rho = cons_undens->rho / W;
  double Y_e = cons_undens->Y_e / cons_undens->rho;
  double T = prims->temperature;
  if(!isfinite(rho) || !isfinite(Y_e) || !isfinite(T)) {
    return ghl_error_table_bisection;
  }
  ghl_tabulated_enforce_bounds_rho_Ye_T(eos, &rho, &Y_e, &T);

  double press = 0.0;
  double eps = 0.0;
  error = ghl_tabulated_compute_P_eps_from_T(eos, rho, Y_e, T, &press, &eps);
  if(error != ghl_success || !isfinite(press) || !isfinite(eps)) {
    return ghl_error_table_bisection;
  }

  harm_aux->Y_e = Y_e;
  harm_aux->temp_guess = T;

  const double w = rho * (1.0 + eps) + press;
  *Z_ptr = w * Wsq;
  *vsq_ptr = vsq;
  if(!isfinite(*Z_ptr) || *Z_ptr <= 0.0) {
    return ghl_error_invalid_Z;
  }

  return ghl_success;
}

/**
 * @ingroup c2p_internal
 * @brief Evaluates the tabulated Noble residual and Newton step.
 */
static void tabulated_func_2D(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double dummy,
      const double x[],
      double dx[],
      double *restrict f,
      double *restrict df) {
  (void)dummy;

  const double Z = x[0];
  const double vsq = x[1];

  // Keep the shared callback signature; report recoverable failures in harm_aux.
  dx[0] = dx[1] = *f = *df = NAN;
  if(!isfinite(Z) || Z <= 0.0 || !isfinite(vsq) || vsq < 0.0 || vsq >= 1.0) {
    harm_aux->error = ghl_error_c2p_singular;
    return;
  }
  const double Wsq = 1.0 / (1.0 - vsq);
  const double W = sqrt(Wsq);
  const double rho = harm_aux->D / W;
  const double Y_e = harm_aux->Y_e;
  const double h = fabs(Z / (rho * Wsq));

  if(!isfinite(rho) || rho <= 0.0 || !isfinite(Y_e) || !isfinite(h) || h <= 0.0
     || !isfinite(harm_aux->temp_guess)) {
    harm_aux->error = ghl_error_c2p_singular;
    return;
  }

  double press = 0.0;
  double eps = 0.0;
  double dPdrho = 0.0;
  double dPdeps = 0.0;
  double T = harm_aux->temp_guess;
  const ghl_error_codes_t error = ghl_tabulated_compute_P_eps_dPdrho_dPdeps_T_from_h(
        eos, rho, Y_e, h, &press, &eps, &dPdrho, &dPdeps, &T);

  if(error != ghl_success) {
    harm_aux->error = error;
    return;
  }
  if(!isfinite(press) || !isfinite(eps) || !isfinite(dPdrho) || !isfinite(dPdeps)
     || !isfinite(T)) {
    harm_aux->error = ghl_error_table_bisection;
    return;
  }
  harm_aux->temp_guess = T;
  const double dPdeps_over_rho = dPdeps / rho;

  const double dPdZ = (dPdeps_over_rho / (1.0 + dPdeps_over_rho)) / Wsq;

  const double dPdvsq_1 = -0.5 * harm_aux->D * W * dPdrho;
  const double dPdvsq_2 = -0.5 * (Z + press * Wsq) / rho;
  const double dPdvsq = (dPdvsq_1 + dPdeps * dPdvsq_2) / (1.0 + dPdeps_over_rho);

  const double t2 = -0.5 * harm_aux->Bsq + dPdvsq;
  const double t3 = harm_aux->Bsq + Z;
  const double t4 = t3 * t3;
  const double t9 = 1.0 / (Z * Z);
  const double t11
        = harm_aux->Qtsq - vsq * t4 + harm_aux->QdotBsq * (harm_aux->Bsq + 2.0 * Z) * t9;
  const double t16 = harm_aux->QdotBsq * t9;
  const double t18
        = -harm_aux->Qdotn - 0.5 * harm_aux->Bsq * (1.0 + vsq) + 0.5 * t16 - Z + press;
  const double t21 = 1.0 / t3;
  const double t23 = 1.0 / Z;
  const double t24 = t16 * t23;
  const double t25 = -1.0 + dPdZ - t24;
  const double t35 = t25 * t3
                     + (harm_aux->Bsq - 2.0 * dPdvsq)
                             * (harm_aux->QdotBsq + vsq * Z * Z * Z) * t9 * t23;
  const double t36 = 1.0 / t35;
  const double t40 = (vsq + t24) * t3;

  dx[0] = -(t2 * t11 + t4 * t18) * t21 * t36;
  dx[1] = -(-t25 * t11 - 2.0 * t40 * t18) * t21 * t36;

  *df = -t11 * t11 - t18 * t18;
  *f = -0.5 * (*df);
}

/**
 * @ingroup c2p_internal
 * @brief Recovers bounded thermodynamics using the final, speed-limited velocity.
 */
static ghl_error_codes_t tabulated_finalize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const harm_aux_vars_struct *restrict harm_aux,
      const double Z,
      const double vsq,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {
  (void)metric_aux;

  const double gtmp = sqrt(1.0 - vsq);
  double W = 1.0 / gtmp;
  double w = Z * (1.0 - vsq);

  const double nU[4] = {
    metric_adm->lapseinv,
    -metric_adm->lapseinv * metric_adm->betaU[0],
    -metric_adm->lapseinv * metric_adm->betaU[1],
    -metric_adm->lapseinv * metric_adm->betaU[2],
  };

  const double g_o_ZBsq = W / (Z + harm_aux->Bsq);
  const double QdB_o_Z = harm_aux->QdotB / Z;
  const double Qtcon[4] = {
    harm_aux->QU[0] + nU[0] * harm_aux->Qdotn,
    harm_aux->QU[1] + nU[1] * harm_aux->Qdotn,
    harm_aux->QU[2] + nU[2] * harm_aux->Qdotn,
    harm_aux->QU[3] + nU[3] * harm_aux->Qdotn,
  };
  double utU[3] = {
    g_o_ZBsq * (Qtcon[1] + QdB_o_Z * prims->BU[0]),
    g_o_ZBsq * (Qtcon[2] + QdB_o_Z * prims->BU[1]),
    g_o_ZBsq * (Qtcon[3] + QdB_o_Z * prims->BU[2]),
  };

  const bool speed_limited
        = ghl_limit_utilde_and_compute_v(params, metric_adm, utU, prims);
  diagnostics->speed_limited = speed_limited;

  if(speed_limited) {
    // Recompute W so it's compatible with new velocity
    W = metric_adm->lapse * prims->u0;
    w = Z / (W * W);
  }
  if(!isfinite(W) || W <= 0.0 || !isfinite(prims->u0) || prims->u0 <= 0.0
     || !isfinite(prims->vU[0]) || !isfinite(prims->vU[1]) || !isfinite(prims->vU[2])) {
    return ghl_error_invalid_utsq;
  }

  prims->rho = harm_aux->D / W;
  prims->Y_e = cons_undens->Y_e / cons_undens->rho;

  double h = w / prims->rho;
  double press = 0.0;
  double eps = 0.0;
  double S = 0.0;
  double T = harm_aux->temp_guess;
  ghl_error_codes_t error = enforce_bounds_rho_Ye_h(eos, &prims->rho, &prims->Y_e, &h);
  if(error != ghl_success) {
    return error;
  }
  error = ghl_tabulated_compute_P_eps_S_T_from_h(
        eos, prims->rho, prims->Y_e, h, &press, &eps, &S, &T);
  if(error != ghl_success || !isfinite(press) || !isfinite(eps) || !isfinite(S)
     || !isfinite(T)) {
    return ghl_error_table_bisection;
  }

  // The enthalpy inversion above is only restricted by the *table* bounds, so
  // it may return a temperature outside the *runtime* bounds [T_min, T_max]
  // (e.g., the table minimum temperature when h is below the tabulated value
  // at T_min). Enforce the runtime bounds and, if the temperature changed,
  // recompute the thermodynamic quantities so they remain consistent with the
  // bounded temperature. This mirrors the finalization of the Palenzuela and
  // Newman routines, which call ghl_tabulated_enforce_bounds_rho_Ye_T before
  // computing the final pressure and specific internal energy.
  const double T_bounded = ghl_clamp(T, eos->T_min, eos->T_max);
  if(T_bounded != T) {
    T = T_bounded;
    const ghl_error_codes_t T_error = ghl_tabulated_compute_P_eps_S_from_T(
          eos, prims->rho, prims->Y_e, T, &press, &eps, &S);
    if(T_error != ghl_success || !isfinite(press) || !isfinite(eps) || !isfinite(S)) {
      return ghl_error_table_bisection;
    }
  }

  prims->press = press;
  prims->eps = eps;
  prims->entropy = S;
  prims->temperature = T;

  return ghl_success;
}

/**
 * @ingroup tab_c2p
 * @brief Recovers tabulated-EOS primitives with the two-dimensional Noble solver.
 *
 * @details Solves for \f$ Z = \rho h W^2 \f$ and the squared Valencia velocity.
 * After speed limiting, thermodynamics use the final Lorentz factor. Density,
 * composition, and enthalpy obey the EOS bounds; temperature also obeys the
 * runtime bounds. Invalid numerical states return errors for caller fallback.
 *
 * @param[in] params recovery settings
 * @param[in] eos initialized tabulated EOS
 * @param[in] metric_adm ADM metric data
 * @param[in] metric_aux auxiliary four-metric data
 * @param[in] cons_undens undensitized conservative variables
 * @param[in,out] prims primitive guess, magnetic field, and recovered state
 * @param[in,out] diagnostics recovery diagnostics
 * @return ghl_success on recovery, or a recoverable numerical/EOS error.
 */
ghl_error_codes_t ghl_tabulated_Noble2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  diagnostics->speed_limited = false;

  double gnr_out[2];
  harm_aux_vars_struct harm_aux;

  ghl_error_codes_t error = tabulated_initialize_Noble(
        params, eos, metric_adm, metric_aux, cons_undens, prims, &harm_aux, &gnr_out[0],
        &gnr_out[1]);
  if(error) {
    return error;
  }

  const ghl_error_codes_t retval = ghl_general_newton_raphson(
        eos, &harm_aux, 2, 0.0, gnr_out, ghl_validate_2D, tabulated_func_2D);
  const double Z = gnr_out[0];

  if(retval != ghl_success) {
    return retval;
  }
  if(!isfinite(Z) || Z <= 0.0 || Z > 1e20) {
    return ghl_error_invalid_Z;
  }

  double vsq = gnr_out[1];
  if(!isfinite(vsq)) {
    return ghl_error_invalid_utsq;
  }
  if(vsq >= 1.0) {
    vsq = 1.0 - 2.e-16;
  }
  else if(vsq < 0.0) {
    return ghl_error_neg_vsq;
  }

  error = tabulated_finalize_Noble(
        params, eos, metric_adm, metric_aux, cons_undens, &harm_aux, Z, vsq, prims,
        diagnostics);
  if(error) {
    return error;
  }
  if(!ghl_Noble_pressure_is_valid(prims->press)) {
    return ghl_error_neg_pressure;
  }

  diagnostics->n_iter = harm_aux.n_iter;
  diagnostics->which_routine = ghl_con2prim_id_Noble2D;

  return ghl_success;
}
