#include "ghl_m1.h"
#include "ghl_m1_utils.h"
#include <float.h>

static unsigned long long closure_counters[6];

/*
 * The low and high 32-bit words hold the failure stage and validation
 * reason, respectively.  Keeping them in one atomic object prevents a
 * getter from observing a torn pair.  Concurrent closure calls may still
 * replace the process-wide diagnostic with one another; the snapshot is
 * therefore a race-free latest-published diagnostic, not per-call state.
 */
#define CLOSURE_FAILURE_COMPONENT_MASK 0xffffffffULL
#define CLOSURE_FAILURE_REASON_SHIFT 32
static unsigned long long closure_failure_snapshot = 0ULL;

static unsigned long long load_closure_failure_snapshot(void) {
  return __atomic_load_n(&closure_failure_snapshot, __ATOMIC_SEQ_CST);
}

static void update_closure_failure_snapshot(
      const bool update_stage, const unsigned int value) {
  unsigned long long observed = load_closure_failure_snapshot();
  do {
    const unsigned long long value_bits = (unsigned long long)value;
    const unsigned long long desired = update_stage
        ? (observed & (CLOSURE_FAILURE_COMPONENT_MASK <<
                       CLOSURE_FAILURE_REASON_SHIFT)) | value_bits
        : (observed & CLOSURE_FAILURE_COMPONENT_MASK) |
          (value_bits << CLOSURE_FAILURE_REASON_SHIFT);
    if(__atomic_compare_exchange_n(&closure_failure_snapshot, &observed,
                                   desired, false, __ATOMIC_SEQ_CST,
                                   __ATOMIC_SEQ_CST))
      return;
  } while(true);
}

static void record_closure_failure_stage(
      const ghl_m1_closure_failure_stage_t stage) {
  update_closure_failure_snapshot(true, (unsigned int)stage);
}

static void increment_counter(const int index) {
  __atomic_fetch_add(&closure_counters[index], 1ULL, __ATOMIC_RELAXED);
}

void ghl_m1_reset_closure_counters(void) {
  /* This is not a stop-the-world barrier: a concurrent closure call may
   * publish a diagnostic after the reset has completed. */
  for(int i = 0; i < 6; ++i)
    __atomic_store_n(&closure_counters[i], 0ULL, __ATOMIC_SEQ_CST);
  __atomic_store_n(&closure_failure_snapshot, 0ULL, __ATOMIC_SEQ_CST);
}

void ghl_m1_get_closure_counters(
      ghl_m1_closure_counters *restrict counters) {
  if(counters == NULL)
    return;
  counters->ordinary_convergence =
      __atomic_load_n(&closure_counters[0], __ATOMIC_SEQ_CST);
  counters->endpoint_fallback =
      __atomic_load_n(&closure_counters[1], __ATOMIC_SEQ_CST);
  counters->iteration_exhaustion =
      __atomic_load_n(&closure_counters[2], __ATOMIC_SEQ_CST);
  counters->invalid_state =
      __atomic_load_n(&closure_counters[3], __ATOMIC_SEQ_CST);
  counters->downstream_repair =
      __atomic_load_n(&closure_counters[4], __ATOMIC_SEQ_CST);
  counters->residual_rejection =
      __atomic_load_n(&closure_counters[5], __ATOMIC_SEQ_CST);
}

void ghl_m1_get_last_closure_failure_stage(
      ghl_m1_closure_failure_stage_t *restrict stage) {
  if(stage != NULL)
    *stage = (ghl_m1_closure_failure_stage_t)(
        load_closure_failure_snapshot() & CLOSURE_FAILURE_COMPONENT_MASK);
}

void ghl_m1_record_closure_validation_failure(const int reason) {
  update_closure_failure_snapshot(false, (unsigned int)reason);
}

void ghl_m1_get_last_closure_validation_reason(int *restrict reason) {
  if(reason != NULL)
    *reason = (int)(unsigned int)(load_closure_failure_snapshot() >>
                                  CLOSURE_FAILURE_REASON_SHIFT);
}

/* Internal hook used by the realizability repair implementation. */
void ghl_m1_record_closure_downstream_repair(void) {
  increment_counter(4);
}

typedef struct minerbo_workspace {
  const ghl_metric_quantities *metric;
  const ghl_primitive_quantities *prims;
  const ghl_m1_rad_state *rad_state;
  double VU[3];
  double VD[3];
  double W;
  double Pthin[3][3];
  double Pthick[3][3];
  double g_dd[4][4];
  double g_uu[4][4];
  double n_d[4];
  double u_u[4];
  double v_d[4];
  double F_d[4];
  double proj_ud[4][4];
  double Pthin_dd[4][4];
  double Pthick_dd[4][4];
} minerbo_workspace;

static double minerbo_chi(const double xi) {
  const double xi2 = xi * xi;
  return 1.0 / 3.0 + xi2 * (6.0 - 2.0 * xi + 6.0 * xi2) / 15.0;
}

static ghl_error_codes_t build_minerbo_workspace(
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      minerbo_workspace *restrict ws) {
  ws->metric = metric;
  ws->prims = prims;
  ws->rad_state = rad_state;
  ghl_error_codes_t error = ghl_m1_compute_eulerian_velocity(
      metric, prims, ws->VU, ws->VD, &ws->W);
  if(error != ghl_success)
    return error;

  double betaD[3];
  ghl_raise_lower_vector_3D(metric->gammaDD, metric->betaU, betaD);
  double beta2 = 0.0;
  for(int i = 0; i < 3; ++i)
    beta2 += betaD[i] * metric->betaU[i];

  const double alpha = metric->lapse;
  ws->g_dd[0][0] = -alpha * alpha + beta2;
  ws->g_uu[0][0] = -1.0 / (alpha * alpha);
  for(int i = 0; i < 3; ++i) {
    ws->g_dd[0][i+1] = betaD[i];
    ws->g_dd[i+1][0] = betaD[i];
    ws->g_uu[0][i+1] = metric->betaU[i] / (alpha * alpha);
    ws->g_uu[i+1][0] = ws->g_uu[0][i+1];
    for(int j = 0; j < 3; ++j) {
      ws->g_dd[i+1][j+1] = metric->gammaDD[i][j];
      ws->g_uu[i+1][j+1] = metric->gammaUU[i][j]
                         - metric->betaU[i] * metric->betaU[j]
                           / (alpha * alpha);
    }
  }

  ws->n_d[0] = -alpha;
  for(int i = 1; i < 4; ++i) ws->n_d[i] = 0.0;
  ws->u_u[0] = ws->W / alpha;
  ws->v_d[0] = 0.0;
  for(int i = 0; i < 3; ++i) {
    ws->u_u[i+1] = ws->W * (ws->VU[i] - metric->betaU[i] / alpha);
    ws->v_d[i+1] = ws->VD[i];
    ws->F_d[i+1] = rad_state->F[i];
    ws->v_d[0] += betaD[i] * ws->VU[i];
  }
  /* F_mu is Eulerian-spatial, so n^mu F_mu = 0 fixes the time component:
   * F_0 = beta^i F_i.  Establish it before any spacetime norm is formed. */
  ws->F_d[0] = 0.0;
  for(int i = 0; i < 3; ++i)
    ws->F_d[0] += metric->betaU[i] * ws->F_d[i+1];
  if(!isfinite(ws->F_d[0])) {
    record_closure_failure_stage(ghl_m1_closure_failure_workspace);
    return ghl_error_m1_invalid_state;
  }

  double u_d[4] = {0.0};
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu)
      u_d[mu] += ws->g_dd[mu][nu] * ws->u_u[nu];
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu)
      ws->proj_ud[mu][nu] = (mu == nu ? 1.0 : 0.0)
                          + ws->u_u[mu] * u_d[nu];

  long double F2_ld = 0.0L;
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu)
      F2_ld += (long double)ws->g_uu[mu][nu]
             * ws->F_d[mu] * ws->F_d[nu];
  const double F2 = (double)F2_ld;
  if(!isfinite(F2) || F2 < 0.0) {
    record_closure_failure_stage(ghl_m1_closure_failure_workspace);
    return ghl_error_m1_invalid_state;
  }

  /* At exact four-dimensional zero flux, the dyadic thin tensor is zero.
   * Do not replace it with an isotropic spatial tensor. */
  if(F2 == 0.0) {
    for(int mu = 0; mu < 4; ++mu)
      for(int nu = 0; nu < 4; ++nu)
        ws->Pthin_dd[mu][nu] = 0.0;
  } else {
    const double thin_scale = rad_state->E / F2;
    for(int mu = 0; mu < 4; ++mu)
      for(int nu = 0; nu < 4; ++nu)
        ws->Pthin_dd[mu][nu] =
            thin_scale * ws->F_d[mu] * ws->F_d[nu];
  }
  for(int mu = 0; mu < 4; ++mu) {
    for(int nu = 0; nu < 4; ++nu) {
      ws->Pthick_dd[mu][nu] = 0.0;
    }
  }

  double v_dot_F = 0.0;
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu)
      v_dot_F += ws->g_uu[mu][nu] * ws->v_d[mu] * ws->F_d[nu];
  const double W2 = ws->W * ws->W;
  const double coef = 1.0 / (2.0 * W2 + 1.0);
  const double Jo3 = coef * ((2.0 * W2 - 1.0) * rad_state->E
                            - 2.0 * W2 * v_dot_F);
  double tH_d[4];
  for(int mu = 0; mu < 4; ++mu)
    tH_d[mu] = ws->F_d[mu] / ws->W
             + coef * ws->W * ws->v_d[mu]
               * ((4.0 * W2 + 1.0) * v_dot_F
                  - 4.0 * W2 * rad_state->E);
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu) {
      ws->Pthick_dd[mu][nu] = Jo3
          * (4.0 * W2 * ws->v_d[mu] * ws->v_d[nu]
             + ws->g_dd[mu][nu] + ws->n_d[mu] * ws->n_d[nu])
          + ws->W * (tH_d[mu] * ws->v_d[nu]
                     + tH_d[nu] * ws->v_d[mu]);
      if(!isfinite(ws->Pthin_dd[mu][nu]) ||
         !isfinite(ws->Pthick_dd[mu][nu])) {
        record_closure_failure_stage(ghl_m1_closure_failure_workspace);
        return ghl_error_m1_invalid_state;
      }
    }

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) {
      ws->Pthin[i][j] = 0.0;
      ws->Pthick[i][j] = 0.0;
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < 3; ++l) {
          ws->Pthin[i][j] += metric->gammaUU[i][k]
                            * metric->gammaUU[j][l]
                            * ws->Pthin_dd[k+1][l+1];
          ws->Pthick[i][j] += metric->gammaUU[i][k]
                             * metric->gammaUU[j][l]
                             * ws->Pthick_dd[k+1][l+1];
        }
    }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_minerbo_decomposition(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      double Pthin[3][3],
      double Pthick[3][3]) {
  if(m1_params == NULL || metric == NULL || prims == NULL ||
     rad_state == NULL || Pthin == NULL || Pthick == NULL)
    return ghl_error_m1_null_pointer;
  ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric, rad_state, 128.0, NULL);
  if(error != ghl_success)
    return error;
  minerbo_workspace ws;
  error = build_minerbo_workspace(
      metric, prims, rad_state, &ws);
  if(error != ghl_success)
    return error;
  double thin_candidate[3][3], thick_candidate[3][3];
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      thin_candidate[i][j] = ws.Pthin[i][j];
      thick_candidate[i][j] = ws.Pthick[i][j];
    }
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      Pthin[i][j] = thin_candidate[i][j];
      Pthick[i][j] = thick_candidate[i][j];
    }
  }
  return ghl_success;
}

static ghl_error_codes_t evaluate_minerbo(
      const minerbo_workspace *restrict ws,
      const double xi,
      double P[3][3],
      double *restrict chi_out,
      double *restrict residual,
      double *restrict normalized_residual,
      double *restrict physical_xi_out) {
  if(!isfinite(xi) || xi < 0.0 || xi > 1.0)
    return ghl_error_m1_invalid_state;
  const double chi = minerbo_chi(xi);
  if(!isfinite(chi) || chi < 1.0 / 3.0 || chi > 1.0)
    return ghl_error_m1_invalid_state;
  const double dthin = 0.5 * (3.0 * chi - 1.0);
  const double dthick = 1.5 * (1.0 - chi);

  double Pdd[4][4] = {{0.0}};
  double rTdd[4][4] = {{0.0}};
  double H_d[4] = {0.0};
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu) {
      Pdd[mu][nu] = dthin * ws->Pthin_dd[mu][nu]
                  + dthick * ws->Pthick_dd[mu][nu];
      rTdd[mu][nu] = ws->rad_state->E * ws->n_d[mu] * ws->n_d[nu]
                   + ws->F_d[mu] * ws->n_d[nu]
                   + ws->n_d[mu] * ws->F_d[nu]
                   + Pdd[mu][nu];
    }

  double J = 0.0;
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu)
      J += rTdd[mu][nu] * ws->u_u[mu] * ws->u_u[nu];
  if(!isfinite(J) || J <= 0.0) {
    record_closure_failure_stage(ghl_m1_closure_failure_comoving_energy);
    return ghl_error_m1_invalid_state;
  }

  for(int mu = 0; mu < 4; ++mu)
    for(int alpha = 0; alpha < 4; ++alpha)
      for(int beta = 0; beta < 4; ++beta)
        H_d[mu] -= ws->proj_ud[alpha][mu] * rTdd[alpha][beta]
                 * ws->u_u[beta];
  double H2 = 0.0;
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu)
      H2 += ws->g_uu[mu][nu] * H_d[mu] * H_d[nu];
  const double scale = ghl_m1_max(
      ws->rad_state->E * ws->rad_state->E,
      ghl_m1_max(J * J * xi * xi, fabs(H2)));
  const double h2_tolerance = 1024.0 * DBL_EPSILON * scale;
  if(!isfinite(H2) || H2 < -h2_tolerance || !isfinite(scale) ||
     scale <= 0.0) {
    record_closure_failure_stage(ghl_m1_closure_failure_comoving_flux_norm);
    return ghl_error_m1_invalid_state;
  }
  if(H2 < 0.0) H2 = 0.0;
  *residual = J * J * xi * xi - H2;
  *normalized_residual = fabs(*residual) / scale;
  *physical_xi_out = sqrt(H2) / J;
  *chi_out = chi;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) {
      P[i][j] = 0.0;
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < 3; ++l)
          P[i][j] += ws->metric->gammaUU[i][k]
                    * ws->metric->gammaUU[j][l] * Pdd[k+1][l+1];
    }
  for(int i = 0; i < 3; ++i)
    for(int j = i + 1; j < 3; ++j) {
      const double symmetric = 0.5 * (P[i][j] + P[j][i]);
      P[i][j] = symmetric;
      P[j][i] = symmetric;
    }
  if(!isfinite(*normalized_residual) || !isfinite(*physical_xi_out)) {
    record_closure_failure_stage(ghl_m1_closure_failure_residual);
    return ghl_error_m1_invalid_state;
  }
  return ghl_success;
}

static ghl_error_codes_t build_eulerian_minerbo_pressure(
      const minerbo_workspace *restrict ws,
      double P[3][3],
      double *restrict chi_out) {
  /* This is an admissibility fallback only. The primary closure remains the
   * covariant four-dimensional Minerbo construction above. */
  double flux_factor;
  ghl_error_codes_t error = ghl_m1_scaled_covector_norm_ratio(
      ws->metric->gammaUU, ws->rad_state->F, ws->rad_state->E,
      &flux_factor);
  if(error != ghl_success || !isfinite(flux_factor) ||
     flux_factor < 0.0 || flux_factor > 1.0 + 128.0 * DBL_EPSILON)
    return ghl_error_m1_invalid_state;
  const double chi = minerbo_chi(ghl_m1_min(flux_factor, 1.0));
  if(!isfinite(chi) || chi < 1.0 / 3.0 || chi > 1.0)
    return ghl_error_m1_invalid_state;

  double direction[3] = {0.0, 0.0, 0.0};
  double F_scale = 0.0;
  for(int i = 0; i < 3; ++i)
    F_scale = ghl_m1_max(F_scale, fabs(ws->rad_state->F[i]));
  if(F_scale > 0.0) {
    double F_scaled[3], F_scaled_con[3];
    for(int i = 0; i < 3; ++i)
      F_scaled[i] = ws->rad_state->F[i] / F_scale;
    ghl_raise_lower_vector_3D(ws->metric->gammaUU,
                              F_scaled, F_scaled_con);
    long double scaled_norm_sq = 0.0L;
    for(int i = 0; i < 3; ++i)
      scaled_norm_sq += (long double)F_scaled[i] * F_scaled_con[i];
    if(!isfinite((double)scaled_norm_sq) || scaled_norm_sq <= 0.0L)
      return ghl_error_m1_invalid_state;
    const double scaled_norm = sqrt((double)scaled_norm_sq);
    for(int i = 0; i < 3; ++i)
      direction[i] = F_scaled_con[i] / scaled_norm;
  }

  const double dthin = 0.5 * (3.0 * chi - 1.0);
  const double dthick = 1.5 * (1.0 - chi);
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) {
      P[i][j] = ws->rad_state->E
          * (dthin * direction[i] * direction[j]
             + dthick * ws->metric->gammaUU[i][j] / 3.0);
      if(!isfinite(P[i][j]))
        return ghl_error_m1_invalid_state;
    }
  *chi_out = chi;
  return ghl_success;
}

static ghl_error_codes_t publish_eulerian_minerbo_fallback(
      const ghl_m1_parameters *restrict m1_params,
      const minerbo_workspace *restrict ws,
      ghl_m1_closure *restrict closure) {
  ghl_m1_closure candidate = {0};
  ghl_error_codes_t error = build_eulerian_minerbo_pressure(
      ws, candidate.P, &candidate.chi);
  if(error != ghl_success)
    return error;

  /* The fallback has no scalar four-dimensional root. The endpoint status and
   * compatibility bit make that fact observable to source/update callers. */
  candidate.root_residual = 0.0;
  candidate.root_iterations = 0;
  candidate.solve_status = ghl_m1_closure_solve_endpoint_fallback;
  candidate.four_point_compatibility = false;

  ghl_m1_comoving comoving;
  double V_con[3], V_cov[3], W;
  error = ghl_m1_compute_comoving_moments_validated(
      m1_params, ws->metric, ws->prims, ws->rad_state, &candidate,
      &comoving, V_con, V_cov, &W);
  if(error != ghl_success || !isfinite(comoving.J) || comoving.J <= 0.0)
    return ghl_error_m1_invalid_state;

  long double H2_ld = -(long double)comoving.Hn * comoving.Hn;
  for(int i = 0; i < 3; ++i)
    H2_ld += (long double)comoving.HD[i] * comoving.HU[i];
  const double H2 = (double)H2_ld;
  const double scale = ghl_m1_max(comoving.J * comoving.J, fabs(H2));
  const double tolerance = 1024.0 * DBL_EPSILON * scale;
  if(!isfinite(H2) || !isfinite(scale) || scale <= 0.0 ||
     H2 < -tolerance)
    return ghl_error_m1_invalid_state;
  const double xi = sqrt(ghl_m1_max(H2, 0.0)) / comoving.J;
  if(!isfinite(xi) || xi < 0.0 || xi > 1.0 + 1024.0 * DBL_EPSILON)
    return ghl_error_m1_invalid_state;
  candidate.xi = ghl_m1_min(ghl_m1_max(xi, 0.0), 1.0);

  error = ghl_m1_validate_closure_tensor(ws->metric, ws->rad_state,
                                         &candidate);
  if(error != ghl_success)
    return error;
  *closure = candidate;
  return ghl_success;
}

static ghl_error_codes_t publish_minerbo(
      const ghl_m1_parameters *restrict m1_params,
      const minerbo_workspace *restrict ws,
      const double xi,
      const int iterations,
      const ghl_m1_closure_solve_status_t status,
      const double residual_tolerance,
      ghl_m1_closure *restrict closure) {
  ghl_m1_closure candidate;
  double signed_residual, physical_xi;
  ghl_error_codes_t error = evaluate_minerbo(
      ws, xi, candidate.P, &candidate.chi, &signed_residual,
      &candidate.root_residual, &physical_xi);
  if(error != ghl_success)
    return error;
  candidate.xi = physical_xi;
  candidate.root_iterations = iterations;
  candidate.solve_status = status;
  candidate.four_point_compatibility = true;
  if(candidate.root_residual > residual_tolerance)
  {
    record_closure_failure_stage(ghl_m1_closure_failure_residual_gate);
    return ghl_error_m1_closure_residual_too_large;
  }
  error = ghl_m1_validate_closure_tensor_psd(ws->metric, &candidate);
  if(error != ghl_success) {
    ghl_m1_record_closure_validation_failure(
        GHL_M1_CLOSURE_VALIDATION_PSD);
    error = publish_eulerian_minerbo_fallback(m1_params, ws, closure);
    if(error == ghl_success)
      return ghl_success;
    record_closure_failure_stage(ghl_m1_closure_failure_tensor_validation);
    return error;
  }
  error = ghl_m1_validate_closure_tensor(ws->metric, ws->rad_state, &candidate);
  if(error != ghl_success) {
    /* With an exactly zero Eulerian flux, the covariant thin dyad is exactly
     * zero.  A moving fluid can then make the primary candidate fail the
     * pressure-trace invariant even though the Eulerian Minerbo tensor is
     * admissible.  Reuse that established admissibility fallback for this
     * exceptional state; no finite-flux cutoff belongs here. */
    if(ws->rad_state->F[0] == 0.0 && ws->rad_state->F[1] == 0.0 &&
       ws->rad_state->F[2] == 0.0) {
      error = publish_eulerian_minerbo_fallback(m1_params, ws, closure);
      if(error == ghl_success)
        return ghl_success;
    }
    record_closure_failure_stage(ghl_m1_closure_failure_tensor_validation);
    return error;
  }
  *closure = candidate;
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_compute_closure_minerbo_internal(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure,
      const bool configuration_validated) {
  if(m1_params == NULL || metric == NULL || prims == NULL ||
     rad_state == NULL || closure == NULL)
    return ghl_error_m1_null_pointer;
  __atomic_store_n(&closure_failure_snapshot, 0ULL, __ATOMIC_SEQ_CST);
  ghl_error_codes_t error = configuration_validated
      ? ghl_m1_validate_realizability_state(m1_params, metric, rad_state, 128.0, NULL)
      : ghl_m1_validate_realizability(m1_params, metric, rad_state, 128.0, NULL);
  if(error != ghl_success) {
    increment_counter(3);
    return error;
  }

  minerbo_workspace ws;
  error = build_minerbo_workspace(
      metric, prims, rad_state, &ws);
  if(error != ghl_success) {
    increment_counter(3);
    return error;
  }

  double Ptmp[3][3], chi0, chi1, g0, g1, nr0, nr1, physical_xi;
  error = evaluate_minerbo(
      &ws, 0.0, Ptmp, &chi0, &g0, &nr0, &physical_xi);
  if(error == ghl_success)
    error = evaluate_minerbo(
        &ws, 1.0, Ptmp, &chi1, &g1, &nr1, &physical_xi);
  if(error != ghl_success) {
    increment_counter(3);
    return error;
  }

  double xi = 0.0;
  int iterations = 0;
  ghl_m1_closure_solve_status_t status;
  const double endpoint_roundoff =
      1024.0 * DBL_EPSILON * ws.rad_state->E * ws.rad_state->E;
  if(g0 == 0.0 || g1 == 0.0 ||
     fabs(g0) <= endpoint_roundoff || fabs(g1) <= endpoint_roundoff) {
    xi = g0 == 0.0 ? 0.0 : 1.0;
    if(fabs(g0) <= endpoint_roundoff)
      xi = 0.0;
    status = ghl_m1_closure_solve_converged;
  } else if(signbit(g0) == signbit(g1)) {
    xi = nr0 <= nr1 ? 0.0 : 1.0;
    status = ghl_m1_closure_solve_endpoint_fallback;
  } else {
    /* Use the canonical bracketed root contract for the full
     * four-dimensional closure. */
    double a = 0.0, b = 1.0, c = 1.0;
    double fa = g0, fb = g1, fc = g1;
    double d = b - a, e = d;
    status = ghl_m1_closure_solve_iteration_exhausted;
    for(iterations = 0;
        iterations < m1_params->closure_root_max_iterations;
        ++iterations) {
      if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
        c = a;
        fc = fa;
        d = b - a;
        e = d;
      }
      if(fabs(fc) < fabs(fb)) {
        /* This is intentionally sequential: after a=b, c=a means that
         * c receives the old b, as required by Brent's invariant. */
        a = b; fa = fb;
        b = c; fb = fc;
        c = a; fc = fa;
      }
      const double tol = 2.0 * DBL_EPSILON * fabs(b)
                       + 0.5 * m1_params->closure_root_tolerance;
      const double midpoint = 0.5 * (c - b);
      if(fabs(midpoint) <= m1_params->closure_root_tolerance || fb == 0.0) {
        xi = b;
        status = ghl_m1_closure_solve_converged;
        ++iterations;
        break;
      }
      if(fabs(e) >= tol && fabs(fa) > fabs(fb)) {
        const double s = fb / fa;
        double p, q;
        if(a == c) {
          p = 2.0 * midpoint * s;
          q = 1.0 - s;
        } else {
          const double q1 = fa / fc;
          const double r = fb / fc;
          p = s * (2.0 * midpoint * q1 * (q1 - r)
                - (b - a) * (r - 1.0));
          q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
        }
        if(p > 0.0) q = -q;
        else p = -p;
        const double accept_bound =
            3.0 * midpoint * q - fabs(tol * q);
        const double interpolation_bound =
            ghl_m1_min(accept_bound, fabs(e * q));
        if(q != 0.0 && 2.0 * p < interpolation_bound) {
          e = d;
          d = p / q;
        } else {
          d = midpoint;
          e = midpoint;
        }
      } else {
        d = midpoint;
        e = midpoint;
      }
      a = b;
      fa = fb;
      b += fabs(d) > tol ? d : copysign(tol, midpoint);
      xi = b;
      double chim, gm, nrm;
      error = evaluate_minerbo(
          &ws, xi, Ptmp, &chim, &gm, &nrm, &physical_xi);
      if(error != ghl_success) {
        increment_counter(3);
        return error;
      }
      fb = gm;
    }
    if(status == ghl_m1_closure_solve_iteration_exhausted)
      xi = b;
  }

  error = publish_minerbo(
      m1_params, &ws, xi, iterations, status,
      m1_params->closure_root_residual_tolerance, closure);
  if(error != ghl_success) {
    if(status == ghl_m1_closure_solve_endpoint_fallback)
      increment_counter(1);
    else if(status == ghl_m1_closure_solve_iteration_exhausted)
      increment_counter(2);
    if(error == ghl_error_m1_closure_residual_too_large) {
      increment_counter(5);
      return error;
    }
    increment_counter(3);
    return error;
  }
  if(closure->solve_status == ghl_m1_closure_solve_converged)
    increment_counter(0);
  else if(closure->solve_status == ghl_m1_closure_solve_endpoint_fallback)
    increment_counter(1);
  else
    increment_counter(2);
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_closure_minerbo(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure) {
  return ghl_m1_compute_closure_minerbo_internal(
      m1_params, metric, prims, rad_state, closure, false);
}

ghl_error_codes_t ghl_m1_compute_closure_minerbo_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure) {
  return ghl_m1_compute_closure_minerbo_internal(
      m1_params, metric, prims, rad_state, closure, true);
}

ghl_error_codes_t ghl_m1_compute_closure_with_primitives(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure) {
  return ghl_m1_compute_closure_minerbo(
      m1_params, metric, prims, rad_state, closure);
}
