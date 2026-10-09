#include "ghl_m1.h"
#include "ghl_m1_closure_private.h"
#include "ghl_m1_utils.h"
#include <float.h>

typedef struct minerbo_workspace {
  const ghl_metric_quantities *metric;
  const ghl_primitive_quantities *prims;
  const ghl_m1_rad_state *rad_state;
  double VU[3];
  double VD[3];
  double W;
  double Pthin[3][3];
  double Pthick[3][3];
  double gDD[4][4];
  double nD[4];
  double vD[4];
  double FD[4];
  double PthinDD[4][4];
  double PthickDD[4][4];
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
  ghl_error_codes_t error
        = ghl_m1_compute_eulerian_velocity(metric, prims, ws->VU, ws->VD, &ws->W);
  if(error != ghl_success) {
    return error;
  }

  double betaD[3];
  ghl_raise_lower_vector_3D(metric->gammaDD, metric->betaU, betaD);
  double beta2 = 0.0;
  for(int i = 0; i < 3; ++i) {
    beta2 += betaD[i] * metric->betaU[i];
  }

  const double alpha = metric->lapse;
  ws->gDD[0][0] = -alpha * alpha + beta2;
  for(int i = 0; i < 3; ++i) {
    ws->gDD[0][i + 1] = betaD[i];
    ws->gDD[i + 1][0] = betaD[i];
    for(int j = 0; j < 3; ++j) {
      ws->gDD[i + 1][j + 1] = metric->gammaDD[i][j];
    }
  }

  ws->nD[0] = -alpha;
  for(int i = 1; i < 4; ++i) {
    ws->nD[i] = 0.0;
  }
  ws->vD[0] = 0.0;
  for(int i = 0; i < 3; ++i) {
    ws->vD[i + 1] = ws->VD[i];
    ws->FD[i + 1] = rad_state->F[i];
    ws->vD[0] += betaD[i] * ws->VU[i];
  }
  /* F_mu is Eulerian-spatial, so n^mu F_mu = 0 fixes the time component:
   * F_0 = beta^i F_i.  Establish it before any spacetime norm is formed. */
  ws->FD[0] = 0.0;
  for(int i = 0; i < 3; ++i) {
    ws->FD[0] += metric->betaU[i] * ws->FD[i + 1];
  }
  if(!isfinite(ws->FD[0])) {
    return ghl_error_m1_invalid_state;
  }

  double F_scale = 0.0;
  for(int mu = 0; mu < 4; ++mu) {
    F_scale = fmax(F_scale, fabs(ws->FD[mu]));
  }

  /* At exact four-dimensional zero flux, the dyadic thin tensor is zero.
   * Test the flux itself so norm underflow cannot select this policy. */
  if(F_scale == 0.0) {
    for(int mu = 0; mu < 4; ++mu) {
      for(int nu = 0; nu < 4; ++nu) {
        ws->PthinDD[mu][nu] = 0.0;
      }
    }
  }
  else {
    /* The thin tensor depends on direction, not flux magnitude. Normalize
     * before squaring to avoid both norm under/overflow and E/F^2 overflow
     * for otherwise representable pressure tensors. */
    double F_shape[4] = { 0.0 };
    for(int mu = 0; mu < 4; ++mu) {
      F_shape[mu] = ws->FD[mu] / F_scale;
    }
    /* F_mu is Eulerian-spatial. Contract only its spatial components with
     * gamma^ij; the equivalent four-dimensional expression subtracts large
     * beta^i beta^j / alpha^2 terms and loses precision at small lapse. */
    double F2_shape_acc = 0.0;
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        F2_shape_acc += metric->gammaUU[i][j] * F_shape[i + 1] * F_shape[j + 1];
      }
    }
    const double F2_shape = F2_shape_acc;
    /* The validated finite SPD metric and nonzero normalized covector make
     * this quadratic form strictly positive and representable. */
    for(int mu = 0; mu < 4; ++mu) {
      for(int nu = 0; nu < 4; ++nu) {
        ws->PthinDD[mu][nu] = rad_state->E * F_shape[mu] * F_shape[nu] / F2_shape;
      }
    }
  }
  for(int mu = 0; mu < 4; ++mu) {
    for(int nu = 0; nu < 4; ++nu) {
      ws->PthickDD[mu][nu] = 0.0;
    }
  }

  const double energy_square_limit = sqrt(DBL_MAX);
  if(rad_state->E <= energy_square_limit) {
    /* v_mu F^mu is V^i F_i for Eulerian-spatial v and F. */
    double v_dot_F_acc = 0.0;
    for(int i = 0; i < 3; ++i) {
      v_dot_F_acc += ws->VU[i] * ws->FD[i + 1];
    }
    const double v_dot_F = v_dot_F_acc;
    const double W2 = ws->W * ws->W;
    const double inverse_W2 = 1.0 / W2;
    const double E_minus_v_dot_F = rad_state->E - v_dot_F;
    /* Divide the thick equations by W^2 before subtracting. Their original
     * O(W^2 E) terms cancel even for boosted isotropic radiation. */
    const double Jo3
          = (2.0 * E_minus_v_dot_F - inverse_W2 * rad_state->E) / (2.0 + inverse_W2);
    double tHD[4];
    for(int mu = 0; mu < 4; ++mu) {
      tHD[mu] = fma(-ws->W * ws->vD[mu], E_minus_v_dot_F + Jo3, ws->FD[mu] / ws->W);
    }
    for(int mu = 0; mu < 4; ++mu) {
      for(int nu = 0; nu < 4; ++nu) {
        ws->PthickDD[mu][nu] = Jo3
                                     * (4.0 * W2 * ws->vD[mu] * ws->vD[nu]
                                        + ws->gDD[mu][nu] + ws->nD[mu] * ws->nD[nu])
                               + ws->W * (tHD[mu] * ws->vD[nu] + tHD[nu] * ws->vD[mu]);
        if(!isfinite(ws->PthinDD[mu][nu]) || !isfinite(ws->PthickDD[mu][nu])) {
          return ghl_error_m1_invalid_state;
        }
      }
    }
  }
  else {
    /* Keep every thick-limit operation normalized by E before restoring the
     * physical pressure scale. This avoids intermediate products such as
     * 4 W^2 E overflowing before they are multiplied by a zero velocity. */
    double F_over_E[4];
    for(int mu = 0; mu < 4; ++mu) {
      /* FD is finite and this branch has E > sqrt(DBL_MAX) > 1. */
      F_over_E[mu] = ws->FD[mu] / rad_state->E;
    }
    double v_dot_F_over_E_acc = 0.0;
    for(int i = 0; i < 3; ++i) {
      v_dot_F_over_E_acc += ws->VU[i] * F_over_E[i + 1];
    }
    const double v_dot_F_over_E = v_dot_F_over_E_acc;
    const double W2 = ws->W * ws->W;
    const double inverse_W2 = 1.0 / W2;
    const double E_minus_v_dot_F_over_E = 1.0 - v_dot_F_over_E;
    const double Jo3_over_E
          = (2.0 * E_minus_v_dot_F_over_E - inverse_W2) / (2.0 + inverse_W2);
    double tHD_over_E[4];
    for(int mu = 0; mu < 4; ++mu) {
      tHD_over_E[mu]
            = fma(-ws->W * ws->vD[mu], E_minus_v_dot_F_over_E + Jo3_over_E,
                  F_over_E[mu] / ws->W);
    }
    for(int mu = 0; mu < 4; ++mu) {
      for(int nu = 0; nu < 4; ++nu) {
        const double Pthick_over_E
              = Jo3_over_E
                      * (4.0 * W2 * ws->vD[mu] * ws->vD[nu] + ws->gDD[mu][nu]
                         + ws->nD[mu] * ws->nD[nu])
                + ws->W * (tHD_over_E[mu] * ws->vD[nu] + tHD_over_E[nu] * ws->vD[mu]);
        ws->PthickDD[mu][nu] = rad_state->E * Pthick_over_E;
        if(!isfinite(ws->PthinDD[mu][nu]) || !isfinite(ws->PthickDD[mu][nu])) {
          return ghl_error_m1_invalid_state;
        }
      }
    }
  }

  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      ws->Pthin[i][j] = 0.0;
      ws->Pthick[i][j] = 0.0;
      for(int k = 0; k < 3; ++k) {
        for(int l = 0; l < 3; ++l) {
          ws->Pthin[i][j] += metric->gammaUU[i][k] * metric->gammaUU[j][l]
                             * ws->PthinDD[k + 1][l + 1];
          ws->Pthick[i][j] += metric->gammaUU[i][k] * metric->gammaUU[j][l]
                              * ws->PthickDD[k + 1][l + 1];
        }
      }
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
  if(m1_params == NULL || metric == NULL || prims == NULL || rad_state == NULL
     || Pthin == NULL || Pthick == NULL) {
    return ghl_error_m1_null_pointer;
  }
  ghl_error_codes_t error
        = ghl_m1_validate_realizability_state(m1_params, metric, rad_state, 128.0, NULL);
  if(error != ghl_success) {
    return error;
  }
  minerbo_workspace ws;
  error = build_minerbo_workspace(metric, prims, rad_state, &ws);
  if(error != ghl_success) {
    return error;
  }
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

ghl_error_codes_t ghl_m1_closure_private_evaluate_invariant(
      const double J,
      const double H2,
      const double scale,
      const double xi,
      double *const H2_clipped,
      double *const residual,
      double *const normalized_residual,
      double *const physical_xi) {
  if(!isfinite(J) || J <= 0.0) {
    return ghl_error_m1_invalid_state;
  }
  const double h2_tolerance = 1024.0 * DBL_EPSILON * scale;
  if(!isfinite(H2) || H2 < -h2_tolerance) {
    return ghl_error_m1_invalid_state;
  }
  const double clipped_H2 = fmax(H2, 0.0);
  if(!isfinite(scale) || !(scale > 0.0)) {
    return ghl_error_m1_invalid_state;
  }

  const double candidate_residual = J * J * xi * xi - clipped_H2;
  /* In production, scale=max(J^2,abs(H2)) and xi is in [0,1], so both
   * residual terms are bounded by scale and this quotient cannot overflow. */
  const double candidate_normalized_residual = fabs(candidate_residual) / scale;
  const double candidate_physical_xi = sqrt(clipped_H2) / J;
  if(!isfinite(candidate_physical_xi)) {
    return ghl_error_m1_invalid_state;
  }
  *H2_clipped = clipped_H2;
  *residual = candidate_residual;
  *normalized_residual = candidate_normalized_residual;
  *physical_xi = candidate_physical_xi;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_closure_private_finish_evaluation(
      const ghl_metric_quantities *const metric,
      const double PDD[4][4],
      const double energy_scale,
      const double W,
      const double J,
      const double H2,
      const double scale,
      const double xi,
      const double chi,
      ghl_m1_closure_evaluation *const evaluation) {
  ghl_m1_closure_evaluation candidate = { 0 };
  double H2_clipped;
  const ghl_error_codes_t error = ghl_m1_closure_private_evaluate_invariant(
        J, H2, scale, xi, &H2_clipped, &candidate.residual,
        &candidate.normalized_residual, &candidate.physical_xi);
  if(error != ghl_success) {
    return error;
  }
  (void)H2_clipped;
  candidate.chi = chi;
  double (*P)[3] = candidate.P;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      P[i][j] = 0.0;
      for(int k = 0; k < 3; ++k) {
        for(int l = 0; l < 3; ++l) {
          P[i][j] += metric->gammaUU[i][k] * metric->gammaUU[j][l] * PDD[k + 1][l + 1];
        }
      }
    }
  }
  /* The thick tensor cancels O(W^2 E) terms to produce its O(E) trace, so
   * its rounding error in gamma_ij P^ij grows like W^2.  Inside that rounding
   * envelope, restore the exact trace with an isotropic correction; a larger
   * discrepancy is left for the tensor validator to reject. */
  double trace_acc = 0.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      trace_acc += metric->gammaDD[i][j] * P[i][j];
    }
  }
  const double trace_error = energy_scale - trace_acc;
  const double trace_envelope = 256.0 * DBL_EPSILON * (1.0 + 4.0 * W * W)
                                * fmax(fabs(trace_acc), energy_scale);
  if(fabs(trace_error) <= trace_envelope) {
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        P[i][j] += (trace_error / 3.0) * metric->gammaUU[i][j];
      }
    }
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = i + 1; j < 3; ++j) {
      const double symmetric = 0.5 * (P[i][j] + P[j][i]);
      P[i][j] = symmetric;
      P[j][i] = symmetric;
    }
  }
  *evaluation = candidate;
  return ghl_success;
}

static ghl_error_codes_t evaluate_minerbo(
      const minerbo_workspace *restrict ws,
      const double xi,
      ghl_m1_closure_evaluation *restrict evaluation) {
  /* Every caller supplies an endpoint or a Brent iterate inside [0,1]. */
  const double chi = minerbo_chi(xi);
  /* minerbo_chi maps that closed interval into [1/3,1]. */
  const double dthin = 0.5 * (3.0 * chi - 1.0);
  const double dthick = 1.5 * (1.0 - chi);

  double PDD[4][4] = { { 0.0 } };
  double rTDD[4][4] = { { 0.0 } };
  const double energy_scale = ws->rad_state->E;
  for(int mu = 0; mu < 4; ++mu) {
    for(int nu = 0; nu < 4; ++nu) {
      PDD[mu][nu] = dthin * ws->PthinDD[mu][nu] + dthick * ws->PthickDD[mu][nu];
      rTDD[mu][nu] = energy_scale * ws->nD[mu] * ws->nD[nu] + ws->FD[mu] * ws->nD[nu]
                     + ws->nD[mu] * ws->FD[nu] + PDD[mu][nu];
      /* A nonfinite PDD component necessarily propagates to rTDD below. */
      if(!isfinite(rTDD[mu][nu])) {
        return ghl_error_m1_invalid_state;
      }
    }
  }

  const bool scaled_energy
        = energy_scale < sqrt(DBL_MIN) || energy_scale > sqrt(DBL_MAX);
  double J = 0.0;
  double H2 = 0.0;
  double scale = 0.0;
  double F_con[3] = { 0.0, 0.0, 0.0 };
  double V_cov[3] = { 0.0, 0.0, 0.0 };
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      F_con[i] += ws->metric->gammaUU[i][j] * ws->FD[j + 1];
      V_cov[i] += ws->metric->gammaDD[i][j] * ws->VU[j];
    }
  }
  if(xi == 0.0) {
    /* Recover the thick moments directly from E/F, not by boosting the
     * rounded pressure back: that contraction amplifies pressure roundoff
     * by W^2 and can miss an otherwise acceptable zero endpoint. Work in
     * units of E so both energy-square limits use the same endpoint rule. */
    double FdotV_over_E = 0.0;
    for(int i = 0; i < 3; ++i) {
      FdotV_over_E += (ws->FD[i + 1] / energy_scale) * ws->VU[i];
    }
    const double inverse_W2 = 1.0 / (ws->W * ws->W);
    const double E_minus_FdotV_over_E = 1.0 - FdotV_over_E;
    const double Jo3_over_E
          = (2.0 * E_minus_FdotV_over_E - inverse_W2) / (2.0 + inverse_W2);
    J = 3.0 * Jo3_over_E;
    double HU_over_E[3];
    for(int i = 0; i < 3; ++i) {
      HU_over_E[i]
            = fma(-ws->W * ws->VU[i], E_minus_FdotV_over_E + Jo3_over_E,
                  F_con[i] / energy_scale / ws->W);
    }
    double Hn_over_E = 0.0;
    for(int i = 0; i < 3; ++i) {
      Hn_over_E -= V_cov[i] * HU_over_E[i];
      for(int j = 0; j < 3; ++j) {
        H2 += ws->metric->gammaDD[i][j] * HU_over_E[i] * HU_over_E[j];
      }
    }
    H2 -= Hn_over_E * Hn_over_E;
    if(!scaled_energy) {
      /* Brent compares signed residuals from different trial points. Keep
       * the same physical units as the non-endpoint evaluations. */
      J *= energy_scale;
      H2 = (H2 * energy_scale) * energy_scale;
    }
    scale = fmax(J * J, fabs(H2));
  }
  else if(!scaled_energy) {
    double FdotV_acc = 0.0;
    double PVV_acc = 0.0;
    for(int i = 0; i < 3; ++i) {
      FdotV_acc += ws->FD[i + 1] * ws->VU[i];
      for(int j = 0; j < 3; ++j) {
        PVV_acc += PDD[i + 1][j + 1] * ws->VU[i] * ws->VU[j];
      }
    }
    const double FdotV = FdotV_acc;
    const double PVV = PVV_acc;
    const double W2 = ws->W * ws->W;
    J = W2 * (energy_scale - 2.0 * FdotV + PVV);

    double HU[3] = { 0.0, 0.0, 0.0 };
    for(int i = 0; i < 3; ++i) {
      double PijVj_acc = 0.0;
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          PijVj_acc += ws->metric->gammaUU[i][k] * PDD[k + 1][j + 1] * ws->VU[j];
        }
      }
      HU[i] = ws->W * (F_con[i] - PijVj_acc - J * ws->VU[i]);
    }
    double Hn_acc = 0.0;
    double H2_acc = 0.0;
    for(int i = 0; i < 3; ++i) {
      Hn_acc -= V_cov[i] * HU[i];
      for(int j = 0; j < 3; ++j) {
        H2_acc += ws->metric->gammaDD[i][j] * HU[i] * HU[j];
      }
    }
    H2_acc -= Hn_acc * Hn_acc;
    H2 = H2_acc;
    /* The residual J^2 xi^2 - H^2 is bounded by J^2 (H <= J), and J can be
     * far below E for flux along a fast flow, so J^2 is the natural scale. */
    scale = fmax(J * J, fabs(H2));
  }
  else {
    /* The comoving consistency equation is homogeneous in radiation energy.
     * Evaluate it after dividing the stress tensor by E so its residual scale
     * remains representable when E^2 underflows or overflows. */
    for(int i = 0; i < 3; ++i) {
      F_con[i] /= energy_scale;
    }
    double FdotV_over_E_acc = 0.0;
    double PVV_over_E_acc = 0.0;
    for(int i = 0; i < 3; ++i) {
      FdotV_over_E_acc += F_con[i] * V_cov[i];
      for(int j = 0; j < 3; ++j) {
        PVV_over_E_acc += (PDD[i + 1][j + 1] / energy_scale) * ws->VU[i] * ws->VU[j];
      }
    }
    const double W2 = ws->W * ws->W;
    J = W2 * (1.0 - 2.0 * FdotV_over_E_acc + PVV_over_E_acc);

    double HU_over_E[3] = { 0.0, 0.0, 0.0 };
    for(int i = 0; i < 3; ++i) {
      double PijVj_over_E_acc = 0.0;
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          PijVj_over_E_acc += ws->metric->gammaUU[i][k]
                              * (PDD[k + 1][j + 1] / energy_scale) * ws->VU[j];
        }
      }
      HU_over_E[i] = ws->W * (F_con[i] - PijVj_over_E_acc - J * ws->VU[i]);
    }
    double Hn_over_E_acc = 0.0;
    double H2_over_E2_acc = 0.0;
    for(int i = 0; i < 3; ++i) {
      Hn_over_E_acc -= V_cov[i] * HU_over_E[i];
      for(int j = 0; j < 3; ++j) {
        H2_over_E2_acc += ws->metric->gammaDD[i][j] * HU_over_E[i] * HU_over_E[j];
      }
    }
    H2_over_E2_acc -= Hn_over_E_acc * Hn_over_E_acc;
    H2 = H2_over_E2_acc;
    scale = fmax(J * J, fabs(H2));
  }
  return ghl_m1_closure_private_finish_evaluation(
        ws->metric, PDD, energy_scale, ws->W, J, H2, scale, xi, chi, evaluation);
}

static ghl_error_codes_t build_eulerian_minerbo_pressure(
      const minerbo_workspace *restrict ws,
      double P[3][3],
      double *restrict chi_out) {
  /* This is an admissibility fallback only. The primary closure remains the
   * covariant four-dimensional Minerbo construction above. */
  double flux_factor = 0.0;
  /* The public closure entry validated this metric and realizability state.
   * The scaled norm therefore succeeds and is bounded by the tolerated cone. */
  (void)ghl_m1_scaled_covector_norm_ratio(
        ws->metric->gammaUU, ws->rad_state->F, ws->rad_state->E, &flux_factor);

  double direction[3] = { 0.0, 0.0, 0.0 };
  double F_scale = 0.0;
  for(int i = 0; i < 3; ++i) {
    F_scale = fmax(F_scale, fabs(ws->rad_state->F[i]));
  }
  if(F_scale > 0.0) {
    double F_scaled[3], F_scaled_con[3];
    for(int i = 0; i < 3; ++i) {
      F_scaled[i] = ws->rad_state->F[i] / F_scale;
    }
    ghl_raise_lower_vector_3D(ws->metric->gammaUU, F_scaled, F_scaled_con);
    double scaled_norm_sq = 0.0;
    for(int i = 0; i < 3; ++i) {
      scaled_norm_sq += F_scaled[i] * F_scaled_con[i];
    }
    /* The validated finite SPD metric makes this normalized norm positive. */
    const double scaled_norm = sqrt(scaled_norm_sq);
    for(int i = 0; i < 3; ++i) {
      direction[i] = F_scaled_con[i] / scaled_norm;
    }
  }

  return ghl_m1_closure_private_construct_eulerian_pressure(
        ws->rad_state->E, ws->metric->gammaUU, flux_factor, direction, P, chi_out);
}

ghl_error_codes_t ghl_m1_closure_private_construct_eulerian_pressure(
      const double energy,
      const double gammaUU[3][3],
      const double flux_factor,
      const double direction[3],
      double P[3][3],
      double *const chi_out) {
  const double chi = minerbo_chi(fmin(flux_factor, 1.0));
  const double dthin = 0.5 * (3.0 * chi - 1.0);
  const double dthick = 1.5 * (1.0 - chi);
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      /* Parenthesize the dyad. Multiplication associates left to right, so
       * dthin*d[i]*d[j] evaluates as (dthin*d[i])*d[j], which does not equal
       * (dthin*d[j])*d[i]. Grouping the commutative factor first keeps this
       * term bit-identical under an index swap. */
      P[i][j] = energy
                * (dthin * (direction[i] * direction[j]) + dthick * gammaUU[i][j] / 3.0);
      if(!isfinite(P[i][j])) {
        return ghl_error_m1_invalid_state;
      }
    }
  }
  /* Symmetrize exactly, as the primary construction in evaluate_minerbo does.
   * The grouping above removes the association-order asymmetry, but the metric
   * validator accepts a gammaUU that is symmetric only to its 64*DBL_EPSILON
   * bound, so a caller metric can still leak asymmetry into the thick term.
   * The published tensor must satisfy the symmetry invariant exactly. */
  for(int i = 0; i < 3; ++i) {
    for(int j = i + 1; j < 3; ++j) {
      const double symmetric = 0.5 * (P[i][j] + P[j][i]);
      P[i][j] = symmetric;
      P[j][i] = symmetric;
    }
  }
  *chi_out = chi;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_closure_private_finalize_fallback(
      const ghl_error_codes_t conversion_error,
      const ghl_m1_comoving *const comoving,
      const ghl_m1_closure *const candidate,
      ghl_m1_closure *const closure) {
  if(conversion_error != ghl_success || comoving->J <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  double xi;
  if(comoving->J >= sqrt(DBL_MIN) && comoving->J <= sqrt(DBL_MAX)) {
    double H2_acc = -comoving->Hn * comoving->Hn;
    for(int i = 0; i < 3; ++i) {
      H2_acc += comoving->HD[i] * comoving->HU[i];
    }
    const double H2 = H2_acc;
    const double scale = fmax(comoving->J * comoving->J, fabs(H2));
    const double tolerance = 1024.0 * DBL_EPSILON * scale;
    /* The bounded J range and finite H2 imply a finite, positive scale. */
    if(!isfinite(H2) || H2 < -tolerance) {
      return ghl_error_m1_invalid_state;
    }
    xi = sqrt(fmax(H2, 0.0)) / comoving->J;
  }
  else {
    /* The fallback uses the same homogeneous normalization as the primary
     * residual, avoiding J^2 underflow or overflow in its admissibility check. */
    const double Hn_over_J = comoving->Hn / comoving->J;
    double H2_scaled_acc = -Hn_over_J * Hn_over_J;
    for(int i = 0; i < 3; ++i) {
      const double HD_over_J = comoving->HD[i] / comoving->J;
      const double HU_over_J = comoving->HU[i] / comoving->J;
      H2_scaled_acc += HD_over_J * HU_over_J;
    }
    const double H2_scaled = H2_scaled_acc;
    const double scale = fmax(1.0, fabs(H2_scaled));
    const double tolerance = 1024.0 * DBL_EPSILON * scale;
    /* Once H2_scaled is finite, max(1,abs(H2_scaled)) is finite and positive. */
    if(!isfinite(Hn_over_J) || !isfinite(H2_scaled) || H2_scaled < -tolerance) {
      return ghl_error_m1_invalid_state;
    }
    xi = sqrt(fmax(H2_scaled, 0.0));
  }
  /* The two guarded square-root paths produce a finite xi: ordinary J is at
   * least sqrt(DBL_MIN), while scaled xi is the square root of finite H2. */
  if(xi > 1.0 + 1024.0 * DBL_EPSILON) {
    return ghl_error_m1_invalid_state;
  }

  /* Keep failed fallback construction transactional for the public caller. */
  ghl_m1_closure result = *candidate;
  result.xi = fmin(xi, 1.0);
  *closure = result;
  return ghl_success;
}

static ghl_error_codes_t publish_eulerian_minerbo_fallback(
      const ghl_m1_parameters *restrict m1_params,
      const minerbo_workspace *restrict ws,
      ghl_m1_closure *restrict closure) {
  ghl_m1_closure candidate = { 0 };
  ghl_error_codes_t error
        = build_eulerian_minerbo_pressure(ws, candidate.P, &candidate.chi);
  if(error != ghl_success) {
    return error;
  }

  /* The fallback has no scalar four-dimensional root. The endpoint status and
   * compatibility bit make that fact observable to source/update callers. */
  candidate.root_residual = 0.0;
  candidate.root_iterations = 0;
  candidate.solve_status = ghl_m1_closure_solve_endpoint_fallback;
  candidate.four_point_compatibility = false;

  ghl_m1_comoving comoving;
  double V_con[3], V_cov[3], W;
  error = ghl_m1_compute_comoving_moments_validated(
        m1_params, ws->metric, ws->prims, ws->rad_state, &candidate, &comoving, V_con,
        V_cov, &W);
  /* The comoving-moment call already validated this unchanged tensor. */
  return ghl_m1_closure_private_finalize_fallback(error, &comoving, &candidate, closure);
}

ghl_error_codes_t ghl_m1_closure_private_check_residual_gate(
      const double normalized_residual,
      const double residual_tolerance) {
  if(normalized_residual > residual_tolerance) {
    return ghl_error_m1_closure_residual_too_large;
  }
  return ghl_success;
}

static ghl_error_codes_t publish_minerbo(
      const ghl_m1_parameters *restrict m1_params,
      const minerbo_workspace *restrict ws,
      const ghl_m1_closure_evaluation *restrict evaluation,
      const int iterations,
      const ghl_m1_closure_solve_status_t status,
      const double residual_tolerance,
      ghl_m1_closure *restrict closure) {
  ghl_m1_closure candidate = { 0 };
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      candidate.P[i][j] = evaluation->P[i][j];
    }
  }
  candidate.chi = evaluation->chi;
  candidate.xi = evaluation->physical_xi;
  candidate.root_residual = evaluation->normalized_residual;
  candidate.root_iterations = iterations;
  candidate.solve_status = status;
  candidate.four_point_compatibility = true;
  ghl_error_codes_t error = ghl_m1_closure_private_check_residual_gate(
        candidate.root_residual, residual_tolerance);
  if(error != ghl_success) {
    return error;
  }
  error = ghl_m1_validate_closure_tensor_psd(ws->metric, &candidate);
  if(error != ghl_success) {
    return publish_eulerian_minerbo_fallback(m1_params, ws, closure);
  }
  error = ghl_m1_validate_closure_tensor(ws->metric, ws->rad_state, &candidate);
  if(error != ghl_success) {
    /* With an exactly zero Eulerian flux, the covariant thin dyad is exactly
     * zero.  A moving fluid can then make the primary candidate fail the
     * pressure-trace invariant even though the Eulerian Minerbo tensor is
     * admissible.  Reuse that established admissibility fallback for this
     * exceptional state; no finite-flux cutoff belongs here. */
    if(ws->rad_state->F[0] == 0.0 && ws->rad_state->F[1] == 0.0
       && ws->rad_state->F[2] == 0.0) {
      error = publish_eulerian_minerbo_fallback(m1_params, ws, closure);
      if(error == ghl_success) {
        return ghl_success;
      }
    }
    return error;
  }
  *closure = candidate;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_closure_private_solve_root(
      const ghl_m1_parameters *const params,
      const ghl_m1_closure_evaluator evaluate,
      void *const context,
      ghl_m1_closure_root_result *const root) {
  ghl_m1_closure_evaluation evaluation0, evaluation1;
  ghl_error_codes_t error = evaluate(context, 0.0, &evaluation0);
  if(error != ghl_success) {
    return error;
  }
  error = evaluate(context, 1.0, &evaluation1);
  if(error != ghl_success) {
    return error;
  }

  ghl_m1_closure_root_result candidate = { 0 };
  /* Endpoint roots are judged on the J^2-normalized residual. */
  const double endpoint_roundoff = 1024.0 * DBL_EPSILON;
  if(evaluation0.normalized_residual <= endpoint_roundoff
     || evaluation1.normalized_residual <= endpoint_roundoff) {
    const bool use_first = evaluation0.normalized_residual <= endpoint_roundoff;
    candidate.xi = use_first ? 0.0 : 1.0;
    candidate.status = ghl_m1_closure_solve_converged;
    candidate.evaluation = use_first ? evaluation0 : evaluation1;
  }
  else if(signbit(evaluation0.residual) == signbit(evaluation1.residual)) {
    const bool use_first
          = evaluation0.normalized_residual <= evaluation1.normalized_residual;
    candidate.xi = use_first ? 0.0 : 1.0;
    candidate.status = ghl_m1_closure_solve_endpoint_fallback;
    candidate.evaluation = use_first ? evaluation0 : evaluation1;
  }
  else {
    /* Use the canonical bracketed root contract for the full
     * four-dimensional closure. */
    double a = 0.0, b = 1.0, c = 1.0;
    double fa = evaluation0.residual;
    double fb = evaluation1.residual;
    double fc = evaluation1.residual;
    ghl_m1_closure_evaluation evaluation_a = evaluation0;
    ghl_m1_closure_evaluation evaluation_b = evaluation1;
    ghl_m1_closure_evaluation evaluation_c = evaluation1;
    double d = b - a, e = d;
    candidate.status = ghl_m1_closure_solve_iteration_exhausted;
    for(candidate.iterations = 0;
        candidate.iterations < params->closure_root_max_iterations;
        ++candidate.iterations) {
      if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
        c = a;
        fc = fa;
        evaluation_c = evaluation_a;
        d = b - a;
        e = d;
      }
      if(fabs(fc) < fabs(fb)) {
        /* Sequential updates preserve Brent's a=b; b=c; c=a invariant. */
        a = b;
        fa = fb;
        evaluation_a = evaluation_b;
        b = c;
        fb = fc;
        evaluation_b = evaluation_c;
        c = a;
        fc = fa;
        evaluation_c = evaluation_a;
      }
      const double tol
            = 2.0 * DBL_EPSILON * fabs(b) + 0.5 * params->closure_root_tolerance;
      const double midpoint = 0.5 * (c - b);
      /* A requested interval tolerance below the double spacing near b
       * cannot be met; accept the bracket once it reaches that spacing. */
      if(fabs(midpoint) <= params->closure_root_tolerance
         || fabs(midpoint) <= 2.0 * DBL_EPSILON * fabs(b) || fb == 0.0) {
        candidate.xi = b;
        candidate.status = ghl_m1_closure_solve_converged;
        ++candidate.iterations;
        candidate.evaluation = evaluation_b;
        break;
      }
      if(fabs(e) >= tol && fabs(fa) > fabs(fb)) {
        const double s = fb / fa;
        double p, q;
        if(a == c) {
          p = 2.0 * midpoint * s;
          q = 1.0 - s;
        }
        else {
          const double q1 = fa / fc;
          const double r = fb / fc;
          p = s * (2.0 * midpoint * q1 * (q1 - r) - (b - a) * (r - 1.0));
          q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
        }
        if(p > 0.0) {
          q = -q;
        }
        else {
          p = -p;
        }
        const double accept_bound = 3.0 * midpoint * q - fabs(tol * q);
        const double interpolation_bound = fmin(accept_bound, fabs(e * q));
        /* q == 0 makes the bound zero; nonnegative 2*p already rejects it. */
        if(2.0 * p < interpolation_bound) {
          e = d;
          d = p / q;
        }
        else {
          d = midpoint;
          e = midpoint;
        }
      }
      else {
        d = midpoint;
        e = midpoint;
      }
      a = b;
      fa = fb;
      evaluation_a = evaluation_b;
      /* Splitting the ternary keeps the branch record honest: both arms add
       * to b in the same order, and the optimized merged ternary edges are
       * not separately executable. */
      if(fabs(d) > tol) {
        b += d;
      }
      else {
        b += copysign(tol, midpoint);
      }
      ghl_m1_closure_evaluation next_evaluation;
      error = evaluate(context, b, &next_evaluation);
      if(error != ghl_success) {
        return error;
      }
      evaluation_b = next_evaluation;
      fb = evaluation_b.residual;
    }
    if(candidate.status == ghl_m1_closure_solve_iteration_exhausted) {
      candidate.xi = b;
      candidate.evaluation = evaluation_b;
    }
  }

  *root = candidate;
  return ghl_success;
}

static ghl_error_codes_t evaluate_workspace_minerbo(
      void *const context,
      const double xi,
      ghl_m1_closure_evaluation *const evaluation) {
  return evaluate_minerbo((const minerbo_workspace *)context, xi, evaluation);
}

static ghl_error_codes_t ghl_m1_compute_closure_minerbo_internal(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure) {
  if(m1_params == NULL || metric == NULL || prims == NULL || rad_state == NULL
     || closure == NULL) {
    return ghl_error_m1_null_pointer;
  }
  ghl_error_codes_t error
        = ghl_m1_validate_realizability_state(m1_params, metric, rad_state, 128.0, NULL);
  if(error != ghl_success) {
    return error;
  }

  minerbo_workspace ws;
  error = build_minerbo_workspace(metric, prims, rad_state, &ws);
  if(error != ghl_success) {
    return error;
  }

  ghl_m1_closure_root_result root;
  error = ghl_m1_closure_private_solve_root(
        m1_params, evaluate_workspace_minerbo, &ws, &root);
  if(error != ghl_success) {
    return error;
  }

  error = publish_minerbo(
        m1_params, &ws, &root.evaluation, root.iterations, root.status,
        m1_params->closure_root_residual_tolerance, closure);
  if(error != ghl_success) {
    return error;
  }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_closure_with_primitives(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure) {
  return ghl_m1_compute_closure_minerbo_internal(
        m1_params, metric, prims, rad_state, closure);
}

ghl_error_codes_t ghl_m1_compute_closure_minerbo_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure) {
  return ghl_m1_compute_closure_minerbo_internal(
        m1_params, metric, prims, rad_state, closure);
}
