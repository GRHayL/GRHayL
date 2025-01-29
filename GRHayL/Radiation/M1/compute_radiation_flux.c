#include "ghl.h"
#include "ghl_radiation.h"

// Eq (5)
/*
 * @param eta   - emissivity
 * @param kabs  - absorption coef
 * @param kscat - scattering coef
 */
void calc_rad_sources(
      const double eta,
      const double kabs,
      const double kscat,
      const double *u4D,
      const double J,
      const ghl_radiation_flux_vector *H4,
      ghl_radiation_con_source_vector *S4) {
  for(int a = 0; a < 4; ++a) {
    S4->U[a] = (eta - kabs * J) * u4D[a] - (kabs + kscat) * H4->D[a];
  }
}

// Eq (28) - 2
double calc_E_flux(
      const ghl_metric_quantities *metric,
      const double E,
      const ghl_radiation_flux_vector *F4,
      const int dir) {
  return metric->lapse * F4->D[dir] - metric->betaU[dir] * E;
}

// Eq (28) - 3
// should I just return 3*3 tensor?
double calc_F_flux(
      const ghl_metric_quantities *metric,
      const ghl_radiation_flux_vector *F4,
      const ghl_radiation_pressure_tensor *P4,
      const int dir,
      const int comp) {

  // Leo: Compute UD from DD if this is the only place we need P^{i}_{j},
  //      otherwise add UD component to struct.

  // FIXME:need to change this to UD
  return metric->lapse * P4->DD[dir][comp] - metric->betaU[dir] * F4->D[comp];
}

// Eq (29) - 2
double calc_rE_source(
      const ghl_metric_quantities *metric,
      const ghl_radiation_con_source_vector *S4) {
  double rE_source = 0.0;
  double nD4[4] = { -metric->lapse, 0.0, 0.0, 0.0 };
  for(int a = 0; a < 4; ++a) {
    rE_source += -metric->lapse * nD4[a] * S4->U[a];
  }
  return rE_source;
}

// Eq (29) - 3
void calc_rF_source(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_radiation_con_source_vector *S4,
      ghl_radiation_con_source_vector *F_src) {
  for(int a = 0; a < 4; ++a) {
    F_src->U[a] = 0.0;
    for(int b = 0; b < 4; ++b) {
      F_src->U[a] += metric->lapse * adm_aux->g4DD[b][a] * S4->U[b];
    }
  }
}

// Eq (30) - 2
double calc_GE_source(
      const ghl_metric_quantities *metric,
      const ghl_metric_quantities *metric_derivs_x,
      const ghl_metric_quantities *metric_derivs_y,
      const ghl_metric_quantities *metric_derivs_z,
      const ghl_radiation_pressure_tensor *P4,
      const ghl_radiation_flux_vector *F4,
      const ghl_extrinsic_curvature *K4) {

  double alpha_dD[3]
        = { metric_derivs_x->lapse, metric_derivs_y->lapse, metric_derivs_z->lapse };

  // FIXME: these are trying to access struct fields that don't exist.
  double GE_source = 0.0;
  for(int i = 0; i < 3; i++) {
    GE_source += -F4->U[i+1] * alpha_dD[i]; // dlog(alpha) = (1/alpha)*dalpha, the 1/alpha
                                          // cancels out with the alpha factor in front.
    for(int j = 0; j < 3; j++) {
      GE_source += metric->lapse * (P4->UU[i+1][j+1] * K4->K[i+1][j+1]);
    }
  }

  return GE_source;
}

// Eq (30) - 3
void calc_GF_source(
      const ghl_metric_quantities *metric,
      const ghl_metric_quantities *metric_derivs_x,
      const ghl_metric_quantities *metric_derivs_y,
      const ghl_metric_quantities *metric_derivs_z,
      const double E,
      const ghl_radiation_flux_vector *F4,
      const ghl_radiation_pressure_tensor *P4,
      ghl_radiation_con_source_vector *F_src) {

  double alpha_dD[3] = { metric_derivs_x->lapse, metric_derivs_y->lapse,
                         metric_derivs_z->lapse }; // index refers to derivative index.

  double betaU_dD[3][3] = {
    { metric_derivs_x->betaU[0], metric_derivs_x->betaU[1], metric_derivs_x->betaU[2] },
    { metric_derivs_y->betaU[0], metric_derivs_y->betaU[1], metric_derivs_y->betaU[2] },
    { metric_derivs_z->betaU[0], metric_derivs_z->betaU[1], metric_derivs_z->betaU[2] }
  }; // dBeta[deriv_index][Beta_index]

  double gammaDD_dD[3][3][3] = {
    { { metric_derivs_x->gammaDD[0][0], metric_derivs_x->gammaDD[0][1],
        metric_derivs_x->gammaDD[0][2] },
      { metric_derivs_x->gammaDD[1][0], metric_derivs_x->gammaDD[1][1],
        metric_derivs_x->gammaDD[1][2] },
      { metric_derivs_x->gammaDD[2][0], metric_derivs_x->gammaDD[2][1],
        metric_derivs_x->gammaDD[2][2] } },
    { { metric_derivs_y->gammaDD[0][0], metric_derivs_y->gammaDD[0][1],
        metric_derivs_y->gammaDD[0][2] },
      { metric_derivs_y->gammaDD[1][0], metric_derivs_y->gammaDD[1][1],
        metric_derivs_y->gammaDD[1][2] },
      { metric_derivs_y->gammaDD[2][0], metric_derivs_y->gammaDD[2][1],
        metric_derivs_y->gammaDD[2][2] } },
    { { metric_derivs_z->gammaDD[0][0], metric_derivs_z->gammaDD[0][1],
        metric_derivs_z->gammaDD[0][2] },
      { metric_derivs_z->gammaDD[1][0], metric_derivs_z->gammaDD[1][1],
        metric_derivs_z->gammaDD[1][2] },
      { metric_derivs_z->gammaDD[2][0], metric_derivs_z->gammaDD[2][1],
        metric_derivs_z->gammaDD[2][2] } }
  }; // dGamma[deriv_index][Gamma_index1][Gamma_index2]

  // NOTE: F_src is unintialized because I assume that the S_source calcs are done before, and the G_source terms add on after. (See Eq (26))
  // To use this function properly, it is likely we will have to call this function (calc_GF source()) after (calc_rF_source())
  for(int i = 0; i < 3; i++) {
    F_src->U[i] += -E * alpha_dD[i];
    for(int j = 0; j < 3; j++) {
      F_src->U[i] += F4->D[j+1] * betaU_dD[i][j];
      for(int k = 0; k < 3; k++) {
        F_src->U[i] += (metric->lapse / 2) * (P4->UU[j+1][k+1] * gammaDD_dD[i][j][k]);
      }
    }
  }
}
