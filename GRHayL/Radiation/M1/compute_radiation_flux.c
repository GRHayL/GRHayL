#include "ghl.h"
#include "ghl_m1.h"

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
        const double * u4D,
        const double J,
        const ghl_radiation_flux_vector *H4,
        ghl_radiation_con_source_vector *S4) {
    for (int a = 0; a < 4; ++a) {
        S4->U[a] = (eta - kabs*J) * u4D[a] - (kabs + kscat) * H4->D[a];
    }
}

// Eq (28) - 2
double calc_E_flux(
        const ghl_metric_quantities *metric,
        const double E,
        const ghl_radiation_flux_vector *F4,
        const int dir) {
    return metric->lapse*F4->D[dir] - metric->betaU[dir]*E;
}

// Eq (28) - 3
// should I just return 3*3 tensor?
double calc_F_flux(
        const ghl_metric_quantities *metric,
        const ghl_radiation_flux_vector *F4,
        const ghl_radiation_pressure_tensor *P4,
        const int dir,
        const int comp) {
        
  // need to change this to UD
  return metric->lapse*P4->DD[dir][comp] - metric->betaU[dir]* F4->D[comp];
}

// Eq (29) - 2
double calc_rE_source(
        const ghl_metric_quantities *metric,
        const ghl_radiation_con_source_vector *S4) {
  double rE_source = 0.0;
  double nD4[4] = {-metric->lapse, 0.0, 0.0, 0.0};
  for (int a = 0; a < 4; ++a) {
        rE_source += - metric->lapse*nD4[a]* S4->U[a];
  }
  return rE_source;
}


// Eq (29) - 3
void calc_rF_source(
        const ghl_metric_quantities *metric,
        const ghl_ADM_aux_quantities *adm_aux,
        const ghl_radiation_con_source_vector *S4,
        ghl_radiation_con_source_vector * F_src) {
    for (int a = 0; a < 4; ++a) {
        F_src->U[a] = 0.0;
        for (int b = 0; b < 4; ++b) {
            F_src->U[a] += metric->lapse*adm_aux->g4DD[b][a]*S4->U[b];
        }
    }
}

