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


// Eq (21)
// f^a = u^a + H^a/J
void assemble_fnu(
        const ghl_ADM_aux_quantities *adm_aux,
        const double *u4U,
        const double J,
        ghl_radiation_flux_vector const * H4,
        ghl_radiation_con_flux_vector * fnu4) {
  double H4U[4];
  for (int a = 0; a < 4; a++){
    for (int b = 0; b < 4; b++){
      H4U[a] += adm_aux->g4UU[a][b] * H4->D[b];
    }
  }
  for (int a = 0; a < 4; a++) {
      fnu4->U[a] = u4U[a] + H4U[a]/J;
  }
}

// TODO: rad_E_floor rad_eps are param
// Eq (24) -> Gamma = W(E/J)(1-f.v)
double compute_Gamma(
        const double W,
        const double *v4U,
        const double J,
        const double E,
        const double rad_E_floor,
        const double rad_eps,
        ghl_radiation_flux_vector const * F4) {
          
    if (E > rad_E_floor && J > rad_E_floor) {
        double f_dot_v;
        double f_dot_v_sum = 0;
        for (int a = 0; a < 4; a++){
          f_dot_v_sum += F4->D[a] * v4U[a];
        }
        f_dot_v = fmin(f_dot_v_sum, 1 - rad_eps);
        return W*(E/J)*(1 - f_dot_v);
    }
    else {
        return 1;
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


//rad_E_floor rad_eps is a param
void apply_floor(
        const ghl_ADM_aux_quantities *adm_aux,
        double * E,
        ghl_radiation_flux_vector * F4,
        const double rad_E_floor,
        const double rad_eps) {
    *E = fmax(rad_E_floor, *E);

    double F2 = 0; // needs a better name
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            F2 += adm_aux->g4UU[a][b] * F4->D[a] * F4->D[b];
        }
    }
    const double lim = (*E)*(*E)*(1 - rad_eps);
    if (F2 > lim) {
        double fac = lim/F2;
        for (int a = 0; a < 4; ++a) {
            F4->D[a] *= fac;
        }
    }
}




// Eq (29) - 2
double calc_rE_source(
        const ghl_metric_quantities *metric,
        const ghl_radiation_con_source_vector *S_src) {
  double rE_source = 0.0;
  double nD4[4] = {-metric->lapse, 0.0, 0.0, 0.0};
  for (int a = 0; a < 4; ++a) {
        rE_source += - metric->lapse*nD4[a]* S_src->U[a];
  }
  return rE_source;
}


// Eq (29) - 3
void calc_rF_source(
        const ghl_metric_quantities *metric,
        const ghl_ADM_aux_quantities *adm_aux,
        const ghl_radiation_con_source_vector *S_src,
        ghl_radiation_con_source_vector * F_src) {
    for (int a = 0; a < 4; ++a) {
        F_src->U[a] = 0.0;
        for (int b = 0; b < 4; ++b) {
            F_src->U[a] += metric->lapse*adm_aux->g4DD[b][a]*S_src->U[b];
        }
    }
}

