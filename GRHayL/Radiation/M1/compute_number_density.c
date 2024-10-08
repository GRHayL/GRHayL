#include "ghl.h"
#include "ghl_m1.h"

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
