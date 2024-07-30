#include "ghl.h"
#include "ghl_m1.h"


// This function computes Eq. (6) of Radice et al. (2022)
void ghl_radiation_compute_pressure_tensor_fluid_frame(
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double J,
      ghl_radiation_pressure_tensor *K) {

  // u^{i} = v^{i} * u^{0}
  const double u4U[4] = { prims->u0, prims->u0 * prims->vU[0], prims->u0 * prims->vU[1],
                          prims->u0 * prims->vU[2] };

  double u4D[4] = { 0, 0, 0, 0 };
  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      u4D[mu] += adm_aux->g4DD[mu][nu] * u4U[nu];
    }
  }

  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      K->DD[mu][nu] = J * (adm_aux->g4DD[mu][nu] + u4D[mu] * u4D[nu]);
    }
  }
}

// Computes radiation pressure in the lab frame in the thick limit (Eq. (8) of Radice et al. (2022))
void ghl_radiation_compute_pressure_tensor_thick(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector* F4,
      ghl_radiation_pressure_tensor *P_thick) {

  // u^{i} = v^{i} * u^{0}
  const double u4U[4] = { prims->u0, prims->u0 * prims->vU[0], prims->u0 * prims->vU[1],
                          prims->u0 * prims->vU[2] };
  const double v4U[4] = { 0, prims->vU[0], prims->vU[1], prims->vU[2] };

  double u4D[4] = { 0, 0, 0, 0 };
  double v4D[4] = { 0, 0, 0, 0 };
  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      u4D[mu] += adm_aux->g4DD[mu][nu] * u4U[nu];
      v4D[mu] += adm_aux->g4DD[mu][nu] * v4U[nu];
    }
  }

  const double n4U[4] = {metric->lapseinv, -metric->betaU[0]*metric->lapseinv
                                   , -metric->betaU[1]*metric->lapseinv
                                   , -metric->betaU[2]*metric->lapseinv};
  double n4D[4] = {metric->lapse, 0, 0, 0};

  double W = 0;

  // fluid Lorentz factor W = -uU nD
  for(int mu = 0; mu < 4; mu++){
    W += -u4D[mu] * n4U[mu];
  }
  const double W2 = W*W;
  const double coef = 1./(2.*W2 + 1.);

  // v.F = vU FD
  double v_dot_F = 0;
  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      v_dot_F += v4U[mu] * F4->D[nu];
    }
  }

  // J/3
  const double Jo3 = coef*((2.*W2 - 1.)*E - 2.*W2*v_dot_F);

  // tH = gamma_UD H_D
  double tH4D[4];
  for (int mu = 0; mu < 4; mu++) {
      tH4D[mu] = F4->D[mu]/W + coef*W*v4D[mu]*((4.*W2 + 1.)*v_dot_F - 4.*W2*E);
  }

  for (int mu = 0; mu < 4; mu++)
  for (int nu= mu; nu < 4; nu++) {
      P_thick->DD[mu][nu]  = Jo3 * (4. * W2 * v4D[mu]*v4D[nu] + adm_aux->g4DD[mu][nu] + n4D[mu]*n4D[nu]);
      P_thick->DD[mu][nu] += W*(tH4D[mu]*v4D[nu] + tH4D[nu]*v4D[mu]);
  }
} // ghl_radiation_compute_pressure_tensor_thick


// Computes radiation pressure in the lab frame in the thin limit (Eq. (15) of Radice et al. (2022))
void ghl_radiation_compute_pressure_tensor_thin(
    const ghl_metric_quantities *metric,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_primitive_quantities *prims,
    const double E,
    const ghl_radiation_flux_vector* F4,
    ghl_radiation_pressure_tensor *P_thin) {
  double F2 = 0;
  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      F2 += F4->D[mu] * F4->D[nu] * adm_aux->g4UU[mu][nu];
    }
  }
  double fac = (F2 > 0 ? E/F2 : 0);
  for (int a = 0; a < 4; a++) {
    for (int b = a; b < 4; b++) {
      P_thin->DD[a][b] = fac * F4->D[a] * F4->D[b];
    }
  }
} // ghl_radiation_compute_pressure_tensor_thin

void ghl_radiation_apply_closure(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector* F4,
      const double chi,
      ghl_radiation_pressure_tensor *P4DD){
  ghl_radiation_pressure_tensor P_thin;
  ghl_radiation_pressure_tensor P_thick;
  ghl_radiation_compute_pressure_tensor_thin(metric, adm_aux, prims, E, F4->D, &P_thin);
  ghl_radiation_compute_pressure_tensor_thick(metric, adm_aux, prims, E, F4->D, &P_thick);

  const double dthick = 3.*(1 - chi)/2.;
  const double dthin = 1. - dthick;

  for (int mu =  0; mu < 4; mu++){
    for (int nu = mu; nu < 4; nu++) {
        P4DD->DD[mu][nu] = dthick*P_thick.DD[mu][nu] + dthin*P_thin.DD[mu][nu];
    }
  }
} //ghl_radiation_apply_closure

double eddington(const double xi) {
    return 1.0/3.0;
}

double kershaw(const double xi) {
    return 1.0/3.0 + 2.0/3.0*xi*xi;
}

double minerbo(const double xi) {
    return 1.0/3.0 + xi*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0;
}

double thin(const double xi) {
    return 1.0;
}

typedef struct {
  double J, H2;
} fparams;

// Function to rootfind in order to determine the closure
static inline double froot(
  const ghl_metric_quantities *metric,
  const ghl_ADM_aux_quantities *adm_aux,
  const ghl_primitive_quantities *prims,
  const double E,
  const ghl_radiation_flux_vector* F4,
  const ghl_radiation_pressure_tensor *P4DD,
  const double xi,
  const ghl_m1_parameters* restrict params,
  void *restrict fparams
) {
    double chi = minerbo(xi);
    ghl_radiation_apply_closure(metric, adm_aux, prims, E, F4->D, chi, P4DD);

    double rT_dd[4][4];
    const double n4D[4] = {metric->lapse, 0, 0, 0};
    assemble_rT(n4D, E, F4->D, P4DD, &rT_dd);

    // u^{i} = v^{i} * u^{0}
    const double u4U[4] = { prims->u0, prims->u0 * prims->vU[0], prims->u0 * prims->vU[1],
                          prims->u0 * prims->vU[2] };
    double proj4UD[4][4];
    for (int a = 0; a < 4; a++){
      for (int b = 0; b < 4; b++){
        for (int c = 0; c < 4; c++){
          proj4UD[a][b] += (a==b) - u4U[a] * u4U[c] * adm_aux->g4DD[b][c];
        }
      }
    }

    const double J = calc_J_from_rT(u4U, proj4UD, rT_dd);

    ghl_radiation_flux_vector H4;
    calc_H_from_rT(rT_dd, u4U, proj4UD, &H4);
    double H2 = 0;
    for (int a = 0; a < 4; a++){
      for (int b = 0; b < 4; b++){
        H2 += adm_aux->g4UU[a][b] * u4U[a] * u4U[b];
      }
    }
  return xi*xi*params->J*params->J - params->H2;
} // static inline double froot

// Computes SE tensor in the fluid frame (Eq. (2) of Radice et al. (2022))
void assemble_rT_fluid_frame(
        const double *  u4D,
        const double J,
        const ghl_radiation_flux_vector *H4,
        const ghl_radiation_pressure_tensor *K4,
        ghl_stress_energy *rT4){
  for (int mu = 0; mu < 4; mu++){
    for (int nu = mu; nu < 4; nu++) {
      rT4->T4[mu][nu] = J*u4D[mu]*u4D[nu] + H4->D[mu]*u4D[nu] + H4->D[nu]*u4D[mu] + K4->DD[mu][nu];
    }
  }
} //assemble_rT_fluid_frame

// Computes SE tensor in the lab frame (Eq. (1) of Radice et al. (2022))
void assemble_rT_lab_frame(
        const double *  n4D,
        const double E,
        const ghl_radiation_flux_vector *F4,
        const ghl_radiation_pressure_tensor *P4,
        ghl_stress_energy *rT4) {
  for (int mu = 0; mu < 4; mu++){
    for (int nu = mu; nu < 4; nu++) {
      rT4->T4[mu][nu] = E*n4D[mu]*n4D[nu] + F4->D[mu]*n4D[nu] + F4->D[nu]*n4D[mu] + P4->DD[mu][nu];
    }
  }
} //assemble_rT_lab_frame

// Computes radiation energy densities in fluid frame (Eq. (2) of Radice et al. (2022))
double calc_J_from_rT(
        const double  *u4U,
        const double **proj4UD, 
        const ghl_stress_energy *rT4) {
  double J = 0;
  for (int mu = 0; mu < 4; mu++){
    for (int nu = 0; nu < 4; nu++) {
        J += rT4->T4[mu][nu] * u4U[mu] *u4U[nu];
    }
  }
} //calc_J_from_rT

// Computes radiation flux in fluid frame (Eq. (2) of Radice et al. (2022))
void* calc_H4D_from_rT(
        const double  *u4U,
        const double **proj4UD, 
        const ghl_radiation_pressure_tensor *rT4,
        ghl_radiation_flux_vector *H4) {
  for (int a = 0; a < 4; a++){
    for (int b = 0; b < 4; b++) {
      for(int c = 0; c < 4; c++) {
        H4->D[a] += proj4UD[b][a] * u4U[c] * rT4->DD[b][c];
      }
    }
  }
} // calc_H4D_from_rT

// Computes radiation pressure tensor in fluid frame (Eq. (2) of Radice et al. (2022))
void calc_K4DD_from_rT(
        const double  *u4U,
        const double **proj4UD, 
        const ghl_radiation_pressure_tensor *rT4,
        ghl_radiation_pressure_tensor *K4) {
  for (int a = 0; a < 4; a++){
    for (int b = 0; b < 4; b++) {
      for(int c = 0; c < 4; c++) {
        for(int d = 0; d < 4; d++){
          K4->DD[a][b] += proj4UD[c][a] * proj4UD[d][b] * rT4->DD[c][d];
        }
      }
    }
  }
} // calc_K4DD_from_rT

int ghl_radiation_rootSolve_closure(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      fparams_struct *restrict fparams
      ) {
  // Step 5: Set specific quantities for this routine (Eq. A7 of [1])
  fparams_struct fparams;

  // Step 6: Bracket x (Eq. A8 of [1])
  double xlow = 0;
  double xup  = 1;

  // Step 8: Call the main function
  roots_params rparams;
  rparams.tol = 1e-15;
  rparams.max_iters = 300;


  ghl_brent(froot, params, eos, &fparams, xlow, xup, &rparams);
  return 0;
}

// Eq (23)
void assemble_fnu(
        const double *u4U,
        const double J,
        ghl_radiation_flux_vector const * H4,
        ghl_radiation_flux_vector * fnu4) {
  for (int a = 0; a < 4; a++) {
      fnu4->D[a] = u4U[a] + H4->D[a]/J;
  }
}

// TODO: rad_E_floor rad_eps are param
// Eq (24)
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
        f_dot_v = min(f_dot_v_sum, 1 - rad_eps);
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
        ghl_radiation_source_vector *S4) {
    for (int a = 0; a < 4; ++a) {
        S4->U[a] = (eta - kabs*J) * u4D[a] - (kabs + kscat) * H4->D[a];
    }
}

// Eq (29) - 2
double calc_rE_source(
        const ghl_metric_quantities *metric,
        const ghl_radiation_source_vector *S_src) {
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
        const ghl_radiation_source_vector *S_src,
        ghl_radiation_source_vector * F_src) {
    for (int a = 0; a < 4; ++a) {
        F_src->U[a] = 0.0;
        for (int b = 0; b < 4; ++b) {
            F_src->U[a] += metric->lapse*adm_aux->g4DD[b][a]*S_src->U[b];
        }
    }
}

//rad_E_floor rad_eps is a param
void apply_floor(
        const ghl_ADM_aux_quantities *adm_aux,
        double * E,
        ghl_radiation_flux_vector * F4,
        const double rad_E_floor,
        const double rad_eps) {
    *E = max(rad_E_floor, *E);

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