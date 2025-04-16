#include "ghl.h"
#include "ghl_radiation.h"
#include "ghl_roots.h"

// Computes radiation pressure in the lab frame in the thick
// limit (Eq. (8) of Radice et al. (2022))
// *** "fluid three velocity" in the paper != fluid three velocity in GRHayL
void ghl_radiation_compute_pressure_tensor_thick(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector *F4,
      ghl_radiation_pressure_tensor *P_thick) {
  // u^{i} = v^{i} * u^{0}
  const double u4U[4] = { prims->u0, prims->u0 * prims->vU[0], prims->u0 * prims->vU[1],
                          prims->u0 * prims->vU[2] };
  const double v4U[4] = { 0, prims->vU[0], prims->vU[1], prims->vU[2] };

  double u4D[4] = { 0, 0, 0, 0 };
  double v4D[4] = { 0, 0, 0, 0 };
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      u4D[a] += adm_aux->g4DD[a][b] * u4U[b];
      v4D[a] += adm_aux->g4DD[a][b] * v4U[b];
    }
  }

  const double n4U[4]
        = { metric->lapseinv, -metric->betaU[0] * metric->lapseinv,
            -metric->betaU[1] * metric->lapseinv, -metric->betaU[2] * metric->lapseinv };
  double n4D[4] = { -metric->lapse, 0, 0, 0 };

  // v.F = vU FD
  double v_dot_F = 0;
  for(int a = 0; a < 4; a++) {
    v_dot_F += v4U[a] * F4->D[a];
  }

  // fluid Lorentz factor W = -uU nD
  double W = 0;
  for(int a = 0; a < 4; a++) {
    W += -u4D[a] * n4U[a];
  }

  // J/3 - Eq (12)
  const double W2 = W * W;
  const double coef = 1. / (2. * W2 + 1.);
  const double Jo3 = coef * ((2. * W2 - 1.) * E - 2. * W2 * v_dot_F);

  // tH = gamma_UD H_D
  // Eq (14)
  double tH4D[4];
  for(int a = 0; a < 4; a++) {
    tH4D[a]
          = F4->D[a] / W + coef * W * v4D[a] * ((4. * W2 + 1.) * v_dot_F - 4. * W2 * E);
  }
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      P_thick->DD[a][b]
            = Jo3 * (4. * W2 * v4D[a] * v4D[b] + adm_aux->g4DD[a][b] + n4D[a] * n4D[b]);
      P_thick->DD[a][b] += W * (tH4D[a] * v4D[b] + tH4D[b] * v4D[a]);
    }
  }
} // ghl_radiation_compute_pressure_tensor_thick

// Computes radiation pressure in the lab frame in the thin
// limit (Eq. (15) of Radice et al. (2022))
void ghl_radiation_compute_pressure_tensor_thin(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector *F4,
      ghl_radiation_pressure_tensor *P_thin) {
  double F2 = 0;
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      F2 += F4->D[a] * F4->D[b] * adm_aux->g4UU[a][b];
    }
  }
  double fac = (F2 > 0 ? E / F2 : 0);
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      P_thin->DD[a][b] = fac * F4->D[a] * F4->D[b];
    }
  }
} // ghl_radiation_compute_pressure_tensor_thin

void ghl_radiation_apply_closure(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector *F4,
      const double chi,
      ghl_radiation_pressure_tensor *P4) {
  ghl_radiation_pressure_tensor P_thin;
  ghl_radiation_pressure_tensor P_thick;
  ghl_radiation_compute_pressure_tensor_thin(metric, adm_aux, prims, E, F4, &P_thin);
  ghl_radiation_compute_pressure_tensor_thick(metric, adm_aux, prims, E, F4, &P_thick);

  const double dthick = 3. * (1 - chi) / 2.;
  const double dthin = 1. - dthick;

  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      P4->DD[a][b] = dthick * P_thick.DD[a][b] + dthin * P_thin.DD[a][b];
    }
  }
  // Leo: we need this -> ghl_radation_compute_P4UU_and_P4UD(adm_aux, P4);
} // ghl_radiation_apply_closure

double eddington(const double xi) { return 1.0 / 3.0; }

double kershaw(const double xi) { return 1.0 / 3.0 + 2.0 / 3.0 * xi * xi; }

double minerbo(const double xi) {
  return 1.0 / 3.0 + xi * xi * (6.0 - 2.0 * xi + 6.0 * xi * xi) / 15.0;
}

double thin(const double xi) { return 1.0; }

// Computes SE tensor in the lab frame (Eq. (1) of Radice et al. (2022))
void assemble_rT_lab_frame(
      const double *n4D,
      const double E,
      const ghl_radiation_flux_vector *F4,
      const ghl_radiation_pressure_tensor *P4,
      ghl_stress_energy *rT4DD) {
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      rT4DD->T4[a][b]
            = E * n4D[a] * n4D[b] + F4->D[a] * n4D[b] + F4->D[b] * n4D[a] + P4->DD[a][b];
    }
  }
} // assemble_rT_lab_frame

// Computes SE tensor in the fluid frame (Eq. (2) of Radice et al. (2022))
// This is not used
void assemble_rT_fluid_frame(
      const double *u4D,
      const double J,
      const ghl_radiation_flux_vector *H4,
      const ghl_radiation_pressure_tensor *K4,
      ghl_stress_energy *rT4DD) {
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      rT4DD->T4[a][b]
            = J * u4D[a] * u4D[b] + H4->D[a] * u4D[b] + H4->D[b] * u4D[a] + K4->DD[a][b];
    }
  }
} // assemble_rT_fluid_frame

// Computes radiation energy densities in fluid frame (Eq. (2) of Radice et al. (2022))
double calc_J_from_rT(
      const double *u4U,
      const ghl_radiation_metric_tensor *proj4,
      const ghl_stress_energy *rT4DD) {
  double J = 0;
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      J += rT4DD->T4[a][b] * u4U[a] * u4U[b];
    }
  }
  return J;
} // calc_J_from_rT

// Computes radiation flux in fluid frame (Eq. (2) of Radice et al. (2022))
void calc_H4D_from_rT(
      const double *u4U,
      const ghl_radiation_metric_tensor *proj4,
      const ghl_stress_energy *rT4DD,
      ghl_radiation_flux_vector *H4) {
  for(int a = 0; a < 4; a++) {
    H4->D[a] = 0.0;
    for(int b = 0; b < 4; b++) {
      for(int c = 0; c < 4; c++) {
        H4->D[a] += -proj4->UD[b][a] * u4U[c] * rT4DD->T4[b][c];
      }
    }
  }
} // calc_H4D_from_rT

// Computes radiation pressure tensor in fluid frame (Eq. (2) of Radice et al. (2022))
void calc_K4DD_from_rT(
      const double *u4U,
      const ghl_radiation_metric_tensor *proj4,
      const ghl_stress_energy *rT4DD,
      ghl_radiation_pressure_tensor *K4) {
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      K4->DD[a][b] = 0.0;
      for(int c = 0; c < 4; c++) {
        for(int d = 0; d < 4; d++) {
          K4->DD[a][b] += proj4->UD[c][a] * proj4->UD[d][b] * rT4DD->T4[c][d];
        }
      }
    }
  }
} // calc_K4DD_from_rT

// Function to rootfind in order to determine the closure
//
static inline double froot(const double xi, void *restrict fparams_in) {
  m1_root_params *fparams = (m1_root_params *)fparams_in;
  const ghl_metric_quantities *metric = fparams->metric; // maybe try this
  const ghl_ADM_aux_quantities *adm_aux = fparams->adm_aux;
  const ghl_primitive_quantities *prims = fparams->prims;
  double E = fparams->E;
  ghl_radiation_flux_vector *F4 = fparams->F4;
  ghl_radiation_pressure_tensor *P4 = fparams->P4;

  // Apply closure (TODO: change this so other functions can be used)
  double chi = ghl_m1_closure(xi);
  ghl_radiation_apply_closure(metric, adm_aux, prims, E, F4, chi, P4);

  ghl_stress_energy rT4DD;
  const double n4D[4] = { -metric->lapse, 0, 0, 0 };
  assemble_rT_lab_frame(n4D, E, F4, P4, &rT4DD);

  // u^{i} = v^{i} * u^{0}
  const double u4U[4] = { prims->u0, prims->u0 * prims->vU[0], prims->u0 * prims->vU[1],
                          prims->u0 * prims->vU[2] };

  double u4D[4] = { 0, 0, 0, 0 };
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      u4D[a] += u4U[b] * adm_aux->g4DD[a][b];
    }
  }
  ghl_radiation_metric_tensor proj4;

  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      proj4.UD[a][b] = (int)(a == b) + u4U[a] * u4D[b];
    }
  }
  const double J = calc_J_from_rT(u4U, &proj4, &rT4DD);

  ghl_radiation_flux_vector H4;
  calc_H4D_from_rT(u4U, &proj4, &rT4DD, &H4);

  double H2 = 0;
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      H2 += adm_aux->g4UU[a][b] * H4.D[a] * H4.D[b];
    }
  }
  double resid = xi * xi * J * J - H2;

  printf("xi = %f, residual = %f\n", xi, resid);
  return resid;
} // static inline double froot

// main root solve closure function
int ghl_radiation_rootSolve_closure(m1_root_params *restrict fparams_in) {
  // Step 1: Set specific quantities for this routine
  // m1_root_params * fparams = (m1_root_params*) fparams_in;

  // Step 2: set min and max for xi
  double xi_low = 0;
  double xi_up = 1;

  // Step 3: Set root finding parameters
  roots_params rparams;
  rparams.tol = 1e-15;
  rparams.max_iters = 300;

  // Step 3: Run Brent root solver and print the root and residual
  ghl_brent(froot, fparams_in, xi_low, xi_up, &rparams);
  double xi_root = rparams.root;
  printf("root = %f\n", xi_root);
  printf("residual = %f\n", rparams.residual);
  fparams_in->chi = ghl_m1_closure(xi_root);

  // TODO: add error message when ghl_brent fails
  return 1;
}
