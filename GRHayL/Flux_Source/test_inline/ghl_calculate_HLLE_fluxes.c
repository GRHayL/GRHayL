#include "ghl_flux_source.h"
/*
 *  * Compute the HLLE-derived fluxes on the left face for all components.
 *   */

/* NRPYSTART
import sympy as sp
import nrpy.indexedexp as ixp
import nrpy.c_codegen as ccg
import ghl_generate_fluxes as ghl

S_flux_rhs = ixp.zerorank1()
rho_flux_rhs, tau_flux_rhs, S_flux_rhs, Ye_flux_rhs, ent_flux_rhs = ghl.calculate_fluxes_rhs()


with open("fluxes_hybrid.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2]], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]"]))

with open("fluxes_hybrid_entropy.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2], ent_flux_rhs], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]", "cons_rhs->entropy"]))

with open("fluxes_tabulated.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2], Ye_flux_rhs], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]", "cons_rhs->Y_e"]))

with open("fluxes_tabulated_entropy.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2], Ye_flux_rhs, ent_flux_rhs], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]", "cons_rhs->Y_e", "cons_rhs->entropy"]))
NRPYEND */

void ghl_calculate_HLLE_fluxes_hybrid(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

const double gammaDD00 = ADM_metric->gammaDD[0][0]
const double gammaDD01 = ADM_metric->gammaDD[0][1]
const double gammaDD02 = ADM_metric->gammaDD[0][2]
const double gammaDD11 = ADM_metric->gammaDD[1][1]
const double gammaDD12 = ADM_metric->gammaDD[1][2]
const double gammaDD22 = ADM_metric->gammaDD[2][2]

const double gammaUU00 = ADM_metric->gammaUU[0][0]
const double gammaUU01 = ADM_metric->gammaUU[0][1]
const double gammaUU02 = ADM_metric->gammaUU[0][2]
const double gammaUU11 = ADM_metric->gammaUU[1][1]
const double gammaUU12 = ADM_metric->gammaUU[1][2]
const double gammaUU22 = ADM_metric->gammaUU[2][2]

const double g4DD00 = metric_aux.g4DD[0][0]
const double g4DD01 = metric_aux.g4DD[0][1]
const double g4DD02 = metric_aux.g4DD[0][2]
const double g4DD03 = metric_aux.g4DD[0][3]
const double g4DD11 = metric_aux.g4DD[1][1]
const double g4DD12 = metric_aux.g4DD[1][2]
const double g4DD13 = metric_aux.g4DD[1][3]
const double g4DD22 = metric_aux.g4DD[2][2]
const double g4DD23 = metric_aux.g4DD[3][3]
const double g4DD33 = metric_aux.g4DD[3][3]

const double g4UU00 = metric_aux.g4UU[0][0]
const double g4UU01 = metric_aux.g4UU[0][1]
const double g4UU02 = metric_aux.g4UU[0][2]
const double g4UU03 = metric_aux.g4UU[0][3]
const double g4UU11 = metric_aux.g4UU[1][1]
const double g4UU12 = metric_aux.g4UU[1][2]
const double g4UU13 = metric_aux.g4UU[1][3]
const double g4UU22 = metric_aux.g4UU[2][2]
const double g4UU23 = metric_aux.g4UU[3][3]
const double g4UU33 = metric_aux.g4UU[3][3]
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

  #include "fluxes_hybrid.h"
}

void ghl_calculate_HLLE_fluxes_hybrid_entropy(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

const double gammaDD00 = ADM_metric->gammaDD[0][0]
const double gammaDD01 = ADM_metric->gammaDD[0][1]
const double gammaDD02 = ADM_metric->gammaDD[0][2]
const double gammaDD11 = ADM_metric->gammaDD[1][1]
const double gammaDD12 = ADM_metric->gammaDD[1][2]
const double gammaDD22 = ADM_metric->gammaDD[2][2]

const double gammaUU00 = ADM_metric->gammaUU[0][0]
const double gammaUU01 = ADM_metric->gammaUU[0][1]
const double gammaUU02 = ADM_metric->gammaUU[0][2]
const double gammaUU11 = ADM_metric->gammaUU[1][1]
const double gammaUU12 = ADM_metric->gammaUU[1][2]
const double gammaUU22 = ADM_metric->gammaUU[2][2]

const double g4DD00 = metric_aux.g4DD[0][0]
const double g4DD01 = metric_aux.g4DD[0][1]
const double g4DD02 = metric_aux.g4DD[0][2]
const double g4DD03 = metric_aux.g4DD[0][3]
const double g4DD11 = metric_aux.g4DD[1][1]
const double g4DD12 = metric_aux.g4DD[1][2]
const double g4DD13 = metric_aux.g4DD[1][3]
const double g4DD22 = metric_aux.g4DD[2][2]
const double g4DD23 = metric_aux.g4DD[3][3]
const double g4DD33 = metric_aux.g4DD[3][3]

const double g4UU00 = metric_aux.g4UU[0][0]
const double g4UU01 = metric_aux.g4UU[0][1]
const double g4UU02 = metric_aux.g4UU[0][2]
const double g4UU03 = metric_aux.g4UU[0][3]
const double g4UU11 = metric_aux.g4UU[1][1]
const double g4UU12 = metric_aux.g4UU[1][2]
const double g4UU13 = metric_aux.g4UU[1][3]
const double g4UU22 = metric_aux.g4UU[2][2]
const double g4UU23 = metric_aux.g4UU[3][3]
const double g4UU33 = metric_aux.g4UU[3][3]
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

  #include "fluxes_hybrid_entropy.h"
}

void ghl_calculate_HLLE_fluxes_tabulated(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

const double gammaDD00 = ADM_metric->gammaDD[0][0]
const double gammaDD01 = ADM_metric->gammaDD[0][1]
const double gammaDD02 = ADM_metric->gammaDD[0][2]
const double gammaDD11 = ADM_metric->gammaDD[1][1]
const double gammaDD12 = ADM_metric->gammaDD[1][2]
const double gammaDD22 = ADM_metric->gammaDD[2][2]

const double gammaUU00 = ADM_metric->gammaUU[0][0]
const double gammaUU01 = ADM_metric->gammaUU[0][1]
const double gammaUU02 = ADM_metric->gammaUU[0][2]
const double gammaUU11 = ADM_metric->gammaUU[1][1]
const double gammaUU12 = ADM_metric->gammaUU[1][2]
const double gammaUU22 = ADM_metric->gammaUU[2][2]

const double g4DD00 = metric_aux.g4DD[0][0]
const double g4DD01 = metric_aux.g4DD[0][1]
const double g4DD02 = metric_aux.g4DD[0][2]
const double g4DD03 = metric_aux.g4DD[0][3]
const double g4DD11 = metric_aux.g4DD[1][1]
const double g4DD12 = metric_aux.g4DD[1][2]
const double g4DD13 = metric_aux.g4DD[1][3]
const double g4DD22 = metric_aux.g4DD[2][2]
const double g4DD23 = metric_aux.g4DD[3][3]
const double g4DD33 = metric_aux.g4DD[3][3]

const double g4UU00 = metric_aux.g4UU[0][0]
const double g4UU01 = metric_aux.g4UU[0][1]
const double g4UU02 = metric_aux.g4UU[0][2]
const double g4UU03 = metric_aux.g4UU[0][3]
const double g4UU11 = metric_aux.g4UU[1][1]
const double g4UU12 = metric_aux.g4UU[1][2]
const double g4UU13 = metric_aux.g4UU[1][3]
const double g4UU22 = metric_aux.g4UU[2][2]
const double g4UU23 = metric_aux.g4UU[3][3]
const double g4UU33 = metric_aux.g4UU[3][3]
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

  #include "fluxes_tabulated.h"
}

void ghl_calculate_HLLE_fluxes_tabulated_entropy(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

const double gammaDD00 = ADM_metric->gammaDD[0][0]
const double gammaDD01 = ADM_metric->gammaDD[0][1]
const double gammaDD02 = ADM_metric->gammaDD[0][2]
const double gammaDD11 = ADM_metric->gammaDD[1][1]
const double gammaDD12 = ADM_metric->gammaDD[1][2]
const double gammaDD22 = ADM_metric->gammaDD[2][2]

const double gammaUU00 = ADM_metric->gammaUU[0][0]
const double gammaUU01 = ADM_metric->gammaUU[0][1]
const double gammaUU02 = ADM_metric->gammaUU[0][2]
const double gammaUU11 = ADM_metric->gammaUU[1][1]
const double gammaUU12 = ADM_metric->gammaUU[1][2]
const double gammaUU22 = ADM_metric->gammaUU[2][2]

const double g4DD00 = metric_aux.g4DD[0][0]
const double g4DD01 = metric_aux.g4DD[0][1]
const double g4DD02 = metric_aux.g4DD[0][2]
const double g4DD03 = metric_aux.g4DD[0][3]
const double g4DD11 = metric_aux.g4DD[1][1]
const double g4DD12 = metric_aux.g4DD[1][2]
const double g4DD13 = metric_aux.g4DD[1][3]
const double g4DD22 = metric_aux.g4DD[2][2]
const double g4DD23 = metric_aux.g4DD[3][3]
const double g4DD33 = metric_aux.g4DD[3][3]

const double g4UU00 = metric_aux.g4UU[0][0]
const double g4UU01 = metric_aux.g4UU[0][1]
const double g4UU02 = metric_aux.g4UU[0][2]
const double g4UU03 = metric_aux.g4UU[0][3]
const double g4UU11 = metric_aux.g4UU[1][1]
const double g4UU12 = metric_aux.g4UU[1][2]
const double g4UU13 = metric_aux.g4UU[1][3]
const double g4UU22 = metric_aux.g4UU[2][2]
const double g4UU23 = metric_aux.g4UU[3][3]
const double g4UU33 = metric_aux.g4UU[3][3]
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

  #include "fluxes_tabulated_entropy.h"
}
