#ifndef GHL_M1_H_
#define GHL_M1_H_

#include "ghl.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// compute_radiation_flux.c functions
double calc_GE_source(
    const ghl_metric_quantities *metric,
    const ghl_metric_quantities *metric_derivs_x,
    const ghl_metric_quantities *metric_derivs_y,
    const ghl_metric_quantities *metric_derivs_z,
    const ghl_radiation_pressure_tensor *P4,
    const ghl_radiation_flux_vector *F4,
    const ghl_extrinsic_curvature *K4);

double calc_E_flux(
    const ghl_metric_quantities *metric,
    const double E,
    const ghl_radiation_flux_vector *F4,
    const int dir);

double calc_F_flux(
    const ghl_metric_quantities *metric,
    const ghl_radiation_flux_vector *F4,
    const ghl_radiation_pressure_tensor *P4,
    const int dir,
    const int comp);

double calc_rE_source(
    const ghl_metric_quantities *metric,
    const ghl_radiation_con_source_vector *S4);

void calc_rF_source(
    const ghl_metric_quantities *metric,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_radiation_con_source_vector *S4,
    ghl_radiation_con_source_vector *rF_source);

double calc_GE_source(
    const ghl_metric_quantities *metric,
    const ghl_metric_quantities *metric_derivs_x,
    const ghl_metric_quantities *metric_derivs_y,
    const ghl_metric_quantities *metric_derivs_z,
    const ghl_radiation_pressure_tensor *P4,
    const ghl_radiation_flux_vector *F4,
    const ghl_extrinsic_curvature *K4);

void calc_GF_source(
    const ghl_metric_quantities *metric,
    const ghl_metric_quantities *metric_derivs_x,
    const ghl_metric_quantities *metric_derivs_y,
    const ghl_metric_quantities *metric_derivs_z,
    const double E,
    const ghl_radiation_flux_vector *F4,
    const ghl_radiation_pressure_tensor *P4,
    ghl_radiation_con_source_vector *GF_source);

#endif // GHL_M1_H_