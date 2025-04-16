#ifndef GHL_M1_H_
#define GHL_M1_H_

#include "ghl.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// compute_radiation_closure.c functions
void assemble_rT_lab_frame(
    const double *n4D,
    const double E,
    const ghl_radiation_flux_vector *F4,
    const ghl_radiation_pressure_tensor *P4,
    ghl_stress_energy *rT4DD);

void assemble_rT_fluid_frame(
    const double *u4D,
    const double J,
    const ghl_radiation_flux_vector *H4,
    const ghl_radiation_pressure_tensor *K4,
    ghl_stress_energy *rT4DD);

double calc_J_from_rT(
    const double *u4U,
    const ghl_radiation_metric_tensor *proj4,
    const ghl_stress_energy *rT4DD);

void calc_H4D_from_rT(
    const double *u4U,
    const ghl_radiation_metric_tensor *proj4,
    const ghl_stress_energy *rT4DD,
    ghl_radiation_flux_vector *H4);

void calc_K4DD_from_rT(
    const double *u4U,
    const ghl_radiation_metric_tensor *proj4,
    const ghl_stress_energy *rT4DD,
    ghl_radiation_pressure_tensor *K4);

    int ghl_radiation_rootSolve_closure(m1_root_params *restrict fparams_in);


// compute_radiation_flux.c functions

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
    const ghl_extrinsic_curvature *curv);

void calc_GF_source(
    const ghl_metric_quantities *metric,
    const ghl_metric_quantities *metric_derivs_x,
    const ghl_metric_quantities *metric_derivs_y,
    const ghl_metric_quantities *metric_derivs_z,
    const double E,
    const ghl_radiation_flux_vector *F4,
    const ghl_radiation_pressure_tensor *P4,
    ghl_radiation_con_source_vector *GF_source);

void calc_rad_sources(
    const double eta,
    const double kabs,
    const double kscat,
    const double *u4U,
    const double J,
    const ghl_radiation_flux_vector *H4,
    ghl_radiation_con_source_vector *S4);

// compute_radiation_source.c functions
void init_params(ghl_m1_powell_params *p);

int ghl_source_update_test(const ghl_m1_thc_params *thc_params,
    ghl_m1_powell_params *p);

int prepare(gsl_vector const * q, ghl_m1_powell_params * p);


#endif // GHL_M1_H_