#ifndef GHL_RADIATION_H_
#define GHL_RADIATION_H_

#include "ghl.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>

// Neutrino quantities
typedef struct ghl_neutrino_luminosities {
  double nue, anue, nux;
} ghl_neutrino_luminosities;

typedef struct ghl_neutrino_opacities {
  double nue[2], anue[2], nux[2];
} ghl_neutrino_opacities;

typedef ghl_neutrino_opacities ghl_neutrino_optical_depths;

typedef struct {
  double UU[4][4];
  double UD[4][4];
  double DD[4][4];
} ghl_radiation_pressure_tensor;

typedef struct {
  double UD[4][4];
} ghl_radiation_metric_tensor;

typedef struct {
  double D[4];
  double U[4];
} ghl_radiation_flux_vector;

// conservative flux for M1 equations
typedef struct {
  double U[4];
  double D[4];
} ghl_radiation_con_flux_vector;

// conservative source for M1 equaitons
typedef struct {
  double U[4];
  double D[4];
} ghl_radiation_con_source_vector;

// conservative source tensor for M1 equaitons
typedef struct {
  double DD[4][4];
} ghl_radiation_con_source_tensor;

typedef struct {
  ghl_metric_quantities *metric;
  ghl_ADM_aux_quantities *adm_aux;
  ghl_primitive_quantities *prims;
  double E;
  ghl_radiation_flux_vector *F4;
  ghl_radiation_pressure_tensor *P4;
  double chi;
} m1_root_params;

// Choice of closure scheme
// Do not move this to ghl_m1.h
typedef enum { Eddington, Kershaw, Minerbo, Thin } ghl_m1_closure_t;
double eddington(const double xi);
double kershaw(const double xi);
double minerbo(const double xi);
double thin(const double xi);
// The function pointer to be used for the closure.
extern double (*ghl_m1_closure)(double);
ghl_error_codes_t ghl_initialize_m1_closure(ghl_m1_closure_t ghl_m1_closure_type);

#include "nrpyleakage.h"

void add_rad_to_Tmunu(
      ghl_stress_energy *T4DD,
      const ghl_stress_energy *rT4DD,
      const ghl_metric_quantities *metric);

void calc_density_ratios(
      const double mb,
      const double rho,
      const double eps,
      const ghl_metric_quantities *metric,
      const double n_e,
      const double n_ae,
      const double n_x,
      const double J_e,
      const double J_ae,
      const double J_x,
      double *number_ratio_e,
      double *number_ratio_ae,
      double *number_ratio_x,
      double *energy_ratio_e,
      double *energy_ratio_ae,
      double *energy_ratio_x);

void flavor_mix_adjustment(
      double *rN_e,
      double *rN_ae,
      double *rN_x,
      double *rN_ax,
      const double Nmin,
      double *rE_e,
      double *rE_ae,
      double *rE_x,
      double *rE_ax,
      ghl_radiation_flux_vector *rF_e,
      ghl_radiation_flux_vector *rF_ae,
      ghl_radiation_flux_vector *rF_x,
      ghl_radiation_flux_vector *rF_ax,
      const bool mix_type);

void ghl_calc_neutrino_densities(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *nudens_0,
      double *nudens_1);

#endif  // GHL_RADIATION_H_