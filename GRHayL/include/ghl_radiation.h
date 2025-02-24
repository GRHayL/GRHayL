#ifndef GHL_RADIATION_H_
#define GHL_RADIATION_H_

#include "ghl.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
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
} ghl_radiation_con_flux_vector;

// conservative source for M1 equaitons
typedef struct {
  double U[4];
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
} m1_root_params;

/*
 * Struct        : ghl_m1_parameters
 * Description   : stores M1 transport GRHayL parameters
 * Documentation :
 */
typedef struct ghl_m1_parameters {
  // ghl_con2prim_method_t main_routine, backup_routine[3];
  double J;
  double H2;

} ghl_m1_parameters;

/*
 * Struct        : ghl_m1_thc_params
 * Description   : THC M1 transport params.ccl parameters
 * Documentation :
 */
typedef struct ghl_m1_thc_params {
  double source_thick_limit;
  double source_scat_limit;
  int source_maxiter;
  double source_epsabs;
  double source_epsrel;
} ghl_m1_thc_params;


// Choice of closure scheme
typedef enum { Eddington, Kershaw, Minerbo, Thin } ghl_m1_closure_t;

#include "nrpyleakage.h"

// Closure function prototypes
double eddington(const double xi);
double kershaw(const double xi);
double minerbo(const double xi);
double thin(const double xi);

// The function pointer to be used for the closure.
extern double (*ghl_m1_closure)(double);

ghl_error_codes_t ghl_initialize_m1_closure(ghl_m1_closure_t ghl_m1_closure_type);

/*
 * Struct        : ghl_m1_powell_params
 * Description   : stores M1 transport GRHayL parameters for Powell's method
 * Documentation :
 */
typedef struct ghl_m1_powell_params {
  // Parameters initialized from init_params
  ghl_m1_closure_t closure;
  gsl_root_fsolver *gsl_solver_1d;
  double cdt;
  ghl_metric_quantities *metric;
  ghl_ADM_aux_quantities *adm_aux;
  ghl_radiation_metric_tensor *proj4;
  ghl_primitive_quantities *prims;
  // Initial E,F value to take in init_params
  double E_star;
  ghl_radiation_flux_vector *F4_star;
  // opacity parameters
  double chi;
  double eta;
  double kabs;
  double kscat;


  // Paraters derived in init_params
  double *u4U;
  double *n4D;
  double W;

  // SE tensor in the lab frame, will be updated in impl solve
  ghl_stress_energy *rT4DD;
  double E;
  ghl_radiation_flux_vector *F4;
  ghl_radiation_pressure_tensor *P4;

  // SE tensor in the fluid frame, will be updated in impl solve
  double J;
  ghl_radiation_flux_vector *H4;
  
  // sources
  ghl_radiation_con_source_vector *S4;
  double rE_source;
  ghl_radiation_con_flux_vector *rF_source;
  double GE_source;
  ghl_radiation_con_flux_vector *F_src;
} ghl_m1_powell_params;

#endif // GHL_RADIATION_H_
