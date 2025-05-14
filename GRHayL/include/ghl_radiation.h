#ifndef GHL_RADIATION_H_
#define GHL_RADIATION_H_

#include "ghl.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// Choice of closure scheme
typedef enum { Eddington, Kershaw, Minerbo, Thin } ghl_m1_closure_t;

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

/*
 * Struct        : ghl_m1_powell_params
 * Description   : stores M1 transport GRHayL parameters for Powell's method
 * Documentation :
 */
typedef struct ghl_m1_powell_params {
  // Parameters initialized from init_params
  ghl_m1_closure_t closure;
  gsl_root_fsolver *gsl_solver_1d;
  gsl_multiroot_fdfsolver *gsl_solver_nd;
  double cdt;
  ghl_metric_quantities *metric;
  ghl_ADM_aux_quantities *adm_aux;
  ghl_primitive_quantities *prims;
  // Initial E,F value to take in init_params
  double E_star;
  ghl_radiation_flux_vector *F4_star;
  double chi;
  // opacity parameters
  double eta;
  double kabs;
  double kscat;


  // Paraters derived in init_params
  double u4U[4];
  double n4D[4];
  double W;
  double vU[4];
  double vD[4];

  // SE tensor in the lab frame, will be updated in impl solve
  // ghl_stress_energy *rT4DD;
  double E;
  ghl_radiation_flux_vector *F4;
  ghl_radiation_pressure_tensor *P4;

  // SE tensor in the fluid frame, will be updated in impl solve
  // double J;
  // ghl_radiation_flux_vector *H4;
  
  // sources
  // ghl_radiation_con_source_vector *S4;
  double rE_source;
  ghl_radiation_con_flux_vector *rF_source;
  double GE_source;
  ghl_radiation_con_flux_vector *GF_source;

  double E_new;
  ghl_radiation_flux_vector *F4_new;

} ghl_m1_powell_params;

#include "nrpyleakage.h"

// Closure function prototypes
double eddington(const double xi);
double kershaw(const double xi);
double minerbo(const double xi);
double thin(const double xi);

// The function pointer to be used for the closure.
extern double (*ghl_m1_closure)(double);

ghl_error_codes_t ghl_initialize_m1_closure(ghl_m1_closure_t ghl_m1_closure_type);

void add_rad_to_Tmunu(ghl_stress_energy *T4DD, 
                      const ghl_stress_energy *rT4DD, 
                      const ghl_metric_quantities *metric);

void calc_density_ratios(const double mb,
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

void flavor_mix_adjustment(double *rN_e,
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

// TODO(LRW): I think this is missing the calculation of rnnu. Please check.
void ghl_m1_set_equilibrium(
    const double rho,
    const double T,
    const double eta_nue,
    const double eta_nuae,
    const double eta_nux,
    const double w_lorentz,
    const double *u4D,
    const double *n4U,
    const ghl_radiation_metric_tensor *proj4,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_metric_quantities *metric,
    double *rN,
    double *rnnu,
    double *rE);

#endif // GHL_RADIATION_H_
