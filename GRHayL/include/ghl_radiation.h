#ifndef GHL_RADIATION_H_
#define GHL_RADIATION_H_

#include "ghl.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>

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
//   ghl_stress_energy rT4DD;

  // SE tensor in the fluid frame, will be updated in impl solve
  // double J;
  // ghl_radiation_flux_vector *H4;

  // sources
  // THC powell param initialize T_dd (rT4DD) H_d (H4) S_d (S4) tS_d (rF_source)
  // but only tS_d is needed for the powell solver
  //   ghl_radiation_con_source_vector *S4;
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

void ghl_m1_set_equilibrium(
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims,
      const ghl_metric_quantities *restrict adm_metric,
      const ghl_ADM_aux_quantities *restrict aux_metric,
      const double nudens_0,
      const double nudens_1,
      double *restrict rN,
      double *restrict rnnu,
      double *restrict rE,
      double *restrict rJ,
      ghl_radiation_flux_vector *restrict rF4,
      ghl_radiation_flux_vector *restrict rH4,
      ghl_radiation_pressure_tensor *restrict rP4);

void ghl_M1_update(
      double cdt,
      ghl_neutrino_optical_depths *restrict tau,
      ghl_neutrino_opacities *restrict kappa,
      ghl_metric_quantities *metric,
      ghl_ADM_aux_quantities *adm_aux,
      ghl_radiation_metric_tensor *proj4,
      ghl_primitive_quantities *prims);

void assemble_fnu(
      const ghl_ADM_aux_quantities *adm_aux,
      const double *u4U,
      const double J,
      const ghl_radiation_flux_vector *H4,
      ghl_radiation_con_flux_vector *fnu4);

double compute_Gamma(
      const double W,
      const double *v4U,
      const double J,
      const double E,
      const double rad_E_floor,
      const double rad_eps,
      const ghl_radiation_flux_vector *F4);

void ghl_radiation_compute_pressure_tensor_thick(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector *F4,
      ghl_radiation_pressure_tensor *P_thick);

void ghl_radiation_compute_pressure_tensor_thin(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector *F4,
      ghl_radiation_pressure_tensor *P_thin);

void ghl_radiation_apply_closure(
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E,
      const ghl_radiation_flux_vector *F4,
      const double chi,
      ghl_radiation_pressure_tensor *P4);

double eddington(const double xi);
double kershaw(const double xi);
double minerbo(const double xi);
double thin(const double xi);

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

void get_FU(const ghl_ADM_aux_quantities *adm_aux, ghl_radiation_flux_vector *F4);

void get_PUD_and_PUU(
      const ghl_ADM_aux_quantities *adm_aux,
      ghl_radiation_pressure_tensor *P4);

void calc_rad_sources(
      const double eta,
      const double kabs,
      const double kscat,
      const double *u4U,
      const double J,
      const ghl_radiation_flux_vector *H4,
      ghl_radiation_con_source_vector *S4);

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

void init_params(ghl_m1_powell_params *p);

void apply_floor(
      const ghl_ADM_aux_quantities *adm_aux,
      double *E,
      ghl_radiation_flux_vector *F4,
      const double rad_E_floor,
      const double rad_eps);

void __source_jacobian_low_level(
      double qpre[4],
      double Fup[4],
      double F2,
      double chi,
      double kapa,
      double kaps,
      double vup[4],
      double vdown[4],
      double v2,
      double W,
      double alpha,
      double cdt,
      gsl_matrix *J);

double dot(double *a, double *b, int length);

int prepare_closure(const gsl_vector *q, ghl_m1_powell_params *p);

int prepare_sources(const gsl_vector *q, ghl_m1_powell_params *p);

int prepare(const gsl_vector *q, ghl_m1_powell_params *p);

void evaluate_zjac(ghl_m1_powell_params *p, gsl_matrix *J);

void evaluate_zfunc(gsl_vector const * q, ghl_m1_powell_params *p, gsl_vector *f);

int impl_func_jac(const gsl_vector *q, void *params, gsl_matrix *J);

int impl_func_val(const gsl_vector *q, void *params, gsl_vector *f);

int impl_func_val_jac(const gsl_vector *q, void *params, gsl_vector *f, gsl_matrix *J);

// void explicit_update(
//       ghl_m1_powell_params *p,
//       double *E_new,
//       ghl_radiation_flux_vector *F4_new);
void explicit_update(
      ghl_m1_powell_params *p,
      double *E_new,
      ghl_radiation_flux_vector *F4_new);

int ghl_source_update_test(const ghl_m1_thc_params *thc_params, ghl_m1_powell_params *p);

int ghl_source_update(
      const ghl_m1_thc_params *thc_params,
      ghl_m1_powell_params *p,
      double *E_new,
      ghl_radiation_flux_vector *F4_new
);

void ghl_calc_neutrino_densities(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *nudens_0,
      double *nudens_1);

#endif // GHL_RADIATION_H_
