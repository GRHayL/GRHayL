#ifndef GHL_M1_H_
#define GHL_M1_H_

#include "ghl.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#define NGROUPS 1
#define NSPECIES 1
#define NCONS_VAR_SIZE 5 // rN, rFx, rFy, rFz, rE

/*
 * Struct        : ghl_m1_conservative_quantities
 * Description   : stores M1 conservative quantities
 * Documentation :
 */
typedef struct ghl_m1_conservative_quantities {
  double rE[NGROUPS*NSPECIES];
  ghl_radiation_flux_vector rF4[NGROUPS*NSPECIES];
  double rN[NGROUPS*NSPECIES];
} ghl_m1_conservative_quantities;

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
  double rad_E_floor;
  double rad_eps;
  double minmod_theta;
  double mindiss;
} ghl_m1_thc_params;


/*
 * Function      : compute_radiation_closure.c functions
 * Description   : Computes the radiation closure
 * Documentation :
 */
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


    /*
    * Struct        : ghl_m1_powell_params
    * Description   : stores M1 transport GRHayL parameters for Powell's method
    * Documentation :
    */
typedef struct ghl_m1_powell_params {
  ghl_m1_closure_t closure;// Parameters initialized from init_params
  gsl_root_fsolver *gsl_solver_1d;
  gsl_multiroot_fdfsolver *gsl_solver_nd;
  double cdt;
  ghl_metric_quantities *metric;
  ghl_ADM_aux_quantities *adm_aux;
  ghl_primitive_quantities *prims;
  double E_star;// Initial E,F value to take in init_params
  ghl_radiation_flux_vector *F4_star;
  double chi;
  double eta; // opacity parameters
  double kabs;
  double kscat;
  double u4U[4]; // Paraters derived in init_params
  double n4D[4];
  double W;
  double vU[4];
  double vD[4];
  double E;
  ghl_radiation_flux_vector *F4;
  ghl_radiation_pressure_tensor *P4;
  double rE_source;
  ghl_radiation_con_flux_vector *rF_source;
  double E_new;
  ghl_radiation_flux_vector *F4_new;
} ghl_m1_powell_params;


/*
 * Function      : compute_number_density.c functions
 * Description   : Computes the number density
 * Documentation :
 */
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



/*
 * Function      : compute_radiation_flux.c functions
 * Description   : Computes the radiation flux
 * Documentation :
 */
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


/*
 * Function      : compute_radiation_source.c functions
 * Description   : Computes the radiation source
 * Documentation :
 */
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

int prepare(gsl_vector const * q, ghl_m1_powell_params * p);

void explicit_update(
      ghl_m1_powell_params *p,
      double *E_new,
      ghl_radiation_flux_vector *F4_new);

int ghl_source_update(
      const ghl_m1_thc_params *thc_params,
      // ghl_m1_powell_params *p,
      double chi,
      double eta,
      double kabs,
      double kscat,
      double cdt,
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_primitive_quantities *prims,
      const double E_star,
      const ghl_radiation_flux_vector *F4_star,
      const ghl_m1_closure_t closure,
      const gsl_multiroot_fdfsolver *gsl_solver_1d,
      const gsl_multiroot_fdfsolver *gsl_solver_nd,
      double *E_new,
      ghl_radiation_flux_vector *F4_new);

int prepare_closure(const gsl_vector *q, ghl_m1_powell_params *p);

int prepare_sources(const gsl_vector *q, ghl_m1_powell_params *p);

int prepare(const gsl_vector *q, ghl_m1_powell_params *p);

void evaluate_zjac(ghl_m1_powell_params *p, gsl_matrix *J);

void evaluate_zfunc(gsl_vector const * q, ghl_m1_powell_params *p, gsl_vector *f);

int impl_func_jac(const gsl_vector *q, void *params, gsl_matrix *J);

int impl_func_val(const gsl_vector *q, void *params, gsl_vector *f);

int impl_func_val_jac(const gsl_vector *q, void *params, gsl_vector *f, gsl_matrix *J);


/*
 * Function      : calc_fluxes.c functions
 * Description   : Calculates the fluxes
 * Documentation :
 */

#define NSTENCIL_CONS 5
#define NSTENCIL_FLUX 3
#define NSTENCIL_CMAX 3

#define STENCIL_CONS_SIZE NGROUPS * NSPECIES * NSTENCIL_CONS
#define CONS_SIZE NGROUPS * NSPECIES * NSTENCIL_CONS * NCONS_VAR_SIZE
#define FLUX_SIZE NGROUPS * NSPECIES * NSTENCIL_FLUX * NCONS_VAR_SIZE
#define CMAX_SIZE NGROUPS * NSPECIES * NSTENCIL_CMAX

int ghl_m1_calc_flux(
    double timestep,
    const ghl_metric_quantities stencil_metric[NSTENCIL_FLUX],
    const ghl_ADM_aux_quantities stencil_adm_aux[NSTENCIL_FLUX],
    const ghl_primitive_quantities stencil_prims[NSTENCIL_FLUX],
    const ghl_m1_thc_params *thc_params,
    const double stencil_rE[STENCIL_CONS_SIZE],
    const ghl_radiation_flux_vector stencil_rF4[STENCIL_CONS_SIZE],
    const ghl_radiation_pressure_tensor stencil_rP4[STENCIL_CONS_SIZE],
    const double stencil_rN[STENCIL_CONS_SIZE],
    const double delta[3], // [dx,dy,dz]
    const double abs_1[2],
    const double scat_1[2],
    const int dir,
    ghl_m1_conservative_quantities *rhs_flux);

/*
 * Function      : calc_grsources.c functions
 * Description   : Calculates the gravitational sources
 * Documentation :
 */
int ghl_m1_calc_grsources(
    double timestep,
    const ghl_metric_quantities *metric,
    const ghl_metric_quantities *metric_derivs_x,
    const ghl_metric_quantities *metric_derivs_y,
    const ghl_metric_quantities *metric_derivs_z,
    const double E,
    const ghl_radiation_pressure_tensor *P4,
    const ghl_radiation_flux_vector *F4,
    const ghl_extrinsic_curvature *curv,
    ghl_m1_conservative_quantities *rhs_gr);

/*
 * Function      : calc_update.c functions
 * Description   : Calculates the update
 * Documentation :
 */
int ghl_m1_update(
    double cdt,
    const ghl_metric_quantities *metric,
    const ghl_extrinsic_curvature *curv,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_primitive_quantities *prims,
    const ghl_m1_thc_params *thc_params,
    const double abs_0[2],
    const double abs_1[2],
    const double eta_0[2],
    const double eta_1[2],
    const double scat_0[2],
    const double scat_1[2],
    ghl_m1_conservative_quantities *rU,
    ghl_m1_conservative_quantities *rU_p,
    ghl_m1_conservative_quantities *rhs_flux,
    ghl_m1_conservative_quantities *rhs_gr);

/*
 * Function      : m1_update_grid.c functions
 * Description   : Updates the grid
 * Documentation :
 */
int ghl_m1_update_grid(
    const int grid_Nx,
    const int grid_Ny,
    const int grid_Nz,
    const int imin,
    const int imax,
    const int jmin,
    const int jmax,
    const int kmin,
    const int kmax,
    const double timestep,
    const ghl_m1_thc_params *restrict thc_params,
    const double delta[3],
    const ghl_metric_quantities *restrict metric_grid,
    const ghl_extrinsic_curvature *restrict curv_grid,
    const ghl_ADM_aux_quantities *restrict adm_aux_grid,
    const ghl_primitive_quantities *restrict prims_grid,
    const double *restrict rE_grid,
    const ghl_radiation_flux_vector *restrict rF4_grid,
    const ghl_radiation_pressure_tensor *restrict rP4_grid,
    const double *restrict rN_grid,
    const double *restrict abs_1_grid,
    const double *restrict scat_1_grid,
    ghl_m1_conservative_quantities *restrict rU_grid,
    ghl_m1_conservative_quantities *restrict rU_p_grid);

/*
 * Function      : set_equilibrium.c functions
 * Description   : Sets the equilibrium state
 * Documentation :
 */
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





#endif // GHL_M1_H_
