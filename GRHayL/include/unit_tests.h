#ifndef UNIT_TESTS_H_
#define UNIT_TESTS_H_

#include "atmosphere.h"
#include "con2prim.h"
#include "induction.h"
#include "reconstruction.h"
#include "flux_source.h"
#include "neutrinos.h"
#include "nrpyeos_tabulated.h"
#include "nrpyeos_hybrid.h"

void ghl_test_compute_A_flux_with_B(
      const int dirlength,
      const int A_dir,
      const double *restrict phi_bssn,
      const double *restrict cmin_1,
      const double *restrict cmax_1,
      const double *restrict cmin_2,
      const double *restrict cmax_2,
      const double *restrict v1rr,
      const double *restrict v1rl,
      const double *restrict v1lr,
      const double *restrict v1ll,
      const double *restrict v2rr,
      const double *restrict v2rl,
      const double *restrict v2lr,
      const double *restrict v2ll,
      const double *restrict B1r,
      const double *restrict B1l,
      const double *restrict B2r,
      const double *restrict B2l,
      double *restrict A_rhs);

void ghl_test_compute_A_flux_with_Btilde(
      const int dirlength,
      const int A_dir,
      const double *restrict cmin_1,
      const double *restrict cmax_1,
      const double *restrict cmin_2,
      const double *restrict cmax_2,
      const double *restrict v1rr,
      const double *restrict v1rl,
      const double *restrict v1lr,
      const double *restrict v1ll,
      const double *restrict v2rr,
      const double *restrict v2rl,
      const double *restrict v2lr,
      const double *restrict v2ll,
      const double *restrict B1r,
      const double *restrict B1l,
      const double *restrict B2r,
      const double *restrict B2l,
      double *restrict A_rhs);

void ghl_test_compute_ccc_ADM(
      const int dirlength,
      const double *restrict lapse,
      const double *restrict betax,
      const double *restrict betay,
      const double *restrict betaz,
      const double *restrict gxx,
      const double *restrict gxy,
      const double *restrict gxz,
      const double *restrict gyy,
      const double *restrict gyz,
      const double *restrict gzz,
      const double *restrict phitilde,
      const double *restrict Ax,
      const double *restrict Ay,
      const double *restrict Az,
      double *restrict alpha_interp,
      double *restrict betax_interp,
      double *restrict betay_interp,
      double *restrict betaz_interp,
      double *restrict alpha_Phi_minus_betaj_A_j_interp,
      double *restrict sqrtg_Ax_interp,
      double *restrict sqrtg_Ay_interp,
      double *restrict sqrtg_Az_interp);
    
void ghl_test_compute_vvv_ADM(
      const int dirlength,
      const double *restrict lapse,
      const double *restrict betax,
      const double *restrict betay,
      const double *restrict betaz,
      const double *restrict gxx,
      const double *restrict gxy,
      const double *restrict gxz,
      const double *restrict gyy,
      const double *restrict gyz,
      const double *restrict gzz,
      const double *restrict phitilde,
      const double *restrict Ax,
      const double *restrict Ay,
      const double *restrict Az,
      double *restrict alpha_Phi_minus_betaj_A_j_interp,
      double *restrict sqrtg_Ax_interp,
      double *restrict sqrtg_Ay_interp,
      double *restrict sqrtg_Az_interp);

void ghl_test_compute_ccc_BSSN(
      const int dirlength,
      const double *restrict lapse,
      const double *restrict betax,
      const double *restrict betay,
      const double *restrict betaz,
      const double *restrict psi,
      const double *restrict gtupxx,
      const double *restrict gtupxy,
      const double *restrict gtupxz,
      const double *restrict gtupyy,
      const double *restrict gtupyz,
      const double *restrict gtupzz,
      const double *restrict phitilde,
      const double *restrict Ax,
      const double *restrict Ay,
      const double *restrict Az,
      double *restrict alpha_interp,
      double *restrict betax_interp,
      double *restrict betay_interp,
      double *restrict betaz_interp,
      double *restrict alpha_Phi_minus_betaj_A_j_interp,
      double *restrict sqrtg_Ax_interp,
      double *restrict sqrtg_Ay_interp,
      double *restrict sqrtg_Az_interp);

void ghl_test_compute_h_and_cs2(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      double *restrict h,
      double *restrict cs2);

// con2prim validation functions
void ghl_pert_test_fail_primitives(
      const bool evolve_entropy,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims_trusted,
      const ghl_primitive_quantities *restrict prims,
      const ghl_primitive_quantities *restrict prims_pert);

void ghl_pert_test_fail_primitives_with_cutoffs(
      const bool evolve_entropy,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims_trusted,
      const ghl_primitive_quantities *restrict prims,
      const ghl_primitive_quantities *restrict prims_pert,
      const double pressure_cutoff,
      const double eps_cutoff);

void ghl_pert_test_fail_conservatives(
      const bool evolve_entropy,
      const ghl_conservative_quantities *restrict cons_trusted,
      const ghl_conservative_quantities *restrict cons,
      const ghl_conservative_quantities *restrict cons_pert);

void ghl_pert_test_fail_stress_energy(
      const ghl_stress_energy *restrict Tmunu_trusted,
      const ghl_stress_energy *restrict Tmunu,
      const ghl_stress_energy *restrict Tmunu_pert);

// con2prim binary input functions
void read_primitive_binary(
      const int eos_type,
      const bool evolve_entropy,
      double *restrict rho,
      double *restrict press,
      double *restrict vx,
      double *restrict vy,
      double *restrict vz,
      double *restrict eps,
      double *restrict Bx,
      double *restrict By,
      double *restrict Bz,
      double *restrict entropy,
      double *restrict Y_e,
      double *restrict temperature,
      FILE *restrict outfile);

void read_primitive_struct_binary(
      const int eos_type,
      const bool evolve_entropy,
      ghl_primitive_quantities *restrict prims,
      FILE *restrict infile);

void read_conservative_binary(
      const bool evolve_entropy,
      double *restrict rho,
      double *restrict tau,
      double *restrict S_x,
      double *restrict S_y,
      double *restrict S_z,
      double *restrict entropy,
      FILE *restrict outfile);

void read_conservative_struct_binary(
      const bool evolve_entropy,
      ghl_conservative_quantities *restrict cons,
      FILE *restrict infile);

void read_stress_energy_binary(
      double *restrict Ttt,
      double *restrict Ttx,
      double *restrict Tty,
      double *restrict Ttz,
      double *restrict Txx,
      double *restrict Txy,
      double *restrict Txz,
      double *restrict Tyy,
      double *restrict Tyz,
      double *restrict Tzz,
      FILE *restrict outfile);

void read_stress_energy_struct_binary(
      ghl_stress_energy *restrict Tmunu,
      FILE *restrict outfile);

void read_metric_binary(
      double *restrict lapse,
      double *restrict gxx,
      double *restrict gxy,
      double *restrict gxz,
      double *restrict gyy,
      double *restrict gyz,
      double *restrict gzz,
      double *restrict betax,
      double *restrict betay,
      double *restrict betaz,
      FILE *restrict infile);

void read_metric_struct_binary(
      ghl_metric_quantities *restrict metric,
      FILE *restrict infile);

// con2prim binary output functions
void write_primitive_binary(
      const int eos_type,
      const bool evolve_entropy,
      const ghl_primitive_quantities *restrict prims,
      FILE *restrict outfile);

void write_conservative_binary(
      const bool evolve_entropy,
      const ghl_conservative_quantities *restrict cons,
      FILE *restrict outfile);

void write_stress_energy_binary(
      const ghl_stress_energy *restrict Tmunu,
      FILE *restrict outfile);

void write_metric_binary(
      const ghl_metric_quantities *restrict metric,
      FILE *restrict outfile);

// Tabulated EOS helper functions
double get_table_quantity(
      const int which_var,
      const double rho,
      const double Y_e,
      const double T);

// Helper functions
static inline double relative_error(const double a, const double b) {
  if      (a != 0) return fabs(1.0-b/a);
  else if (b != 0) return fabs(1.0-a/b);
  else             return 0.0;
}

static inline int indexf(const int gridmax, const int i, const int j, const int k) {
  return i + j*gridmax + k*gridmax*gridmax;
}

static inline bool ghl_pert_test_fail_with_tolerance(
      const double trusted,
      const double computed,
      const double perturbed,
      const double rel_tol,
      const double abs_tol) {
  if (isnan(computed) && isfinite(trusted)) return true; // NaN failure
  if (fabs(trusted - computed) < abs_tol) return false;  // Absolute tolerance success
  if (isnan(perturbed)) return false; // NaN "success"
  return relative_error(trusted, computed) > fmax(4.0*relative_error(trusted, perturbed), rel_tol);
}

static inline bool ghl_pert_test_fail(
      const double trusted,
      const double computed,
      const double perturbed) {
  const double min_rel = 8.0e-14;
  const double min_abs = 1.0e-30;
  return ghl_pert_test_fail_with_tolerance(trusted, computed, perturbed, min_rel, min_abs);
}

static inline double randf(double low,double high) {
  return (rand()/(double)(RAND_MAX))*(high-low)+low;
}

void ghl_initial_random_data(
      const double xrho,
      const double xpress,
      ghl_metric_quantities *restrict metric,
      ghl_primitive_quantities *restrict prims);

void ghl_randomize_metric(
      double *restrict lapse,
      double *restrict gxx_ptr,
      double *restrict gxy_ptr,
      double *restrict gxz_ptr,
      double *restrict gyy_ptr,
      double *restrict gyz_ptr,
      double *restrict gzz_ptr,
      double *restrict betax,
      double *restrict betay,
      double *restrict betaz);

void ghl_randomize_primitives(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double press,
      double *restrict vx,
      double *restrict vy,
      double *restrict vz,
      double *restrict Bx,
      double *restrict By,
      double *restrict Bz);

#define check_file_was_successfully_open(fp, filename) \
  if(!fp) ghl_error("Could not open file %s.\n", filename);

static inline
FILE *
fopen_with_check(const char *filename, const char *mode) {
  FILE *fp = fopen(filename, mode);
  if(!fp) ghl_error("Could not open file %s.\n", filename);
  return fp;
}
#endif // UNIT_TESTS_H_
