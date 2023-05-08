#ifndef UNIT_TESTS_H_
#define UNIT_TESTS_H_

#include "con2prim.h"
#include "induction.h"
#include "reconstruction.h"
#include "flux_source.h"

// con2prim validation functions
void validate_primitives(
      const bool evolve_entropy,
      const eos_parameters *restrict eos,
      const primitive_quantities *restrict prims_trusted,
      const primitive_quantities *restrict prims,
      const primitive_quantities *restrict prims_pert);

void validate_conservatives(
      const bool evolve_entropy,
      const conservative_quantities *restrict cons_trusted,
      const conservative_quantities *restrict cons,
      const conservative_quantities *restrict cons_pert);

void validate_stress_energy(
      const stress_energy *restrict Tmunu_trusted,
      const stress_energy *restrict Tmunu,
      const stress_energy *restrict Tmunu_pert);

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
      primitive_quantities *restrict prims,
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
      conservative_quantities *restrict cons,
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
      stress_energy *restrict Tmunu,
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
      metric_quantities *restrict metric,
      FILE *restrict infile);

// con2prim binary output functions
void write_primitive_binary(
      const int eos_type,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims,
      FILE *restrict outfile);

void write_conservative_binary(
      const bool evolve_entropy,
      const conservative_quantities *restrict cons,
      FILE *restrict outfile);

void write_stress_energy_binary(
      const stress_energy *restrict Tmunu,
      FILE *restrict outfile);

void write_metric_binary(
      const metric_quantities *restrict metric,
      FILE *restrict outfile);

// Tabulated EOS helper functions
double get_table_quantity(
      const int which_var,
      const double rho,
      const double Y_e,
      const double T );

// Helper functions
static inline double relative_error( const double a, const double b ) {
  if     ( a != 0 ) return( fabs(1.0-b/a) );
  else if( b != 0 ) return( fabs(1.0-a/b) );
  else              return( 0.0 );
}

inline int indexf(const int gridmax, const int i, const int j, const int k) {
  return i + j*gridmax + k*gridmax*gridmax;
}

static inline bool validate_with_tolerance(
      const double trusted,
      const double computed,
      const double perturbed,
      const double rel_tol,
      const double abs_tol) {
  if (fabs(trusted - computed) < abs_tol) return false;
  return relative_error(trusted, computed) > fmax(4.0*relative_error(trusted, perturbed), rel_tol);
}

static inline bool validate(
      const double trusted,
      const double computed,
      const double perturbed) {
  const double min_rel = 8.0e-14;
  const double min_abs = 1.0e-30;
  return validate_with_tolerance(trusted, computed, perturbed, min_rel, min_abs);
}

static inline double randf(double low,double high) {
  return (rand()/(double)(RAND_MAX))*(high-low)+low;
}

void initial_random_data(
      const double xrho,
      const double xpress,
      metric_quantities *restrict metric,
      primitive_quantities *restrict prims);

void randomize_metric(
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

void randomize_primitives(
      const eos_parameters *restrict eos,
      const double rho,
      const double press,
      double *restrict eps,
      double *restrict vx,
      double *restrict vy,
      double *restrict vz,
      double *restrict Bx,
      double *restrict By,
      double *restrict Bz);

#define check_file_was_successfully_open(fp, filename) \
  if(!fp) grhayl_error("Could not open file %s.\n", filename);

static inline
FILE *
fopen_with_check(const char *filename, const char *mode) {
  FILE *fp = fopen(filename, mode);
  if(!fp) grhayl_error("Could not open file %s.\n", filename);
  return fp;
}
#endif // UNIT_TESTS_H_
