#ifndef UNIT_TESTS_H_
#define UNIT_TESTS_H_

#include "stdio.h"
#include "stdlib.h"
#include "con2prim.h"
#include "induction.h"

// con2prim validation functions
void validate_primitives(
      const int eos_type,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims,
      const primitive_quantities *restrict prims_trusted,
      const primitive_quantities *restrict prims_pert);

void validate_conservatives(
      const bool evolve_entropy,
      const conservative_quantities *restrict cons,
      const conservative_quantities *restrict cons_trusted,
      const conservative_quantities *restrict cons_pert);

void validate_stress_energy(
      const stress_energy *restrict Tmunu,
      const stress_energy *restrict Tmunu_trusted,
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
      double *restrict Bx,
      double *restrict By,
      double *restrict Bz,
      double *restrict entropy,
      double *restrict Y_e,
      double *restrict temperature,
      FILE *restrict outfile);

void read_conservative_binary(
      const bool evolve_entropy,
      double *restrict rho,
      double *restrict tau,
      double *restrict S_x,
      double *restrict S_y,
      double *restrict S_z,
      double *restrict entropy,
      FILE *restrict outfile);

void read_stress_energy_binary(
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

// Helper functions
static inline double relative_error( const double a, const double b ) {
  if     ( a != 0 ) return( fabs(1.0-b/a) );
  else if( b != 0 ) return( fabs(1.0-a/b) );
  else              return( 0.0 );
}

static inline double randf(double low,double high) {
  return (rand()/(double)(RAND_MAX))*(high-low)+low;
}

void initial_random_data(
      const double xrho,
      const double xpress,
      metric_quantities *restrict metric,
      primitive_quantities *restrict prims);

void randomize_metric(metric_quantities *restrict metric);

#define check_file_was_successfully_open(fp, filename) \
  if(!fp) grhayl_error("Could not open file %s.\n", filename);

#endif // UNIT_TESTS_H_
