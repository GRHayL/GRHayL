#ifndef UNIT_TESTS_H_
#define UNIT_TESTS_H_

#include "stdio.h"
#include "stdlib.h"
#include "con2prim.h"
#include "induction.h"

// con2prim validation functions
int validate_primitives(
      const double tolerance,
      const int eos_type,
      const bool velocity_only,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims_orig,
      const primitive_quantities *restrict prims,
      FILE *restrict infile);

int validate_conservatives(
      const double tolerance,
      const bool evolve_entropy,
      const conservative_quantities *restrict cons_orig,
      const conservative_quantities *restrict cons,
      FILE *restrict infile);

int validate_stress_energy(
      const double tolerance,
      const stress_energy *restrict Tmunu_orig,
      const stress_energy *restrict Tmunu,
      FILE *restrict infile);

// con2prim binary input functions
void read_primitive_binary(
      const int eos_type,
      const bool velocity_only,
      const bool evolve_entropy,
      primitive_quantities *restrict prims, 
      FILE *restrict outfile);

void read_conservative_binary(
      const bool evolve_entropy,
      conservative_quantities *restrict cons_orig,
      conservative_quantities *restrict cons,
      FILE *restrict outfile);

void read_stress_energy_binary(
      stress_energy *restrict Tmunu_orig,
      stress_energy *restrict Tmunu,
      FILE *restrict outfile);

void read_metric_binary(
      metric_quantities *restrict metric,
      FILE *restrict infile);

// con2prim binary output functions
void write_primitive_binary(
      const int eos_type,
      const bool velocity_only,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims, 
      FILE *restrict outfile);

void write_conservative_binary(
      const bool evolve_entropy,
      const conservative_quantities *restrict cons_orig,
      const conservative_quantities *restrict cons,
      FILE *restrict outfile);

void write_stress_energy_binary(
      const stress_energy *restrict Tmunu_orig,
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

static int rel_tol(const double tolerance, const double x1, const double x2) {
  const double rel_diff = relative_error(x1, x2);
  if(rel_diff > tolerance) return 1;
  return 0;
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
  if( fp == NULL ) { \
    fprintf(stderr, "(GRHayL) ERROR: Could not open file %s. Terminating.\n", filename); \
    exit(1); \
  }

#endif // UNIT_TESTS_H_
