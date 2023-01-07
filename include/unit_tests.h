#include "stdio.h"
#include "stdlib.h"
#include "con2prim.h"

void con2prim_unit_test( );

// con2prim binary output functions
void write_primitive_binary(
      const int eos_type,
      const bool velocity_only,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims_orig, 
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
      const bool random_metric,
      metric_quantities *restrict metric,
      primitive_quantities *restrict prims);
