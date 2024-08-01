#ifndef GHL_M1_H_
#define GHL_M1_H_

#include "ghl.h"
#include "radiation.h"

typedef struct {
  double DD[4][4];
} ghl_radiation_pressure_tensor; 

typedef struct {
  double UD[4][4];
} ghl_radiation_metric_tensor; 

typedef struct {
  double D[4];
} ghl_radiation_flux_vector; 

// conservative flux for M1 equations
typedef struct {
  double U[4];
} ghl_radiation_con_flux_vector; 

// conservative source for M1 equaitons
typedef struct {
  double U[4];
} ghl_radiation_con_source_vector; 

typedef struct {
  ghl_metric_quantities *metric;
  ghl_ADM_aux_quantities *adm_aux;
  ghl_primitive_quantities *prims;
  double E;
  ghl_radiation_flux_vector* F4;
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

#endif // GHL_M1_H_
