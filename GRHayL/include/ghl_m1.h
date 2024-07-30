#ifndef GHL_M1_H_
#define GHL_M1_H_

#include "ghl.h"
#include "radiation.h"

typedef struct {
  double DD[4][4];
  double UD[4][4];
} ghl_radiation_pressure_tensor; 

typedef struct {
  double UD[4][4];
} ghl_radiation_metric_tensor; 

typedef struct {
  double D[4];
} ghl_radiation_flux_vector; 

typedef struct {
  double U[4];
} ghl_radiation_source_vector; 

typedef struct {
  const ghl_metric_quantities *metric;
  const ghl_ADM_aux_quantities *adm_aux;
  const ghl_primitive_quantities *prims;
  const double E;
  const ghl_radiation_flux_vector* F4;
  const ghl_radiation_pressure_tensor *P4DD;
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
