#include "ghl.h"
#include "ghl_radiation.h"
#include "ghl_m1.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#define NGROUPS 1
#define NSPECIES 1

// This file translates from thc_M1_calc_grsources.cc
// The goal of this function is to input metric deriatives, E, F4, P4 and output the GR sources
// using functions provided in compute_radiation_flux.c for rhs from gr sources
int ghl_m1_calc_grsources(
    double timestep,
    const ghl_metric_quantities *metric,
    const ghl_metric_quantities *metric_derivs_x,
    const ghl_metric_quantities *metric_derivs_y,
    const ghl_metric_quantities *metric_derivs_z,
    const double E,
    const ghl_radiation_flux_vector *F4,
    const ghl_radiation_pressure_tensor *P4,
    const ghl_extrinsic_curvature *curv,
    ghl_m1_conservative_quantities *rhs_gr // output
    ){
    for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
        //int i4D_rhs = ig*5 + iv;
        double GE_source = 0.0;
        ghl_radiation_con_source_vector GF_source = {0};

        // Compute GR sources using functions from compute_radiation_flux.c 
        GE_source = calc_GE_source(metric, metric_derivs_x, metric_derivs_y, metric_derivs_z, P4, F4, curv);
        calc_GF_source(metric, metric_derivs_x, metric_derivs_y, metric_derivs_z, E, F4, P4, &GF_source);
        rhs_gr->rE[ig]  = GE_source;
        rhs_gr->rF4[ig].D[0] = 0.0; // No source for rF4[ig].D[0] in GR sources
        rhs_gr->rF4[ig].D[1] = GF_source.D[1];
        rhs_gr->rF4[ig].D[2] = GF_source.D[2];
        rhs_gr->rF4[ig].D[3] = GF_source.D[3];
        rhs_gr->rN[ig]  = 0.0; // No source for rN in GR sources
        
    } // for ig
    return 1;
} // ghl_m1_calc_grsources