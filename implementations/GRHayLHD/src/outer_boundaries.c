/*******************************************************
 * Outer boundaries are handled as follows:
 * (-1) Update RHS quantities, leave RHS quantities zero on all outer ghostzones (including outer AMR refinement, processor, and outer boundaries)
 * ( 0) Let MoL update all evolution variables
 * ( 1) Apply outer boundary conditions (BCs) on A_{\mu}
 * ( 2) Compute B^i from A_i everywhere, synchronize B^i
 * ( 3) Call con2prim to get primitives on interior pts
 * ( 4) Apply outer BCs on {P,rho_b,vx,vy,vz}.
 * ( 5) (optional) set conservatives on outer boundary.
 *******************************************************/

#include "GRHayLHD.h"

void GRHayLHD_enforce_primitive_limits_and_compute_conservs(const cGH* cctkGH, const int index);

/*******************************************************
 * Apply outer boundary conditions on {P,rho_b,vx,vy,vz}
 * It is better to apply BCs on primitives than conservs,
 * because small errors in conservs can be greatly
 * amplified in con2prim, sometimes leading to unphysical
 * primitives & unnecessary fixes.
 *******************************************************/
void GRHayLHD_outer_boundaries_on_P_rho_b_vx_vy_vz(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_outer_boundaries_on_P_rho_b_vx_vy_vz;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Matter_BC,"frozen")) return;

  bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  int levelnumber = GetRefinementLevel(cctkGH);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VERROR("ERROR: GRHayLHD outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {

    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) {
      const int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int indm1 = CCTK_GFINDEX3D(cctkGH,imax-1, j, k);
          const int index = CCTK_GFINDEX3D(cctkGH,imax, j, k);
          pressure[index] = pressure[indm1];
          rho_b[index]    = rho_b[indm1];
          vx[index]       = vx[indm1];
          vy[index]       = vy[indm1];
          vz[index]       = vz[indm1]; 
          if(vx[index]<0.) vx[index] = 0.0;
          GRHayLHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index);
        }
      }
    }
    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      const int imin=cctk_nghostzones[0]-which_bdry_pt-1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH, imin, j, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, imin+1, j, k);
          pressure[index] = pressure[indp1];
          rho_b[index]    = rho_b[indp1];
          vx[index]       = vx[indp1];
          vy[index]       = vy[indp1];
          vz[index]       = vz[indp1]; 
          if(vx[index]>0.) vx[index] = 0.0;
          GRHayLHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index);
        }
      }
    }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) {
      const int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, jmax-1, k);
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmax, k);
          pressure[index] = pressure[indm1];
          rho_b[index]    = rho_b[indm1];
          vx[index]       = vx[indm1];
          vy[index]       = vy[indm1];
          vz[index]       = vz[indm1]; 
          if(vx[index]<0.) vx[index] = 0.0;
          GRHayLHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index);
        }
      }
    }
    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      const int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
#pragma omp parallel for
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmin, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, jmin+1, k);
          pressure[index] = pressure[indp1];
          rho_b[index]    = rho_b[indp1];
          vx[index]       = vx[indp1];
          vy[index]       = vy[indp1];
          vz[index]       = vz[indp1]; 
          if(vx[index]>0.) vx[index] = 0.0;
          GRHayLHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index);
        }
      }
    }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) {
      const int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, j, kmax-1);
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmax);
          pressure[index] = pressure[indm1];
          rho_b[index]    = rho_b[indm1];
          vx[index]       = vx[indm1];
          vy[index]       = vy[indm1];
          vz[index]       = vz[indm1]; 
          if(vx[index]<0.) vx[index] = 0.0;
          GRHayLHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index);
        }
      }
    }
    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      const int kmin=cctk_nghostzones[2]-which_bdry_pt-1;
#pragma omp parallel for
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmin);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, j, kmin+1);
          pressure[index] = pressure[indp1];
          rho_b[index]    = rho_b[indp1];
          vx[index]       = vx[indp1];
          vy[index]       = vy[indp1];
          vz[index]       = vz[indp1]; 
          if(vx[index]>0.) vx[index] = 0.0;
          GRHayLHD_enforce_primitive_limits_and_compute_conservs(cctkGH, index);
        }
      }
    }
  }
}

void GRHayLHD_enforce_primitive_limits_and_compute_conservs(const cGH* cctkGH, const int index) {
  // We cheat here by using the argument list of the scheduled function
  // instead of explicitly passing all these variables.
  DECLARE_CCTK_ARGUMENTS_GRHayLHD_outer_boundaries_on_P_rho_b_vx_vy_vz;

  const double poison = 0.0/0.0;
  double dummy1, dummy2, dummy3;
  double dummy4, dummy5, dummy6;
  int speed_limited = 0;

  metric_quantities ADM_metric;
  ghl_enforce_detgtij_and_initialize_ADM_metric(
        alp[index],
        betax[index], betay[index], betaz[index],
        gxx[index], gxy[index], gxz[index],
        gyy[index], gyz[index], gzz[index],
        &ADM_metric);

  ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

  primitive_quantities prims;
  ghl_initialize_primitives(
        rho_b[index], pressure[index], eps[index],
        vx[index], vy[index], vz[index],
        0.0, 0.0, 0.0,
        poison, poison, poison, &prims);

  conservative_quantities cons;
  ghl_enforce_primitive_limits_and_compute_u0(
        ghl_params, ghl_eos, &ADM_metric,
        &prims, &speed_limited);

  ghl_compute_conservs(
        &ADM_metric, &metric_aux, &prims, &cons);

  ghl_return_primitives(
        &prims,
        &rho_b[index], &pressure[index], &eps[index],
        &vx[index], &vy[index], &vz[index],
        &dummy1, &dummy2, &dummy3,
        &dummy4, &dummy5, &dummy6);

  ghl_return_conservatives(
        &cons,
        &rho_star[index], &tau[index],
        &Stildex[index], &Stildey[index], &Stildez[index],
        &dummy1, &dummy2);
}
