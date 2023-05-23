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

#include "GRHayLHDX.h"

/*******************************************************
 * Apply outer boundary conditions on {P,rho_b,vx,vy,vz}
 * It is better to apply BCs on primitives than conservs,
 * because small errors in conservs can be greatly
 * amplified in con2prim, sometimes leading to unphysical
 * primitives & unnecessary fixes.
 *******************************************************/
void GRHayLHDX_outer_boundaries_on_P_rho_b_vx_vy_vz(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHayLHDX_outer_boundaries_on_P_rho_b_vx_vy_vz;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Matter_BC,"frozen")) return;

  int levelnumber = std::ilogb(cctk_levfac[0]);
  printf("cctk_levfac %d reflevel %d\n", cctk_levfac[0], std::ilogb(cctk_levfac[0]));

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VERROR("ERROR: GRHayLHDX outer BC driver does not support unequal number of ghostzones in different directions!");

  // This sets the outer boundary for cell-centered variables, which have 1 less gridpoint.
  const int ccc_last[3] = {cctk_lsh[0]-1,
                           cctk_lsh[1]-1,
                           cctk_lsh[2]-1};
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {

    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) {
      const int imax=ccc_last[0]-cctk_nghostzones[0]+which_bdry_pt;
      for(int k=0; k<ccc_last[2]; k++)
        for(int j=0; j<ccc_last[1]; j++) {
          pressure[CCTK_GFINDEX3D(cctkGH,imax,j,k)] = pressure[CCTK_GFINDEX3D(cctkGH,imax-1,j,k)];
          rho_b[CCTK_GFINDEX3D(cctkGH,imax,j,k)]   = rho_b[CCTK_GFINDEX3D(cctkGH,imax-1,j,k)];
          vx[CCTK_GFINDEX3D(cctkGH,imax,j,k)]    = vx[CCTK_GFINDEX3D(cctkGH,imax-1,j,k)];
          vy[CCTK_GFINDEX3D(cctkGH,imax,j,k)]    = vy[CCTK_GFINDEX3D(cctkGH,imax-1,j,k)];
          vz[CCTK_GFINDEX3D(cctkGH,imax,j,k)]    = vz[CCTK_GFINDEX3D(cctkGH,imax-1,j,k)]; 
          if(vx[CCTK_GFINDEX3D(cctkGH,imax,j,k)]<0.) vx[CCTK_GFINDEX3D(cctkGH,imax,j,k)] = 0.0;
      }
    }
    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      const int imin=cctk_nghostzones[0]-which_bdry_pt-1;
      for(int k=0; k<ccc_last[2]; k++)
        for(int j=0; j<ccc_last[1]; j++) {
          pressure[CCTK_GFINDEX3D(cctkGH,imin,j,k)] = pressure[CCTK_GFINDEX3D(cctkGH,imin+1,j,k)];
          rho_b[CCTK_GFINDEX3D(cctkGH,imin,j,k)]   = rho_b[CCTK_GFINDEX3D(cctkGH,imin+1,j,k)];
          vx[CCTK_GFINDEX3D(cctkGH,imin,j,k)]    = vx[CCTK_GFINDEX3D(cctkGH,imin+1,j,k)];
          vy[CCTK_GFINDEX3D(cctkGH,imin,j,k)]    = vy[CCTK_GFINDEX3D(cctkGH,imin+1,j,k)];
          vz[CCTK_GFINDEX3D(cctkGH,imin,j,k)]    = vz[CCTK_GFINDEX3D(cctkGH,imin+1,j,k)]; 
          if(vx[CCTK_GFINDEX3D(cctkGH,imin,j,k)]>0.) vx[CCTK_GFINDEX3D(cctkGH,imin,j,k)] = 0.0;
      }
    }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) {
      const int jmax=ccc_last[1]-cctk_nghostzones[1]+which_bdry_pt;
      for(int k=0; k<ccc_last[2]; k++)
        for(int i=0; i<ccc_last[0]; i++) {
          pressure[CCTK_GFINDEX3D(cctkGH,i,jmax,k)] = pressure[CCTK_GFINDEX3D(cctkGH,i,jmax-1,k)];
          rho_b[CCTK_GFINDEX3D(cctkGH,i,jmax,k)]   = rho_b[CCTK_GFINDEX3D(cctkGH,i,jmax-1,k)];
          vx[CCTK_GFINDEX3D(cctkGH,i,jmax,k)]    = vx[CCTK_GFINDEX3D(cctkGH,i,jmax-1,k)];
          vy[CCTK_GFINDEX3D(cctkGH,i,jmax,k)]    = vy[CCTK_GFINDEX3D(cctkGH,i,jmax-1,k)];
          vz[CCTK_GFINDEX3D(cctkGH,i,jmax,k)]    = vz[CCTK_GFINDEX3D(cctkGH,i,jmax-1,k)]; 
          if(vx[CCTK_GFINDEX3D(cctkGH,i,jmax,k)]<0.) vx[CCTK_GFINDEX3D(cctkGH,i,jmax,k)] = 0.0;
      }
    }
    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      const int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
      for(int k=0; k<ccc_last[2]; k++)
        for(int i=0; i<ccc_last[0]; i++) {
          pressure[CCTK_GFINDEX3D(cctkGH,i,jmin,k)] = pressure[CCTK_GFINDEX3D(cctkGH,i,jmin+1,k)];
          rho_b[CCTK_GFINDEX3D(cctkGH,i,jmin,k)]   = rho_b[CCTK_GFINDEX3D(cctkGH,i,jmin+1,k)];
          vx[CCTK_GFINDEX3D(cctkGH,i,jmin,k)]    = vx[CCTK_GFINDEX3D(cctkGH,i,jmin+1,k)];
          vy[CCTK_GFINDEX3D(cctkGH,i,jmin,k)]    = vy[CCTK_GFINDEX3D(cctkGH,i,jmin+1,k)];
          vz[CCTK_GFINDEX3D(cctkGH,i,jmin,k)]    = vz[CCTK_GFINDEX3D(cctkGH,i,jmin+1,k)]; 
          if(vx[CCTK_GFINDEX3D(cctkGH,i,jmin,k)]>0.) vx[CCTK_GFINDEX3D(cctkGH,i,jmin,k)] = 0.0;
      }
    }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) {
      const int kmax=ccc_last[2]-cctk_nghostzones[2]+which_bdry_pt;
      for(int j=0; j<ccc_last[1]; j++)
        for(int i=0; i<ccc_last[0]; i++) {
          pressure[CCTK_GFINDEX3D(cctkGH,i,j,kmax)] = pressure[CCTK_GFINDEX3D(cctkGH,i,j,kmax-1)];
          rho_b[CCTK_GFINDEX3D(cctkGH,i,j,kmax)]   = rho_b[CCTK_GFINDEX3D(cctkGH,i,j,kmax-1)];
          vx[CCTK_GFINDEX3D(cctkGH,i,j,kmax)]    = vx[CCTK_GFINDEX3D(cctkGH,i,j,kmax-1)];
          vy[CCTK_GFINDEX3D(cctkGH,i,j,kmax)]    = vy[CCTK_GFINDEX3D(cctkGH,i,j,kmax-1)];
          vz[CCTK_GFINDEX3D(cctkGH,i,j,kmax)]    = vz[CCTK_GFINDEX3D(cctkGH,i,j,kmax-1)]; 
          if(vx[CCTK_GFINDEX3D(cctkGH,i,j,kmax)]<0.) vx[CCTK_GFINDEX3D(cctkGH,i,j,kmax)] = 0.0;
      }
    }
    // k=kmin=outer boundary
    if(cctk_bbox[4]) {
      const int kmin=cctk_nghostzones[2]-which_bdry_pt-1;
      for(int j=0; j<ccc_last[1]; j++)
        for(int i=0; i<ccc_last[0]; i++) {
          pressure[CCTK_GFINDEX3D(cctkGH,i,j,kmin)] = pressure[CCTK_GFINDEX3D(cctkGH,i,j,kmin+1)];
          rho_b[CCTK_GFINDEX3D(cctkGH,i,j,kmin)]   = rho_b[CCTK_GFINDEX3D(cctkGH,i,j,kmin+1)];
          vx[CCTK_GFINDEX3D(cctkGH,i,j,kmin)]    = vx[CCTK_GFINDEX3D(cctkGH,i,j,kmin+1)];
          vy[CCTK_GFINDEX3D(cctkGH,i,j,kmin)]    = vy[CCTK_GFINDEX3D(cctkGH,i,j,kmin+1)];
          vz[CCTK_GFINDEX3D(cctkGH,i,j,kmin)]    = vz[CCTK_GFINDEX3D(cctkGH,i,j,kmin+1)]; 
          if(vx[CCTK_GFINDEX3D(cctkGH,i,j,kmin)]>0.) vx[CCTK_GFINDEX3D(cctkGH,i,j,kmin)] = 0.0;
      }
    }
  }

  const double poison = 0.0/0.0;

  for(int k=0; k<ccc_last[2]; k++)
    for(int j=0; j<ccc_last[1]; j++)
      for(int i=0; i<ccc_last[0]; i++) {
        if(((cctk_bbox[0]) && i<cctk_nghostzones[0]) ||
           ((cctk_bbox[1]) && i>=ccc_last[0]-cctk_nghostzones[0]) ||
           ((cctk_bbox[2]) && j<cctk_nghostzones[1]) ||
           ((cctk_bbox[3]) && j>=ccc_last[1]-cctk_nghostzones[1]) ||
           ((cctk_bbox[4]) && k<cctk_nghostzones[2]) ||
           ((cctk_bbox[5]) && k>=ccc_last[2]-cctk_nghostzones[2])) {
          const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          metric_quantities ADM_metric;
          ghl_initialize_metric(
                ccc_lapse[index],
                ccc_betax[index], ccc_betay[index], ccc_betaz[index],
                ccc_gxx[index], ccc_gxy[index], ccc_gxz[index],
                ccc_gyy[index], ccc_gyz[index], ccc_gzz[index],
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
          stress_energy Tmunu;
          int speed_limited = 0;
          ghl_enforce_primitive_limits_and_compute_u0(
                ghl_params, ghl_eos, &ADM_metric,
                &metric_aux, &prims, &speed_limited);
          ghl_compute_conservs_and_Tmunu(
                &ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

          double dummy1, dummy2, dummy3;
          double dummy4, dummy5, dummy6;
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

          ghl_return_stress_energy(
                &Tmunu, &ccc_Ttt[index], &ccc_Ttx[index],
                &ccc_Tty[index], &ccc_Ttz[index], &ccc_Txx[index],
                &ccc_Txy[index], &ccc_Txz[index], &ccc_Tyy[index],
                &ccc_Tyz[index], &ccc_Tzz[index]);
        }
  }
}
