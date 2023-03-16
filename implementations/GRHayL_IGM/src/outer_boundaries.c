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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "IGM.h"

#define IDX(i,j,k) CCTK_GFINDEX3D(cctkGH,(i),(j),(k))

#define XMAX_OB_LINEAR_EXTRAP(FUNC,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imax,j,k)] = 2.0 * FUNC[IDX(imax-1,j,k)] - FUNC[IDX(imax-2,j,k)];
#define YMAX_OB_LINEAR_EXTRAP(FUNC,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmax,k)] = 2.0 * FUNC[IDX(i,jmax-1,k)] - FUNC[IDX(i,jmax-2,k)];
#define ZMAX_OB_LINEAR_EXTRAP(FUNC,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmax)] = 2.0 * FUNC[IDX(i,j,kmax-1)] - FUNC[IDX(i,j,kmax-2)];

#define XMIN_OB_LINEAR_EXTRAP(FUNC,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imin,j,k)] = 2.0 * FUNC[IDX(imin+1,j,k)] - FUNC[IDX(imin+2,j,k)];
#define YMIN_OB_LINEAR_EXTRAP(FUNC,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmin,k)] = 2.0 * FUNC[IDX(i,jmin+1,k)] - FUNC[IDX(i,jmin+2,k)];
#define ZMIN_OB_LINEAR_EXTRAP(FUNC,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmin)] = 2.0 * FUNC[IDX(i,j,kmin+1)] - FUNC[IDX(i,j,kmin+2)];

/*********************************************
 * Apply outer boundary conditions on A_{\mu}
 ********************************************/
void GRHayL_IGM_outer_boundaries_on_A_mu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGM_outer_boundaries_on_A_mu;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(EM_BC,"frozen")) return;

  bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  const int levelnumber = GetRefinementLevel(cctkGH);

  //GRHayL_convert_ADM_to_BSSN(cctkGH, gxx, gxy, gxz, gyy, gyz, gzz,
  //                               phi_bssn,psi_bssn,
  //                               gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,
  //                               gtupxx, gtupxy, gtupxz, gtupyy, gtupyz, gtupzz);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VERROR("ERROR: GRHayL_IGM outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0; which_bdry_pt<cctk_nghostzones[0]; which_bdry_pt++) {
    const int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt; // for cctk_nghostzones==3, this goes {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    const int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
    const int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;

    const int imin=cctk_nghostzones[0]-which_bdry_pt-1; // for cctk_nghostzones==3, this goes {2,1,0}
    const int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
    const int kmin=cctk_nghostzones[2]-which_bdry_pt-1;

    if(cctk_bbox[1]) { XMAX_OB_LINEAR_EXTRAP(Ax,imax); XMAX_OB_LINEAR_EXTRAP(Ay,imax); XMAX_OB_LINEAR_EXTRAP(Az,imax); XMAX_OB_LINEAR_EXTRAP(phitilde,imax); }
    if(cctk_bbox[3]) { YMAX_OB_LINEAR_EXTRAP(Ax,jmax); YMAX_OB_LINEAR_EXTRAP(Ay,jmax); YMAX_OB_LINEAR_EXTRAP(Az,jmax); YMAX_OB_LINEAR_EXTRAP(phitilde,jmax); }
    if(cctk_bbox[5]) { ZMAX_OB_LINEAR_EXTRAP(Ax,kmax); ZMAX_OB_LINEAR_EXTRAP(Ay,kmax); ZMAX_OB_LINEAR_EXTRAP(Az,kmax); ZMAX_OB_LINEAR_EXTRAP(phitilde,kmax); }

    if(cctk_bbox[0]) {                    XMIN_OB_LINEAR_EXTRAP(Ax,imin); XMIN_OB_LINEAR_EXTRAP(Ay,imin); XMIN_OB_LINEAR_EXTRAP(Az,imin); XMIN_OB_LINEAR_EXTRAP(phitilde,imin); }
    if(cctk_bbox[2]) {                    YMIN_OB_LINEAR_EXTRAP(Ax,jmin); YMIN_OB_LINEAR_EXTRAP(Ay,jmin); YMIN_OB_LINEAR_EXTRAP(Az,jmin); YMIN_OB_LINEAR_EXTRAP(phitilde,jmin); }
    if((cctk_bbox[4]) && Symmetry_none) { ZMIN_OB_LINEAR_EXTRAP(Ax,kmin); ZMIN_OB_LINEAR_EXTRAP(Ay,kmin); ZMIN_OB_LINEAR_EXTRAP(Az,kmin); ZMIN_OB_LINEAR_EXTRAP(phitilde,kmin); }
  }
}

#define XMAX_OB_SIMPLE_COPY(FUNC,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imax,j,k)] = FUNC[IDX(imax-1,j,k)];
#define YMAX_OB_SIMPLE_COPY(FUNC,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmax,k)] = FUNC[IDX(i,jmax-1,k)];
#define ZMAX_OB_SIMPLE_COPY(FUNC,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmax)] = FUNC[IDX(i,j,kmax-1)];

#define XMIN_OB_SIMPLE_COPY(FUNC,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imin,j,k)] = FUNC[IDX(imin+1,j,k)];
#define YMIN_OB_SIMPLE_COPY(FUNC,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmin,k)] = FUNC[IDX(i,jmin+1,k)];
#define ZMIN_OB_SIMPLE_COPY(FUNC,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmin)] = FUNC[IDX(i,j,kmin+1)];


#define XMAX_INFLOW_CHECK(vx,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) if(vx[IDX(imax,j,k)]<0.) vx[IDX(imax,j,k)]=0.;
#define YMAX_INFLOW_CHECK(vy,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) if(vy[IDX(i,jmax,k)]<0.) vy[IDX(i,jmax,k)]=0.;
#define ZMAX_INFLOW_CHECK(vz,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) if(vz[IDX(i,j,kmax)]<0.) vz[IDX(i,j,kmax)]=0.;

#define XMIN_INFLOW_CHECK(vx,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) if(vx[IDX(imin,j,k)]>0.) vx[IDX(imin,j,k)]=0.;
#define YMIN_INFLOW_CHECK(vy,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) if(vy[IDX(i,jmin,k)]>0.) vy[IDX(i,jmin,k)]=0.;
#define ZMIN_INFLOW_CHECK(vz,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) if(vz[IDX(i,j,kmin)]>0.) vz[IDX(i,j,kmin)]=0.;


/*******************************************************
 * Apply outer boundary conditions on {P,rho_b,vx,vy,vz}
 * It is better to apply BCs on primitives than conservs,
 * because small errors in conservs can be greatly
 * amplified in con2prim, sometimes leading to unphysical
 * primitives & unnecessary fixes.
 *******************************************************/
void GRHayL_IGM_outer_boundaries_on_P_rho_b_vx_vy_vz(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHayL_IGM_outer_boundaries_on_P_rho_b_vx_vy_vz;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Matter_BC,"frozen")) return;

  bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  int levelnumber = GetRefinementLevel(cctkGH);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  int ENABLE=1;

  //GRHayL_IGM_convert_ADM_to_BSSN(cctkGH, gxx, gxy, gxz, gyy, gyz, gzz,
  //                               phi_bssn,psi_bssn,
  //                               gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,
  //                               gtupxx, gtupxy, gtupxz, gtupyy, gtupyz, gtupzz);

  //if(levelnumber<=11110) {
  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VERROR("ERROR: GRHayL_IGM outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {
    int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt; // for cctk_nghostzones==3, this goes {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
    int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;

    int imin=cctk_nghostzones[0]-which_bdry_pt-1; // for cctk_nghostzones==3, this goes {2,1,0}
    int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
    int kmin=cctk_nghostzones[2]-which_bdry_pt-1;

    // Order here is for compatibility with old version of this code.
    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) { XMAX_OB_SIMPLE_COPY(press,imax); XMAX_OB_SIMPLE_COPY(rho,imax); XMAX_OB_SIMPLE_COPY(vx,imax); XMAX_OB_SIMPLE_COPY(vy,imax); XMAX_OB_SIMPLE_COPY(vz,imax); if(ENABLE) XMAX_INFLOW_CHECK(vx,imax); }
    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      XMIN_OB_SIMPLE_COPY(press,imin); XMIN_OB_SIMPLE_COPY(rho,imin); XMIN_OB_SIMPLE_COPY(vx,imin); XMIN_OB_SIMPLE_COPY(vy,imin); XMIN_OB_SIMPLE_COPY(vz,imin); if(ENABLE) XMIN_INFLOW_CHECK(vx,imin); }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) { YMAX_OB_SIMPLE_COPY(press,jmax); YMAX_OB_SIMPLE_COPY(rho,jmax); YMAX_OB_SIMPLE_COPY(vx,jmax); YMAX_OB_SIMPLE_COPY(vy,jmax); YMAX_OB_SIMPLE_COPY(vz,jmax); if(ENABLE) YMAX_INFLOW_CHECK(vy,jmax); }
    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      YMIN_OB_SIMPLE_COPY(press,jmin); YMIN_OB_SIMPLE_COPY(rho,jmin); YMIN_OB_SIMPLE_COPY(vx,jmin); YMIN_OB_SIMPLE_COPY(vy,jmin); YMIN_OB_SIMPLE_COPY(vz,jmin); if(ENABLE) YMIN_INFLOW_CHECK(vy,jmin); }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) { ZMAX_OB_SIMPLE_COPY(press,kmax); ZMAX_OB_SIMPLE_COPY(rho,kmax); ZMAX_OB_SIMPLE_COPY(vx,kmax); ZMAX_OB_SIMPLE_COPY(vy,kmax); ZMAX_OB_SIMPLE_COPY(vz,kmax); if(ENABLE) ZMAX_INFLOW_CHECK(vz,kmax); }
    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      ZMIN_OB_SIMPLE_COPY(press,kmin); ZMIN_OB_SIMPLE_COPY(rho,kmin); ZMIN_OB_SIMPLE_COPY(vx,kmin); ZMIN_OB_SIMPLE_COPY(vy,kmin); ZMIN_OB_SIMPLE_COPY(vz,kmin); if(ENABLE) ZMIN_INFLOW_CHECK(vz,kmin); }
  }

  const double poison = 0.0/0.0;
  double dummy1, dummy2, dummy3;

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++)
    for(int j=0; j<cctk_lsh[1]; j++)
      for(int i=0; i<cctk_lsh[0]; i++) {
        if(((cctk_bbox[0]) && i<cctk_nghostzones[0]) ||
           ((cctk_bbox[1]) && i>=cctk_lsh[0]-cctk_nghostzones[0]) ||
           ((cctk_bbox[2]) && j<cctk_nghostzones[1]) ||
           ((cctk_bbox[3]) && j>=cctk_lsh[1]-cctk_nghostzones[1]) ||
           ((cctk_bbox[4]) && k<cctk_nghostzones[2] && CCTK_EQUALS(Symmetry,"none")) ||
           ((cctk_bbox[5]) && k>=cctk_lsh[2]-cctk_nghostzones[2])) {
          const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          metric_quantities metric;
          GRHayL_enforce_detgtij_and_initialize_metric(
                alp[index],
                gxx[index], gxy[index], gxz[index],
                gyy[index], gyz[index], gzz[index],
                betax[index], betay[index], betaz[index],
                &metric);
      
          primitive_quantities prims;
          initialize_primitives(rho[index],
                press[index], eps[index],
                vx[index], vy[index], vz[index],
                Bx_center[index], By_center[index], Bz_center[index],
                poison, poison, poison, &prims);
                //wont have storage for these vars for hybrid
                //entropy[index], Y_e[index], temperature[index], &prims);
      
          conservative_quantities cons;
          stress_energy Tmunu;
          int speed_limited = 0;
          enforce_primitive_limits_and_compute_u0(grhayl_params, grhayl_eos, &metric, &prims, &speed_limited);
          compute_conservs_and_Tmunu(grhayl_params, grhayl_eos, &metric, &prims, &cons, &Tmunu);
      
          return_primitives(&prims,
                &rho[index], &press[index], &eps[index],
                &vx[index], &vy[index], &vz[index],
                //&Bvec[index0], &Bvec[index1], &Bvec[index2],
                &Bx_center[index], &By_center[index], &Bz_center[index],
                &dummy1, &dummy2, &dummy3);
                //wont have storage for these vars for hybrid
                //&entropy[index], &Y_e[index], &temperature[index]);
      
          return_conservatives(&cons,
                &rho_star[index], &tau[index],
                &Stildex[index], &Stildey[index], &Stildez[index],
                &dummy1, &dummy2);
                //&Y_e[index], &entropy[index]);
      
          if(grhayl_params->update_Tmunu) {
            return_stress_energy(&Tmunu, &eTtt[index], &eTtx[index],
                  &eTty[index], &eTtz[index], &eTxx[index],
                  &eTxy[index], &eTxz[index], &eTyy[index],
                  &eTyz[index], &eTzz[index]);
          }

          //if(i==5 && j==5 && k==5) CCTK_VInfo(CCTK_THORNSTRING,"%e %e %e %e",eTtt[index],eTtx[index],eTty[index],eTxy[index]);
          //CCTK_VInfo(CCTK_THORNSTRING,"YAY: "); for(ww=0;ww<10;ww++) CCTK_VInfo(CCTK_THORNSTRING,"%e ",TDNMUNU[ww]); CCTK_VInfo(CCTK_THORNSTRING,"");
        }
  }
}
