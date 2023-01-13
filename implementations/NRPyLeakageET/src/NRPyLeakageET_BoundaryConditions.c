#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define IDX(i,j,k) CCTK_GFINDEX3D(cctkGH,(i),(j),(k))

#define XMAX_OB_ZERO(FUNC,imax) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imax,j,k)] = 0.0;
#define YMAX_OB_ZERO(FUNC,jmax) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmax,k)] = 0.0;
#define ZMAX_OB_ZERO(FUNC,kmax) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmax)] = 0.0;

#define XMIN_OB_ZERO(FUNC,imin) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imin,j,k)] = 0.0;
#define YMIN_OB_ZERO(FUNC,jmin) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmin,k)] = 0.0;
#define ZMIN_OB_ZERO(FUNC,kmin) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmin)] = 0.0;


/*********************************************
 * Apply outer boundary conditions on A_{\mu}
 ********************************************/
void NRPyLeakageET_BoundaryConditions(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int levelnumber = GetRefinementLevel(cctkGH);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VERROR("ERROR: NRPyLeakageET outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {
    int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt; // for cctk_nghostzones==3, this goes {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
    int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;

    int imin=cctk_nghostzones[0]-which_bdry_pt-1; // for cctk_nghostzones==3, this goes {2,1,0}
    int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
    int kmin=cctk_nghostzones[2]-which_bdry_pt-1;

    if(cctk_bbox[1]) {
      XMAX_OB_ZERO(tau_0_nue  ,imax);XMAX_OB_ZERO(tau_0_anue  ,imax);XMAX_OB_ZERO(tau_0_nux  ,imax);XMAX_OB_ZERO(tau_1_nue  ,imax);XMAX_OB_ZERO(tau_1_anue  ,imax);XMAX_OB_ZERO(tau_1_nux  ,imax);
      XMAX_OB_ZERO(kappa_0_nue,imax);XMAX_OB_ZERO(kappa_0_anue,imax);XMAX_OB_ZERO(kappa_0_nux,imax);XMAX_OB_ZERO(kappa_1_nue,imax);XMAX_OB_ZERO(kappa_1_anue,imax);XMAX_OB_ZERO(kappa_1_nux,imax);
    }
    if(cctk_bbox[3]) {
      YMAX_OB_ZERO(tau_0_nue  ,jmax);YMAX_OB_ZERO(tau_0_anue  ,jmax);YMAX_OB_ZERO(tau_0_nux  ,jmax);YMAX_OB_ZERO(tau_1_nue  ,jmax);YMAX_OB_ZERO(tau_1_anue  ,jmax);YMAX_OB_ZERO(tau_1_nux  ,jmax);
      YMAX_OB_ZERO(kappa_0_nue,jmax);YMAX_OB_ZERO(kappa_0_anue,jmax);YMAX_OB_ZERO(kappa_0_nux,jmax);YMAX_OB_ZERO(kappa_1_nue,jmax);YMAX_OB_ZERO(kappa_1_anue,jmax);YMAX_OB_ZERO(kappa_1_nux,jmax);
    }
    if(cctk_bbox[5]) {
      ZMAX_OB_ZERO(tau_0_nue  ,kmax);ZMAX_OB_ZERO(tau_0_anue  ,kmax);ZMAX_OB_ZERO(tau_0_nux  ,kmax);ZMAX_OB_ZERO(tau_1_nue  ,kmax);ZMAX_OB_ZERO(tau_1_anue  ,kmax);ZMAX_OB_ZERO(tau_1_nux  ,kmax);
      ZMAX_OB_ZERO(kappa_0_nue,kmax);ZMAX_OB_ZERO(kappa_0_anue,kmax);ZMAX_OB_ZERO(kappa_0_nux,kmax);ZMAX_OB_ZERO(kappa_1_nue,kmax);ZMAX_OB_ZERO(kappa_1_anue,kmax);ZMAX_OB_ZERO(kappa_1_nux,kmax);
    }

    if(cctk_bbox[0]) {
      XMIN_OB_ZERO(tau_0_nue  ,imin);XMIN_OB_ZERO(tau_0_anue  ,imin);XMIN_OB_ZERO(tau_0_nux  ,imin);XMIN_OB_ZERO(tau_1_nue  ,imin);XMIN_OB_ZERO(tau_1_anue  ,imin);XMIN_OB_ZERO(tau_1_nux  ,imin);
      XMIN_OB_ZERO(kappa_0_nue,imin);XMIN_OB_ZERO(kappa_0_anue,imin);XMIN_OB_ZERO(kappa_0_nux,imin);XMIN_OB_ZERO(kappa_1_nue,imin);XMIN_OB_ZERO(kappa_1_anue,imin);XMIN_OB_ZERO(kappa_1_nux,imin);
    }
    if(cctk_bbox[2]) {
      YMIN_OB_ZERO(tau_0_nue  ,jmin);YMIN_OB_ZERO(tau_0_anue  ,jmin);YMIN_OB_ZERO(tau_0_nux  ,jmin);YMIN_OB_ZERO(tau_1_nue  ,jmin);YMIN_OB_ZERO(tau_1_anue  ,jmin);YMIN_OB_ZERO(tau_1_nux  ,jmin);
      YMIN_OB_ZERO(kappa_0_nue,jmin);YMIN_OB_ZERO(kappa_0_anue,jmin);YMIN_OB_ZERO(kappa_0_nux,jmin);YMIN_OB_ZERO(kappa_1_nue,jmin);YMIN_OB_ZERO(kappa_1_anue,jmin);YMIN_OB_ZERO(kappa_1_nux,jmin);
    }
    if(cctk_bbox[4]) {
      ZMIN_OB_ZERO(tau_0_nue  ,kmin);ZMIN_OB_ZERO(tau_0_anue  ,kmin);ZMIN_OB_ZERO(tau_0_nux  ,kmin);ZMIN_OB_ZERO(tau_1_nue  ,kmin);ZMIN_OB_ZERO(tau_1_anue  ,kmin);ZMIN_OB_ZERO(tau_1_nux  ,kmin);
      ZMIN_OB_ZERO(kappa_0_nue,kmin);ZMIN_OB_ZERO(kappa_0_anue,kmin);ZMIN_OB_ZERO(kappa_0_nux,kmin);ZMIN_OB_ZERO(kappa_1_nue,kmin);ZMIN_OB_ZERO(kappa_1_anue,kmin);ZMIN_OB_ZERO(kappa_1_nux,kmin);
    }
  }
}
