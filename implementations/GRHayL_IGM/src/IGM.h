#include "GRHayLET.h"

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in GRHayLET_evaluate_MHD_rhs.C.
static const int RHOB=0,PRESSURE=1,VX=2,VY=3,VZ=4,
  BX_CENTER=5,BY_CENTER=6,BZ_CENTER=7,BX_STAGGER=8,BY_STAGGER=9,BZ_STAGGER=10,
  VXR=11,VYR=12,VZR=13,VXL=14,VYL=15,VZL=16,MAXNUMVARS=17;  //<-- Be _sure_ to define MAXNUMVARS appropriately!

