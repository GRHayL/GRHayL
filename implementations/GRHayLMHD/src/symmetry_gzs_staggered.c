#include "GRHayLMHD.h"

void GRHayLMHD_set_symmetry_gzs_staggered(
      const cGH *cctkGH,
      const double *X,
      const double *Y,
      const double *Z,
      double *gridfunc,
      const double *gridfunc_syms,
      const int stagger_x,  //TODO: unused
      const int stagger_y,  //TODO: unused
      const int stagger_z) {

  DECLARE_CCTK_PARAMETERS;

  //const int imax = cctkGH->cctk_lsh[0];
  //const int jmax = cctkGH->cctk_lsh[1];
  //const int kmax = cctkGH->cctk_lsh[2];
  const int lsh[3] = {cctkGH->cctk_lsh[0], cctkGH->cctk_lsh[1], cctkGH->cctk_lsh[2]};

  if(CCTK_EQUALS(Symmetry, "equatorial"))
    CCTK_VERROR("Warning: Symmetry==equatorial not supported! USE AT YOUR OWN RISK. You will need to comment this error message out.");

  // No symmetries -> return.
  if(CCTK_EQUALS(Symmetry, "none")) return;

  CCTK_REAL dz = Z[CCTK_GFINDEX3D(cctkGH,0,0,1)] - Z[CCTK_GFINDEX3D(cctkGH,0,0,0)];

  CCTK_REAL z_offset = dz*0.5*stagger_z;

  int num_gzs=0;
  //FIXME: Might want to use cctk_nghostzones instead...
  while( (Z[CCTK_GFINDEX3D(cctkGH,0,0,num_gzs)]+z_offset) < -dz*0.1 && num_gzs<lsh[2]) num_gzs++;
  if(num_gzs*2>=lsh[2]) CCTK_VERROR("ERROR in GRHayLMHD_set_symmetry_gzs_staggered.c");
  //while( (Z[CCTK_GFINDEX3D(cctkGH,0,0,num_gzs)]+z_offset) < -dz*0.1 && num_gzs<kmax) num_gzs++;
  //if(num_gzs*2>=kmax) CCTK_VERROR("ERROR in GRHayLMHD_set_symmetry_gzs_staggered.c");

#pragma omp parallel for
  for(int k=0; k<num_gzs; k++) {
    for(int j=0; j<lsh[1]; j++) {
      for(int i=0; i<lsh[0]; i++) {
        const int index_inside__sym_gz = CCTK_GFINDEX3D(cctkGH,i,j,k);

        /* This loop sets symmetry ghostzones, regardless of how the gridfunction is staggered.
         *
         * STAGGERED PATTERN:
         * if num_gzs==1 && stagger_z==1:
         * z[] = {-dz/2,dz/2,3dz/2, etc} -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 1]
         *
         * if num_gzs==2 && stagger_z==1:
         * z[] = {-3dz/2,-dz/2,dz/2,3dz/2 etc}
         * -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 3]
         * -> gridfunc[index 1] = gridfunc_syms[2]*gridfunc[index 2]
         * .
         * .
         * .
         * -> gridfunc[i] = gridfunc_syms[2]*gridfunc[(num_gz*2-1)-i]
         *
         * UNSTAGGERED PATTERN:
         * if num_gzs==1 && stagger_z==0:
         * z[] = {-dz,0,dz, etc} -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 2]
         *
         * if num_gzs==2 && stagger_z==0:
         * z[] = {-2dz,-dz,0,dz,2dz, etc} -> gridfunc[index 0] = gridfunc_syms[2]*gridfunc[index 4]
         * z[] = {-2dz,-dz,0,dz,2dz, etc} -> gridfunc[index 1] = gridfunc_syms[2]*gridfunc[index 3]
         * .
         * .
         * .
         * -> gridfunc[i] = gridfunc_syms[2]*gridfunc[(num_gz*2)-i]
         *
         * OVERALL PATTERN: gridfunc[i] = gridfunc_syms[2]*gridfunc[(num_gz*2-stagger_z)-i] */

        const int matching_index_outside_sym_gz = CCTK_GFINDEX3D(cctkGH,i,j,(num_gzs*2-stagger_z)-k);

        gridfunc[index_inside__sym_gz] = gridfunc_syms[2]*gridfunc[matching_index_outside_sym_gz];
      }
    }
  }
}
