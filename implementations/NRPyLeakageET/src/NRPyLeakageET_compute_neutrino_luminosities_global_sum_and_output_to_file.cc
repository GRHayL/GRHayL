#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"
#include "NRPyLeakageET.h"

extern "C"
void NRPyLeakageET_compute_neutrino_luminosities_global_sum_and_output_to_file(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( cctk_iteration%compute_luminosities_every==0 ) {
    if( verbosity_level > 1 ) CCTK_INFO("Computing luminosities on all refinement levels");
    // Step 1: Compute the neutrino luminosities
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
        BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {          
          NRPyLeakageET_compute_neutrino_luminosities(CCTK_PASS_CTOC);
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
    } END_REFLEVEL_LOOP;
    if( verbosity_level > 1 ) CCTK_INFO("Finished computing luminosities on all refinement levels");

    // Step 2: Perform the global integration
    // Step 2.a: Get sum operation handle
    const int op_sum = CCTK_ReductionHandle("sum");

    // Step 2.b: Set variables we want to reduce
    const char varnames[3][24] = {"NRPyLeakageET::lum_nue ",
                                  "NRPyLeakageET::lum_anue",
                                  "NRPyLeakageET::lum_nux "};


    // Step 2.c: Get indices for the variables we want to reduce
    const int varindices[3] = {CCTK_VarIndex("NRPyLeakageET::lum_nue"),
                               CCTK_VarIndex("NRPyLeakageET::lum_anue"),
                               CCTK_VarIndex("NRPyLeakageET::lum_nux")};

    // Step 2.d: Compute the global sum of the luminosities
    const CCTK_REAL d3x = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];
    CCTK_REAL global_luminosities[3];
    for(int i=0;i<3;i++) {

      if( verbosity_level > 1 ) CCTK_VINFO("Performing global reduction of gridfunction %s",varnames[i]);

      int varindex = varindices[i];
      assert(varindex>=0);
      CCTK_REAL global_luminosity = 0.0;
      const int ierr = CCTK_Reduce(cctkGH,
                                   -1, // target processors; -1 -> all
                                   op_sum,
                                   1,  // Number of outputs
                                   CCTK_VARIABLE_REAL,
                                   &global_luminosity,
                                   1,  // Number of inputs
                                   varindex);
      if( ierr ) CCTK_VERROR("Error in reduction of %s",varnames[i]);
      global_luminosities[i] = d3x * global_luminosity;

      if( verbosity_level > 1 ) CCTK_VINFO("Success! Reduction value: %e",global_luminosities[i]);
    }

    // Step 3: Now output the luminosities to file
    if( CCTK_MyProc(cctkGH) == 0 ) {
      char filename[512];
      sprintf(filename,"%s/%s",out_dir,luminosities_outfile);
      if( verbosity_level > 0 ) CCTK_VINFO("Outputting luminosities to file %s at iteration %d",filename,cctk_iteration);
      FILE *fp = fopen(filename,"a+");
      if( !fp ) CCTK_VERROR("Could not open file %s",filename);

      if( cctk_iteration == 0 ) {
        fprintf(fp,"# NRPyLeakageET output: Neutrino luminosities integrated over entire grid\n");
        fprintf(fp,"# Luminosities are given in geometrized units. Convert to\n");
        fprintf(fp,"# cgs units [erg/s] by multiplying by c^5 / G.\n");
        fprintf(fp,"# Column 1: cctk_iteration\n");
        fprintf(fp,"# Column 2: cctk_time\n");
        fprintf(fp,"# Column 3: Electron neutrino luminosity\n");
        fprintf(fp,"# Column 4: Electron antineutrino luminosity\n");
        fprintf(fp,"# Column 5: Heavy lepton neutrino (single species) luminosity (mult. by 4 to obtain total)\n");
      }

      fprintf(fp,"%d %.15e %.15e %.15e %.15e\n",cctk_iteration,cctk_time,global_luminosities[0],global_luminosities[1],global_luminosities[2]);
      fclose(fp);
      if( verbosity_level > 1 ) CCTK_INFO("Completed luminosities output");
    }
  }
}
