#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"

// mini NoMPI
#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#define BCAST(buffer, size) MPI_Bcast(buffer, size, MPI_BYTE, my_reader_process, MPI_COMM_WORLD)
#else
#define BCAST(buffer, size) do { /* do nothing */ } while(0)
#endif

// If on the IO proc (doIO == True) actually perform HDF5 IO, catch possible
// HDF5 errors
#define HDF5_DO_IO(fn_call)                                             \
  {                                                                     \
    int error_code = fn_call;                                           \
    if (error_code < 0) {                                               \
      fprintf(stderr,"(NRPyEOS) HDF5 call '%s' returned error code %d", \
                  #fn_call, error_code);                                \
    }                                                                   \
  }

/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_readtable_set_EOS_params(const char *nuceos_table_name, NRPyEOS_params *restrict eos_params) {


  fprintf(stderr,"(NRPyEOS) *******************************\n");
  fprintf(stderr,"(NRPyEOS) Reading EOS table from file:\n");
  fprintf(stderr,"(NRPyEOS) %s\n",nuceos_table_name);
  fprintf(stderr,"(NRPyEOS) *******************************\n");

  hid_t file;
  HDF5_DO_IO(file = H5Fopen(nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_BCAST_EOS_HDF5(NAME,VAR,TYPE,MEM,NELEMS)                   \
  do {                                                                  \
    hid_t dataset;                                                      \
    HDF5_DO_IO(dataset = H5Dopen(file, NAME, H5P_DEFAULT));             \
    HDF5_DO_IO(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR)); \
    BCAST (VAR, sizeof(*(VAR))*(NELEMS));                               \
    HDF5_DO_IO(H5Dclose(dataset));                                      \
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp array
#define READ_BCAST_EOSTABLE_HDF5(NAME,OFF,DIMS)                          \
  do {                                                                   \
    READ_BCAST_EOS_HDF5(NAME,&alltables_temp[(OFF)*(DIMS)[1]],H5T_NATIVE_DOUBLE,H5S_ALL,(DIMS)[1]); \
  } while (0)

  // Read size of tables
  READ_BCAST_EOS_HDF5("pointsrho",  &eos_params->nrho,  H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointstemp", &eos_params->ntemp, H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointsye",   &eos_params->nye,   H5T_NATIVE_INT, H5S_ALL, 1);

  // Allocate memory for tables
  double* alltables_temp;
  if (!(alltables_temp = (double*)malloc(eos_params->nrho * eos_params->ntemp * eos_params->nye * NRPyEOS_ntablekeys * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  if (!(eos_params->logrho = (double*)malloc(eos_params->nrho * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  if (!(eos_params->logtemp = (double*)malloc(eos_params->ntemp * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  if (!(eos_params->yes = (double*)malloc(eos_params->nye * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NRPyEOS_ntablekeys, (hsize_t)eos_params->nrho * eos_params->ntemp * eos_params->nye};
  hid_t mem3 =  H5Screate_simple(2, table_dims, NULL);

  // Read alltables_temp
  READ_BCAST_EOSTABLE_HDF5("logpress" , NRPyEOS_press_key  , table_dims);
  READ_BCAST_EOSTABLE_HDF5("logenergy", NRPyEOS_eps_key    , table_dims);
  READ_BCAST_EOSTABLE_HDF5("entropy"  , NRPyEOS_entropy_key, table_dims);
  READ_BCAST_EOSTABLE_HDF5("munu"     , NRPyEOS_munu_key   , table_dims);
  READ_BCAST_EOSTABLE_HDF5("cs2"      , NRPyEOS_cs2_key    , table_dims);
  READ_BCAST_EOSTABLE_HDF5("dedt"     , NRPyEOS_depsdT_key , table_dims);
  READ_BCAST_EOSTABLE_HDF5("dpdrhoe"  , NRPyEOS_dPdrho_key , table_dims);
  READ_BCAST_EOSTABLE_HDF5("dpderho"  , NRPyEOS_dPdeps_key , table_dims);
  // chemical potentials
  READ_BCAST_EOSTABLE_HDF5("muhat"    , NRPyEOS_muhat_key  , table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_e"     , NRPyEOS_mu_e_key   , table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_p"     , NRPyEOS_mu_p_key   , table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_n"     , NRPyEOS_mu_n_key   , table_dims);
  // compositions
  READ_BCAST_EOSTABLE_HDF5("Xa"       , NRPyEOS_X_a_key    , table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xh"       , NRPyEOS_X_h_key    , table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xn"       , NRPyEOS_X_n_key    , table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xp"       , NRPyEOS_X_p_key    , table_dims);
  // average nucleus
  READ_BCAST_EOSTABLE_HDF5("Abar"     , NRPyEOS_Abar_key   , table_dims);
  READ_BCAST_EOSTABLE_HDF5("Zbar"     , NRPyEOS_Zbar_key   , table_dims);
  // Gamma
  READ_BCAST_EOSTABLE_HDF5("gamma"    , NRPyEOS_Gamma_key  , table_dims);

  // Read additional tables and variables
  READ_BCAST_EOS_HDF5("logrho"      , eos_params->logrho       , H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->nrho);
  READ_BCAST_EOS_HDF5("logtemp"     , eos_params->logtemp      , H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->ntemp);
  READ_BCAST_EOS_HDF5("ye"          , eos_params->yes          , H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->nye);
  READ_BCAST_EOS_HDF5("energy_shift", &eos_params->energy_shift, H5T_NATIVE_DOUBLE, H5S_ALL, 1);

  HDF5_DO_IO(H5Sclose(mem3));
  HDF5_DO_IO(H5Fclose(file));

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  if (!(eos_params->alltables = (double*)malloc(eos_params->nrho * eos_params->ntemp * eos_params->nye * NRPyEOS_ntablekeys
                                                * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for EOS table");
  }
  for(int iv = 0;iv<NRPyEOS_ntablekeys;iv++)
    for(int k = 0; k<eos_params->nye;k++)
      for(int j = 0; j<eos_params->ntemp; j++)
        for(int i = 0; i<eos_params->nrho; i++) {
          int indold = i + eos_params->nrho*(j + eos_params->ntemp*(k + eos_params->nye*iv));
          int indnew = iv + NRPyEOS_ntablekeys*(i + eos_params->nrho*(j + eos_params->ntemp*k));
          eos_params->alltables[indnew] = alltables_temp[indold];
        }

  // free memory of temporary array
  free(alltables_temp);

  // convert units, convert logs to natural log
  // The latter is great, because exp() is way faster than pow()
  // pressure
  eos_params->energy_shift = eos_params->energy_shift * CGS_TO_CODE_ENERGY;
  for(int i=0;i<eos_params->nrho;i++) {
    // rewrite:
    // logrho[i] = log(pow(10.0,logrho[i]) * CGS_TO_CODE_DENSITY);
    // by using log(a^b*c) = b*log(a)+log(c)
    eos_params->logrho[i] = eos_params->logrho[i] * log(10.) + log(CGS_TO_CODE_DENSITY);
  }

  for(int i=0;i<eos_params->ntemp;i++) {
    // logtemp[i] = log(pow(10.0,logtemp[i]));
    eos_params->logtemp[i] = eos_params->logtemp[i]*log(10.0);
  }

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(eos_params->epstable = (double*)malloc(eos_params->nrho * eos_params->ntemp * eos_params->nye
                                               * sizeof(double)))) {
    fprintf(stderr,"(NRPyEOS) Cannot allocate memory for eps table\n");
  }

  // convert units
  int idx;
  for(int i=0;i<eos_params->nrho*eos_params->ntemp*eos_params->nye;i++) {
    // pressure
    idx = NRPyEOS_press_key + NRPyEOS_ntablekeys*i;
    eos_params->alltables[idx] = eos_params->alltables[idx] * log(10.0) + log(CGS_TO_CODE_PRESSURE);

    // eps
    idx = NRPyEOS_eps_key + NRPyEOS_ntablekeys*i;
    eos_params->alltables[idx] = eos_params->alltables[idx] * log(10.0) + log(CGS_TO_CODE_ENERGY);
    eos_params->epstable[i] = exp(eos_params->alltables[idx]);

    // cs2
    idx = NRPyEOS_cs2_key + NRPyEOS_ntablekeys*i;
    eos_params->alltables[idx] *= CGS_TO_CODE_LENGTH*CGS_TO_CODE_LENGTH/CGS_TO_CODE_TIME/CGS_TO_CODE_TIME;

    // dedT
    idx = NRPyEOS_depsdT_key + NRPyEOS_ntablekeys*i;
    eos_params->alltables[idx] *= CGS_TO_CODE_ENERGY;

    // dpdrhoe
    idx = NRPyEOS_dPdrho_key + NRPyEOS_ntablekeys*i;
    eos_params->alltables[idx] *= CGS_TO_CODE_PRESSURE/CGS_TO_CODE_DENSITY;

    // dpderho
    idx = NRPyEOS_dPdeps_key + NRPyEOS_ntablekeys*i;
    eos_params->alltables[idx] *= CGS_TO_CODE_PRESSURE/CGS_TO_CODE_ENERGY;
  }

  eos_params->temp0 = exp(eos_params->logtemp[0]);
  eos_params->temp1 = exp(eos_params->logtemp[1]);

  // set up some vars
  eos_params->dtemp  = (eos_params->logtemp[eos_params->ntemp-1] - eos_params->logtemp[0]) / (1.0*(eos_params->ntemp-1));
  eos_params->dtempi = 1.0/eos_params->dtemp;

  eos_params->dlintemp = eos_params->temp1-eos_params->temp0;
  eos_params->dlintempi = 1.0/eos_params->dlintemp;

  eos_params->drho  = (eos_params->logrho[eos_params->nrho-1] - eos_params->logrho[0]) / (1.0*(eos_params->nrho-1));
  eos_params->drhoi = 1.0/eos_params->drho;

  eos_params->dye  = (eos_params->yes[eos_params->nye-1] - eos_params->yes[0]) / (1.0*(eos_params->nye-1));
  eos_params->dyei = 1.0/eos_params->dye;

  eos_params->drhotempi      = eos_params->drhoi     * eos_params->dtempi;
  eos_params->drholintempi   = eos_params->drhoi     * eos_params->dlintempi;
  eos_params->drhoyei        = eos_params->drhoi     * eos_params->dyei;
  eos_params->dtempyei       = eos_params->dtempi    * eos_params->dyei;
  eos_params->dlintempyei    = eos_params->dlintempi * eos_params->dyei;
  eos_params->drhotempyei    = eos_params->drhoi     * eos_params->dtempi    * eos_params->dyei;
  eos_params->drholintempyei = eos_params->drhoi     * eos_params->dlintempi * eos_params->dyei;

  eos_params->eos_rhomax = exp(eos_params->logrho[eos_params->nrho-1]);
  eos_params->eos_rhomin = exp(eos_params->logrho[0]);

  eos_params->eos_tempmax = exp(eos_params->logtemp[eos_params->ntemp-1]);
  eos_params->eos_tempmin = exp(eos_params->logtemp[0]);

  eos_params->eos_yemax = eos_params->yes[eos_params->nye-1];
  eos_params->eos_yemin = eos_params->yes[0];
}
