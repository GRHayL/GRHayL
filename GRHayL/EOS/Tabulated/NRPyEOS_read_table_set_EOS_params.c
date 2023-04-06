#include "NRPyEOS_Tabulated.h"

// mini NoMPI
#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#define BCAST(buffer, size) MPI_Bcast(buffer, size, MPI_BYTE, my_reader_process, MPI_COMM_WORLD)
#else
#define BCAST(buffer, size) do { /* do nothing */ } while(0)
#endif

// If on the IO proc (doIO == True) actually perform HDF5 IO, catch possible
// HDF5 errors
#define HDF5_DO_IO(fn_call)                                   \
  {                                                           \
    int error_code = fn_call;                                 \
    if (error_code < 0) {                                     \
      grhayl_error("HDF5 call '%s' returned error code %d\n", \
                   #fn_call, error_code);                     \
    }                                                         \
  }

#define check_if_file_exists(filename) {                     \
  FILE *fp = fopen(filename, "r");                           \
  if( !fp )                                                  \
    grhayl_error("Could not open EOS file: %s\n", filename); \
  fclose(fp);                                                \
}

static inline double get_EOS_table_max(
      const eos_parameters *restrict eos,
      const int var_key ) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Loop over the table, searching for the maximum value
  int totalsize        = eos->N_rho * eos->N_Ye * eos->N_T;
  double var_max_value = eos->table_all[var_key];

  for(int i=0;i<totalsize;i++) {
    double var_value = eos->table_all[var_key + NRPyEOS_ntablekeys*i];
    if( var_value > var_max_value ) var_max_value = var_value;
  }
  return var_max_value;
#endif
}

static inline double get_EOS_table_min(
      const eos_parameters *restrict eos,
      const int var_key ) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Loop over the table, searching for the minimum value
  int totalsize        = eos->N_rho * eos->N_Ye * eos->N_T;
  double var_min_value = eos->table_all[var_key];

  for(int i=0;i<totalsize;i++) {
    double var_value = eos->table_all[var_key + NRPyEOS_ntablekeys*i];
    if( var_value < var_min_value ) var_min_value = var_value;
  }
  return var_min_value;
#endif
}
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_read_table_set_EOS_params(const char *EOS_tablename, eos_parameters *restrict eos_params) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  check_if_file_exists(EOS_tablename);

  grhayl_info("*******************************\n");
  grhayl_info("Reading EOS table from file:\n");
  grhayl_info("%s\n", EOS_tablename);
  grhayl_info("*******************************\n");

  hid_t file = H5Fopen(EOS_tablename, H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_DO_IO(file);

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
  READ_BCAST_EOS_HDF5("pointsrho",  &eos_params->N_rho, H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointstemp", &eos_params->N_T,   H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointsye",   &eos_params->N_Ye,  H5T_NATIVE_INT, H5S_ALL, 1);

  // Allocate memory for tables
  double* alltables_temp;
  if (!(alltables_temp = (double*)malloc(eos_params->N_rho * eos_params->N_T * eos_params->N_Ye * NRPyEOS_ntablekeys * sizeof(double))))
    grhayl_error("Cannot allocate memory for EOS table\n");

  if (!(eos_params->table_logrho = (double*)malloc(eos_params->N_rho * sizeof(double))))
    grhayl_error("Cannot allocate memory for EOS table\n");

  if (!(eos_params->table_logT = (double*)malloc(eos_params->N_T * sizeof(double))))
    grhayl_error("Cannot allocate memory for EOS table\n");

  if (!(eos_params->table_Ye = (double*)malloc(eos_params->N_Ye * sizeof(double))))
    grhayl_error("Cannot allocate memory for EOS table\n");


  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NRPyEOS_ntablekeys, (hsize_t)eos_params->N_rho * eos_params->N_T * eos_params->N_Ye};
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
  READ_BCAST_EOS_HDF5("logrho"      , eos_params->table_logrho , H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->N_rho);
  READ_BCAST_EOS_HDF5("logtemp"     , eos_params->table_logT   , H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->N_T);
  READ_BCAST_EOS_HDF5("ye"          , eos_params->table_Ye     , H5T_NATIVE_DOUBLE, H5S_ALL, eos_params->N_Ye);
  READ_BCAST_EOS_HDF5("energy_shift", &eos_params->energy_shift, H5T_NATIVE_DOUBLE, H5S_ALL, 1);

  HDF5_DO_IO(H5Sclose(mem3));
  HDF5_DO_IO(H5Fclose(file));

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  if (!(eos_params->table_all = (double*)malloc(eos_params->N_rho * eos_params->N_T * eos_params->N_Ye * NRPyEOS_ntablekeys
                                                * sizeof(double)))) {
    grhayl_info("Cannot allocate memory for EOS table");
  }
  for(int iv = 0;iv<NRPyEOS_ntablekeys;iv++)
    for(int k = 0; k<eos_params->N_Ye;k++)
      for(int j = 0; j<eos_params->N_T; j++)
        for(int i = 0; i<eos_params->N_rho; i++) {
          int indold = i + eos_params->N_rho*(j + eos_params->N_T*(k + eos_params->N_Ye*iv));
          int indnew = iv + NRPyEOS_ntablekeys*(i + eos_params->N_rho*(j + eos_params->N_T*k));
          eos_params->table_all[indnew] = alltables_temp[indold];
        }

  // free memory of temporary array
  free(alltables_temp);

  // convert units, convert logs to natural log
  // The latter is great, because exp() is way faster than pow()
  // pressure
  eos_params->energy_shift = eos_params->energy_shift * CGS_TO_CODE_ENERGY;
  for(int i=0;i<eos_params->N_rho;i++) {
    // rewrite:
    // logrho[i] = log(pow(10.0,logrho[i]) * CGS_TO_CODE_DENSITY);
    // by using log(a^b*c) = b*log(a)+log(c)
    eos_params->table_logrho[i] = eos_params->table_logrho[i] * log(10.) + log(CGS_TO_CODE_DENSITY);
  }

  for(int i=0;i<eos_params->N_T;i++) {
    // logtemp[i] = log(pow(10.0,logtemp[i]));
    eos_params->table_logT[i] = eos_params->table_logT[i]*log(10.0);
  }

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(eos_params->table_eps = (double*)malloc(eos_params->N_rho * eos_params->N_T * eos_params->N_Ye
                                               * sizeof(double))))
    grhayl_error("Cannot allocate memory for eps table\n");

  // convert units
  int idx;
  for(int i=0;i<eos_params->N_rho*eos_params->N_T*eos_params->N_Ye;i++) {
    // pressure
    idx = NRPyEOS_press_key + NRPyEOS_ntablekeys*i;
    eos_params->table_all[idx] = eos_params->table_all[idx] * log(10.0) + log(CGS_TO_CODE_PRESSURE);

    // eps
    idx = NRPyEOS_eps_key + NRPyEOS_ntablekeys*i;
    eos_params->table_all[idx] = eos_params->table_all[idx] * log(10.0) + log(CGS_TO_CODE_ENERGY);
    eos_params->table_eps[i] = exp(eos_params->table_all[idx]);

    // cs2
    idx = NRPyEOS_cs2_key + NRPyEOS_ntablekeys*i;
    eos_params->table_all[idx] *= CGS_TO_CODE_LENGTH*CGS_TO_CODE_LENGTH/CGS_TO_CODE_TIME/CGS_TO_CODE_TIME;

    // dedT
    idx = NRPyEOS_depsdT_key + NRPyEOS_ntablekeys*i;
    eos_params->table_all[idx] *= CGS_TO_CODE_ENERGY;

    // dpdrhoe
    idx = NRPyEOS_dPdrho_key + NRPyEOS_ntablekeys*i;
    eos_params->table_all[idx] *= CGS_TO_CODE_PRESSURE/CGS_TO_CODE_DENSITY;

    // dpderho
    idx = NRPyEOS_dPdeps_key + NRPyEOS_ntablekeys*i;
    eos_params->table_all[idx] *= CGS_TO_CODE_PRESSURE/CGS_TO_CODE_ENERGY;
  }

  eos_params->temp0 = exp(eos_params->table_logT[0]);
  eos_params->temp1 = exp(eos_params->table_logT[1]);

  // set up some vars
  eos_params->dtemp  = (eos_params->table_logT[eos_params->N_T-1] - eos_params->table_logT[0]) / (1.0*(eos_params->N_T-1));
  eos_params->dtempi = 1.0/eos_params->dtemp;

  eos_params->dlintemp = eos_params->temp1-eos_params->temp0;
  eos_params->dlintempi = 1.0/eos_params->dlintemp;

  eos_params->drho  = (eos_params->table_logrho[eos_params->N_rho-1] - eos_params->table_logrho[0]) / (1.0*(eos_params->N_rho-1));
  eos_params->drhoi = 1.0/eos_params->drho;

  eos_params->dye  = (eos_params->table_Ye[eos_params->N_Ye-1] - eos_params->table_Ye[0]) / (1.0*(eos_params->N_Ye-1));
  eos_params->dyei = 1.0/eos_params->dye;

  eos_params->drhotempi      = eos_params->drhoi     * eos_params->dtempi;
  eos_params->drholintempi   = eos_params->drhoi     * eos_params->dlintempi;
  eos_params->drhoyei        = eos_params->drhoi     * eos_params->dyei;
  eos_params->dtempyei       = eos_params->dtempi    * eos_params->dyei;
  eos_params->dlintempyei    = eos_params->dlintempi * eos_params->dyei;
  eos_params->drhotempyei    = eos_params->drhoi     * eos_params->dtempi    * eos_params->dyei;
  eos_params->drholintempyei = eos_params->drhoi     * eos_params->dlintempi * eos_params->dyei;

  eos_params->table_rho_max  = exp(eos_params->table_logrho[eos_params->N_rho-1]);
  eos_params->table_rho_min  = exp(eos_params->table_logrho[0]);
  eos_params->table_T_max    = exp(eos_params->table_logT[eos_params->N_T-1]);
  eos_params->table_T_min    = exp(eos_params->table_logT[0]);
  eos_params->table_Ye_max   = eos_params->table_Ye[eos_params->N_Ye-1];
  eos_params->table_Ye_min   = eos_params->table_Ye[0];

  // Table bounds for useful quantities
  eos_params->table_P_min   = exp(get_EOS_table_min(eos_params, NRPyEOS_press_key));
  eos_params->table_P_max   = exp(get_EOS_table_max(eos_params, NRPyEOS_press_key));
  eos_params->table_eps_min = exp(get_EOS_table_min(eos_params, NRPyEOS_eps_key)) - eos_params->energy_shift;
  eos_params->table_eps_max = exp(get_EOS_table_max(eos_params, NRPyEOS_eps_key)) - eos_params->energy_shift;
  eos_params->table_ent_min = get_EOS_table_min(eos_params, NRPyEOS_entropy_key);
  eos_params->table_ent_max = get_EOS_table_max(eos_params, NRPyEOS_entropy_key);
#endif
}
