/**
 * @file NRPyEOS_stellarcollapse.h
 * @author Leo Werneck
 *
 * @brief Defines structures and functions for handling stellar collapse equation of
 * state (EOS) tables.
 */
#ifndef NRPYEOS_STELLARCOLLAPSE_H
#define NRPYEOS_STELLARCOLLAPSE_H

#include "ghl.h"
#include <stdbool.h>

#define SPEED_OF_LIGHT_SI          (299792458.0)
#define SPEED_OF_LIGHT_CGS         (SPEED_OF_LIGHT_SI * 100.0)
#define SPEED_OF_LIGHT_SQUARED_CGS (SPEED_OF_LIGHT_CGS * SPEED_OF_LIGHT_CGS)

/**
 * @brief Enum defining the various quantities available in the EOS table.
 */
typedef enum {
  NRPyEOS_sc_Abar,        ///< Average mass number A_h of heavy nuclei
  NRPyEOS_sc_Xa,          ///< Mass fraction of alpha particles
  NRPyEOS_sc_Xh,          ///< Mass fraction of heavy nuclei
  NRPyEOS_sc_Xn,          ///< Mass fraction of free neutrons
  NRPyEOS_sc_Xp,          ///< Mass fraction of free protons
  NRPyEOS_sc_Zbar,        ///< Average charge number Z_h of heavy nuclei
  NRPyEOS_sc_cs2,         ///< Sound speed squared, c_s^2
  NRPyEOS_sc_dedt,        ///< Derivative of specific internal energy w.r.t. temperature
  NRPyEOS_sc_dpderho,     ///< dP/deps = dpress/denergy, at constant density
  NRPyEOS_sc_dpdrhoe,     ///< dP/drho = dpress/drho, constant specific internal energy
  NRPyEOS_sc_entropy,     ///< Specific entropy per baryon
  NRPyEOS_sc_gamma,       ///< Adiabatic index
  NRPyEOS_sc_logenergy,   ///< Log10 of specific internal energy, log10(eps+energy_shift)
  NRPyEOS_sc_logpress,    ///< Log10 of pressure, log10(P)
  NRPyEOS_sc_mu_e,        ///< Electron chemical potential, mu_e
  NRPyEOS_sc_mu_n,        ///< Neutron chemical potential, mu_n
  NRPyEOS_sc_mu_p,        ///< Proton chemical potential, mu_p
  NRPyEOS_sc_muhat,       ///< Chemical potential difference, mu_hat = mu_n - mu_p
  NRPyEOS_sc_munu,        ///< Neutrino chemical potential, mu_nu = mu_n - mu_p + mu_e
  NRPyEOS_sc_n_quantities ///< Total number of EOS quantities
} NRPyEOS_stellarcollapse_quantity;

/**
 * @brief Structure representing a stellar collapse EOS table.
 *
 * Stores the grid dimensions, grid points, and data arrays for all EOS quantities.
 */
typedef struct {
  bool cs2_is_relativistic;  ///< Whether or not the sound speed is relativistic
  int n_rho;                 ///< Number of grid points in baryonic density, rho.
  int n_temperature;         ///< Number of grid points in temperature, T.
  int n_ye;                  ///< Number of grid points in electron fraction, Y_e.
  double *log10_rho;         ///< Array storing log10(rho).
  double *log10_temperature; ///< Array storing log10(T).
  double *ye;                ///< Array storing Y_e.
  double energy_shift;       ///< Energy shift applied to the specific internal energy.
  double *data[NRPyEOS_sc_n_quantities]; ///< Tabulated data.
} NRPyEOS_stellarcollapse_t;

/**
 * @brief Reads a stellar collapse EOS table from a file.
 *
 * Allocates memory for the table structure and reads the data from the specified HDF5
 * file.
 *
 * @param filepath Path to the EOS table file (HDF5 format).
 *
 * @return Pointer to the allocated and populated NRPyEOS_stellarcollapse_t structure, or
 * NULL on failure.
 */
NRPyEOS_stellarcollapse_t *NRPyEOS_stellarcollapse_read_table(const char *filepath);

/**
 * @brief Frees the memory allocated for a stellar collapse EOS table.
 *
 * Deallocates all memory associated with the table structure, including grid point
 * arrays and data arrays.
 *
 * @param table Pointer to the NRPyEOS_stellarcollapse_t structure to free.
 */
void NRPyEOS_stellarcollapse_free_table(NRPyEOS_stellarcollapse_t *table);

/**
 * @brief Returns a string version of the enumeration value.
 *
 * @param qty NRPyEOS_stellarcollapse_quantity enumeration value.
 */
char *NRPyEOS_stellarcollapse_qty_to_str(NRPyEOS_stellarcollapse_quantity qty);

/**
 * @brief Convert stellar collapse EOS struct to GRHayL EOS struct.
 *
 * @param sc Pointer to input NRPyEOS_stellarcollapse_t struct.
 * @param eos Pointer to output ghl_eos_parameters struct.
 */
void NRPyEOS_stellarcollapse_to_ghl(NRPyEOS_stellarcollapse_t *restrict sc, ghl_eos_parameters *restrict eos);

#endif // NRPYEOS_STELLARCOLLAPSE_H
