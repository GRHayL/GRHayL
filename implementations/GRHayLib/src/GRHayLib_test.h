/**
 * @file GRHayLib_test.h
 * @author Leo Werneck
 * @email wernecklr@gmail.com
 * @brief Testing utilities and macros for GRHayLib
 *
 * This header file provides a comprehensive set of macros for outputting quantities
 * in any Einstein Toolkit thorn in a format that is GRHayL-friendly. Functionality
 * includes logging, random number generation, and metric/primitive variable manipulation
 * for unit tests and validation.
 *
 * @section macros Defined Macros
 *
 * **Must have macros**
 * - GHL_TEST_LOG_START - Initialize binary log file for current function
 * - GHL_TEST_LOG_END - Close and finalize log file
 * - GHL_TEST_LOG_WRITE_METRIC_AND_PRIMS(index_) - Log both metric and primitives
 * - GHL_TEST_LOG_WRITE(var) - Write variable to binary log file
 *
 * **Logging and Output:**
 * - GHL_PRINTF(...) - Printf wrapper with GHL_TEST prefix
 * - GHL_TEST_LOG_START - Initialize binary log file for current function
 * - GHL_TEST_LOG_END - Close and finalize log file
 * - GHL_TEST_LOG_WRITE(var) - Write variable to binary log file
 * - GHL_TEST_LOG_WRITE_METRIC_AND_PRIMS(index_) - Log both metric and primitives
 * - GHL_TEST_LOG_METRIC(index_) - Log metric quantities to file
 * - GHL_TEST_LOG_PRIMS(index_) - Log primitive variables to file
 *
 * **Random Number Generation:**
 * - GHL_TEST_RAND_IN_RANGE(a, b) - Generate random double between a and b
 *
 * **Metric Manipulation:**
 * - GHL_TEST_SET_METRIC_TO_FLAT(index_) - Set metric to flat spacetime
 * - GHL_TEST_SET_METRIC_TO_RANDOM(index_) - Set metric to random values
 *
 * **Primitive Variable Manipulation:**
 * - GHL_TEST_SET_PRIMS_TO_RANDOM(index_) - Set primitives to random values
 *
 * @note Random number generation requires proper seeding with srand().
 * @warning The random metric generation may not produce physically reasonable
 *          extrinsic curvature values in all cases.
 */

#ifndef GRHAYLIB_TEST_H
#define GRHAYLIB_TEST_H

#include <stdio.h>
#include <stdlib.h>

#include "GRHayLib.h"

#define GHL_TEST_MAGIC_NUMBER (4224)

#define GHL_PRINTF(...)    \
    printf("(GHL TEST) "); \
    printf(__VA_ARGS__);

// Create log file. This must appear at the beginning of every tested function.
#define GHL_TEST_LOG_START                                          \
    char ghl_test_filename_[1024];                                  \
    sprintf(ghl_test_filename_, "ghl_test_%s.bin", __func__);       \
    FILE *ghl_test_fp_ = fopen(ghl_test_filename_, "wb");           \
    if(ghl_test_fp_ == NULL) {                                      \
        CCTK_VERROR("Failed to open file %s", ghl_test_filename_);  \
    }                                                               \
    GHL_PRINTF("Succesfully opened file %s\n", ghl_test_filename_);

// Close log file. This must appear at the end of every tested function.
#define GHL_TEST_LOG_END                                             \
    fclose(ghl_test_fp_);                                            \
    GHL_PRINTF("Successfully closed file %s\n", ghl_test_filename_); \
    GHL_PRINTF("Terminating run\n");                                 \
    exit(0);

// Generic macro to output "anything" to file. Uses decltype instead of typeof.
#define GHL_TEST_LOG_WRITE(var) fwrite(&var, sizeof(decltype(var)), 1, ghl_test_fp_);

// Write a magic number, for sanity checks
#define GHL_TEST_LOG_WRITE_MAGIC_NUMBER                           \
    {                                                             \
        const int ghl_test_magic_number_ = GHL_TEST_MAGIC_NUMBER; \
        GHL_TEST_LOG_WRITE(ghl_test_magic_number_);               \
    }

// Generate random number between a and b
#define GHL_TEST_RAND_IN_RANGE(a, b) ((a) + ((double)rand() / RAND_MAX) * ((b) - (a)))

// Begin ET-specific stuff
// Log GRHayL structs for metric quantities, auxiliary quantities, and extrinsic curvature
#define GHL_TEST_LOG_METRIC(index_)                                            \
    ghl_metric_quantities   ghl_test_adm_metric_ = {0};                        \
    ghl_ADM_aux_quantities  ghl_test_aux_metric_ = {0};                        \
    ghl_extrinsic_curvature ghl_test_curv_       = {0};                        \
    ghl_initialize_metric(                                                     \
        alp[index_],                                                           \
        betax[index_],                                                         \
        betay[index_],                                                         \
        betaz[index_],                                                         \
        gxx[index_],                                                           \
        gxy[index_],                                                           \
        gxz[index_],                                                           \
        gyy[index_],                                                           \
        gyz[index_],                                                           \
        gzz[index_],                                                           \
        &ghl_test_adm_metric_                                                  \
    );                                                                         \
    ghl_compute_ADM_auxiliaries(&ghl_test_adm_metric_, &ghl_test_aux_metric_); \
    ghl_initialize_extrinsic_curvature(                                        \
        kxx[index_],                                                           \
        kxy[index_],                                                           \
        kxz[index_],                                                           \
        kyy[index_],                                                           \
        kyz[index_],                                                           \
        kzz[index_],                                                           \
        &ghl_test_curv_                                                        \
    );                                                                         \
    GHL_TEST_LOG_WRITE(ghl_test_adm_metric_);                                  \
    GHL_TEST_LOG_WRITE(ghl_test_aux_metric_);                                  \
    GHL_TEST_LOG_WRITE(ghl_test_curv_);

// Log GRHayL struct for primitives
#define GHL_TEST_LOG_PRIMS(index_)                                                                             \
    const CCTK_INT           ghl_test_size_  = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];                        \
    ghl_primitive_quantities ghl_test_prims_ = {0};                                                            \
    ghl_initialize_primitives(                                                                                 \
        rho[index_],                                                                                           \
        press[index_],                                                                                         \
        eps[index_],                                                                                           \
        vel[index_ + 0 * ghl_test_size_],                                                                      \
        vel[index_ + 1 * ghl_test_size_],                                                                      \
        vel[index_ + 2 * ghl_test_size_],                                                                      \
        Bvec[index_ + 0 * ghl_test_size_],                                                                     \
        Bvec[index_ + 1 * ghl_test_size_],                                                                     \
        Bvec[index_ + 2 * ghl_test_size_],                                                                     \
        entropy[index_],                                                                                       \
        Y_e[index_],                                                                                           \
        temperature[index_],                                                                                   \
        &ghl_test_prims_                                                                                       \
    );                                                                                                         \
    bool ghl_test_speed_limited_ = false;                                                                      \
    ghl_limit_v_and_compute_u0(ghl_params, &ghl_test_adm_metric_, &ghl_test_prims_, &ghl_test_speed_limited_); \
    GHL_TEST_LOG_WRITE(ghl_test_prims_);

// Set metric to flat space
#define GHL_TEST_SET_METRIC_TO_FLAT(index_) \
    alp[index_]   = 1.0;                    \
    betax[index_] = 0.0;                    \
    betay[index_] = 0.0;                    \
    betaz[index_] = 0.0;                    \
    gxx[index_]   = 1.0;                    \
    gxy[index_]   = 0.0;                    \
    gxz[index_]   = 0.0;                    \
    gyy[index_]   = 1.0;                    \
    gyz[index_]   = 0.0;                    \
    gzz[index_]   = 1.0;                    \
    kxx[index_]   = 0.0;                    \
    kxy[index_]   = 0.0;                    \
    kxz[index_]   = 0.0;                    \
    kyy[index_]   = 0.0;                    \
    kyz[index_]   = 0.0;                    \
    kzz[index_]   = 0.0;

// Randomize metric. Note(LRW): we compute gxx enforcing det(gamma_ij) = e^(12 phi).
// Warning(LRW): not sure the extrinsic curvature values are sane.
#define GHL_TEST_SET_METRIC_TO_RANDOM(index_)                                                                          \
    const double ghl_test_gyy_ = 1.0 + GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                              \
    const double ghl_test_gzz_ = 1.0 + GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                              \
    const double ghl_test_gxy_ = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                    \
    const double ghl_test_gxz_ = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                    \
    const double ghl_test_gyz_ = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                    \
    const double ghl_test_phi_ = GHL_TEST_RAND_IN_RANGE(0.0, 0.2);                                                     \
    const double ghl_test_gxx_ =                                                                                       \
        (exp(12.0 * ghl_test_phi_) - ghl_test_gxy_ * (ghl_test_gyz_ * ghl_test_gxz_ - ghl_test_gxy_ * ghl_test_gzz_) - \
         ghl_test_gxz_ * (ghl_test_gxy_ * ghl_test_gyz_ - ghl_test_gyy_ * ghl_test_gxz_)) /                            \
        (ghl_test_gyy_ * ghl_test_gzz_ - ghl_test_gyz_ * ghl_test_gyz_);                                               \
    alp[index_]   = GHL_TEST_RAND_IN_RANGE(0.1, 1);                                                                    \
    betax[index_] = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    betay[index_] = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    betaz[index_] = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    gxx[index_]   = ghl_test_gxx_;                                                                                     \
    gxy[index_]   = ghl_test_gxy_;                                                                                     \
    gxz[index_]   = ghl_test_gxz_;                                                                                     \
    gyy[index_]   = ghl_test_gyy_;                                                                                     \
    gyz[index_]   = ghl_test_gyz_;                                                                                     \
    gzz[index_]   = ghl_test_gzz_;                                                                                     \
    kxx[index_]   = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    kxy[index_]   = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    kxz[index_]   = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    kyy[index_]   = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    kyz[index_]   = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);                                                                 \
    kzz[index_]   = GHL_TEST_RAND_IN_RANGE(-0.1, 0.1);

// Randomize primitives. Note only rho, T, and Y_e are randomized, the rest are
// computed using the EOS to ensure thermodynamic sanity.
#define GHL_TEST_SET_PRIMS_TO_RANDOM(index_)                          \
    rho[index_]         = pow(10.0, GHL_TEST_RAND_IN_RANGE(-12, -3)); \
    temperature[index_] = pow(10.0, GHL_TEST_RAND_IN_RANGE(-2, 2));   \
    Y_e[index_]         = GHL_TEST_RAND_IN_RANGE(0.1, 0.5);           \
    ghl_tabulated_compute_P_eps_S_from_T(                             \
        ghl_eos,                                                      \
        rho[index_],                                                  \
        Y_e[index_],                                                  \
        temperature[index_],                                          \
        &press[index_],                                               \
        &eps[index_],                                                 \
        &entropy[index_]                                              \
    );

#define GHL_TEST_LOG_WRITE_METRIC_AND_PRIMS(index_) \
    GHL_TEST_LOG_METRIC(index_);                    \
    GHL_TEST_LOG_PRIMS(index_);

// End ET-specific stuff

#endif // GRHAYLIB_TEST_H
