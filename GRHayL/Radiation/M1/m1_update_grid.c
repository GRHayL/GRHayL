/* 
 * m1_update_grid.c
 * Grid-level M1 driver. This is a temporary device for testing the M1 update on a grid.
 * For integration into Einstein Toolkit, this will be replaced by gridfunctions and scheduled functions.
 * Boundary conditions for ghost points and metric derivatives are written with LLM, proceed with caution...
 */

#include "ghl.h"
#include "ghl_m1.h"
#include "ghl_radiation.h"

#include <omp.h>
#include <stdlib.h>

#define NGROUPS 1
#define NSPECIES 1

#define NSTENCIL_CONS 5
#define NSTENCIL_FLUX 3
#define NCONS_VAR_SIZE 5

#define STENCIL_CONS_SIZE (NGROUPS * NSPECIES * NSTENCIL_CONS)
#define M1_GHOST_WIDTH 2
#define M1_GRID_INDEX(ig, ix, iy, iz, Nx, Ny, Nz) (((ig) * (Nx) * (Ny) * (Nz)) + ((ix) * (Ny) * (Nz)) + ((iy) * (Nz)) + (iz))
 
enum {
  M1_BC_PERIODIC = 0,
  M1_BC_REFLECTIVE = 1
};

#ifndef M1_GHOST_BOUNDARY_MODE
#define M1_GHOST_BOUNDARY_MODE M1_BC_REFLECTIVE
#endif

static inline int m1_ghost_index(
      const int ig,
      const int ix,
      const int iy,
      const int iz,
      const int nx,
      const int ny,
      const int nz) {
  return M1_GRID_INDEX(ig,
                       ix + M1_GHOST_WIDTH,
                       iy + M1_GHOST_WIDTH,
                       iz + M1_GHOST_WIDTH,
                       nx + 2 * M1_GHOST_WIDTH,
                       ny + 2 * M1_GHOST_WIDTH,
                       nz + 2 * M1_GHOST_WIDTH);
}

static int m1_ghost_coord_periodic(const int coord, const int npts) {
  int wrapped = coord % npts;
  if(wrapped < 0) wrapped += npts;
  return wrapped;
}

static int m1_ghost_coord_reflective(const int coord, const int npts) {
  if(coord < 0) return -coord - 1;
  if(coord >= npts) return 2 * npts - coord - 2;
  return coord;
}

static int m1_ghost_coord(const int coord, const int npts) {
#if M1_GHOST_BOUNDARY_MODE == M1_BC_REFLECTIVE
  return m1_ghost_coord_reflective(coord, npts);
#else
  return m1_ghost_coord_periodic(coord, npts);
#endif
}



/*
 * The intended flow is:
 *   1. calc_fluxes
 *   2. calc_grsources from compute_radiation_flux
 *   3. calc_update
 */
int ghl_m1_update_grid(
      const int grid_Nx,
      const int grid_Ny,
      const int grid_Nz,
      const int imin,
      const int imax,
      const int jmin,
      const int jmax,
      const int kmin,
      const int kmax,
      const double timestep,
      const ghl_m1_thc_params *restrict thc_params,
      const double delta[3],
      const ghl_metric_quantities *restrict metric_grid,
      const ghl_extrinsic_curvature *restrict curv_grid,
      const ghl_ADM_aux_quantities *restrict adm_aux_grid,
      const ghl_primitive_quantities *restrict prims_grid,
      // const double *restrict rE_grid,
      // const ghl_radiation_flux_vector *restrict rF4_grid,
      // const double *restrict rN_grid,
      const double *restrict abs_0_grid,
      const double *restrict abs_1_grid,
      const double *restrict scat_0_grid,
      const double *restrict eta_0_grid,
      const double *restrict eta_1_grid,
      
      const double *restrict scat_1_grid,
      ghl_m1_conservative_quantities *restrict rU_grid,
      ghl_m1_conservative_quantities *restrict rU_p_grid) {

  const int npoints = grid_Nx * grid_Ny * grid_Nz;
  const int ghost_nx = grid_Nx + 2 * M1_GHOST_WIDTH;
  const int ghost_ny = grid_Ny + 2 * M1_GHOST_WIDTH;
  const int ghost_nz = grid_Nz + 2 * M1_GHOST_WIDTH;
  const int ghost_npoints = ghost_nx * ghost_ny * ghost_nz;
  
  // we need to perform closure first
  ghl_radiation_pressure_tensor *rP4_grid = calloc((size_t) npoints, sizeof(*rP4_grid));

  
  ghl_m1_closure_t closure = Minerbo;
#pragma omp parallel for
  for (int ig = 0; ig < NGROUPS * NSPECIES; ++ig) {
    for (int ix = 0; ix < grid_Nx; ++ix) {
      for (int iy = 0; iy < grid_Ny; ++iy) {
        for (int iz = 0; iz < grid_Nz; ++iz) {
          m1_root_params closure_params = {0};
          int idx = M1_GRID_INDEX(ig, ix, iy, iz, grid_Nx, grid_Ny, grid_Nz);
          for (int a = 0; a < 4; a++){
            for (int b = 0; b < 4; b++){
              rP4_grid[idx].DD[a][b] = 0.0;
            }
          }
          closure_params.metric   = &metric_grid[idx];
          closure_params.adm_aux  = &adm_aux_grid[idx];
          closure_params.prims    = &prims_grid[idx];
          closure_params.E        = *rU_p_grid[idx].rE;
          closure_params.F4       = &rU_p_grid[idx].rF4;
          closure_params.P4       = &rP4_grid[idx];
          ghl_radiation_rootSolve_closure(&closure_params);
        }
      }
    }
  }





  ghl_m1_conservative_quantities *rhs_flux_grid = calloc((size_t) npoints, sizeof(*rhs_flux_grid));
  ghl_m1_conservative_quantities *rhs_gr_grid = calloc((size_t) npoints, sizeof(*rhs_gr_grid));

  double *rE_ghost = calloc((size_t) ghost_npoints, sizeof(*rE_ghost));
  ghl_radiation_flux_vector *rF4_ghost = calloc((size_t) ghost_npoints, sizeof(*rF4_ghost));
  ghl_radiation_pressure_tensor *rP4_ghost = calloc((size_t) ghost_npoints, sizeof(*rP4_ghost));
  double *rN_ghost = calloc((size_t) ghost_npoints, sizeof(*rN_ghost));
  ghl_metric_quantities *metric_ghost = calloc((size_t) ghost_npoints, sizeof(*metric_ghost));
  ghl_ADM_aux_quantities *adm_aux_ghost = calloc((size_t) ghost_npoints, sizeof(*adm_aux_ghost));
  ghl_primitive_quantities *prims_ghost = calloc((size_t) ghost_npoints, sizeof(*prims_ghost));
  double *abs_1_ghost = calloc((size_t) ghost_npoints, sizeof(*abs_1_ghost));
  double *scat_1_ghost = calloc((size_t) ghost_npoints, sizeof(*scat_1_ghost));
  for(int ig = 0; ig < NGROUPS * NSPECIES; ++ig) {
    for(int ix = 0; ix < grid_Nx; ++ix) {
      for(int iy = 0; iy < grid_Ny; ++iy) {
        for(int iz = 0; iz < grid_Nz; ++iz) {
          const int index = M1_GRID_INDEX(ig, ix, iy, iz, grid_Nx, grid_Ny, grid_Nz);
          const int ghost_index = m1_ghost_index(ig, ix, iy, iz, grid_Nx, grid_Ny, grid_Nz);
          rE_ghost[ghost_index] = *rU_p_grid[index].rE;
          rF4_ghost[ghost_index] = *rU_p_grid[index].rF4;
          rP4_ghost[ghost_index] = rP4_grid[index];
          rN_ghost[ghost_index] = *rU_p_grid[index].rN;
          metric_ghost[ghost_index] = metric_grid[index];
          adm_aux_ghost[ghost_index] = adm_aux_grid[index];
          prims_ghost[ghost_index] = prims_grid[index];
          abs_1_ghost[ghost_index] = abs_1_grid[index];
          scat_1_ghost[ghost_index] = scat_1_grid[index];
        }//iz
      }//iy
    }//ix
  }//ig
#pragma omp parallel for
  for(int ig = 0; ig < NGROUPS * NSPECIES; ++ig) {
    for(int iz = -M1_GHOST_WIDTH; iz < grid_Nz + M1_GHOST_WIDTH; ++iz) {
      for(int iy = -M1_GHOST_WIDTH; iy < grid_Ny + M1_GHOST_WIDTH; ++iy) {
        for(int ix = -M1_GHOST_WIDTH; ix < grid_Nx + M1_GHOST_WIDTH; ++ix) {
          const int in_domain_x = (ix >= 0 && ix < grid_Nx);
          const int in_domain_y = (iy >= 0 && iy < grid_Ny);
          const int in_domain_z = (iz >= 0 && iz < grid_Nz);
          if(in_domain_x && in_domain_y && in_domain_z) continue;

          const int ghost_ix = ix + M1_GHOST_WIDTH;
          const int ghost_iy = iy + M1_GHOST_WIDTH;
          const int ghost_iz = iz + M1_GHOST_WIDTH;
          const int ghost_index = M1_GRID_INDEX(ig, ghost_ix, ghost_iy, ghost_iz, ghost_nx, ghost_ny, ghost_nz);

          const int src_ix = m1_ghost_coord(ix, grid_Nx);
          const int src_iy = m1_ghost_coord(iy, grid_Ny);
          const int src_iz = m1_ghost_coord(iz, grid_Nz);
          const int src_index = M1_GRID_INDEX(ig, src_ix, src_iy, src_iz, grid_Nx, grid_Ny, grid_Nz);

          rE_ghost[ghost_index] = *rU_p_grid[src_index].rE;
          rF4_ghost[ghost_index] = *rU_p_grid[src_index].rF4;
          rP4_ghost[ghost_index] = rP4_grid[src_index];
          rN_ghost[ghost_index] = *rU_p_grid[src_index].rN;
          metric_ghost[ghost_index] = metric_grid[src_index];
          adm_aux_ghost[ghost_index] = adm_aux_grid[src_index];
          prims_ghost[ghost_index] = prims_grid[src_index];
          abs_1_ghost[ghost_index] = abs_1_grid[src_index];
          scat_1_ghost[ghost_index] = scat_1_grid[src_index];
        }//ix
      }//iy
    }//iz
  }//ig

  double rE_rhs[NGROUPS * NSPECIES] = {0.0};
  double rFx_rhs[NGROUPS * NSPECIES] = {0.0};
  double rFy_rhs[NGROUPS * NSPECIES] = {0.0};
  double rFz_rhs[NGROUPS * NSPECIES] = {0.0};
  double rN_rhs[NGROUPS * NSPECIES] = {0.0};

  /* Stage 1: centered finite-difference flux stencils on the ghosted grid. */
  // printf("m1_update_grid Stage 1: compute the directional flux RHS for this grid point.\n");
#pragma omp parallel for
  for(int ig = 0; ig < NGROUPS * NSPECIES; ++ig) {
    for(int iz = kmin; iz < kmax + 1; ++iz) {
      for(int iy = jmin; iy < jmax + 1; ++iy) {
        for(int ix = imin; ix < imax + 1; ++ix) {
          const int index = M1_GRID_INDEX(ig, ix, iy, iz, grid_Nx, grid_Ny, grid_Nz);
          double stencil_rE[STENCIL_CONS_SIZE] = {0.0};
          ghl_radiation_flux_vector stencil_rF4[STENCIL_CONS_SIZE] = {0};
          ghl_radiation_pressure_tensor stencil_rP4[STENCIL_CONS_SIZE] = {0};
          double stencil_rN[STENCIL_CONS_SIZE] = {0.0};
          ghl_metric_quantities stencil_metric[NSTENCIL_FLUX] = {0};
          ghl_ADM_aux_quantities stencil_adm_aux[NSTENCIL_FLUX] = {0};
          ghl_primitive_quantities stencil_prims[NSTENCIL_FLUX] = {0};
          double stencil_abs_1[NSTENCIL_FLUX] = {0.0};
          double stencil_scat_1[NSTENCIL_FLUX] = {0.0};
          

          ghl_m1_conservative_quantities rhs_flux = {0};

          for(int stage_dir = 0; stage_dir < 3; ++stage_dir) {
            for(int k_stencil = 0; k_stencil < NSTENCIL_CONS; ++k_stencil) {
              const int stencil_index = ig * NSTENCIL_CONS + k_stencil;
              int stencil_ix = ix + (stage_dir == 0 ? (k_stencil - 2) : 0);
              int stencil_iy = iy + (stage_dir == 1 ? (k_stencil - 2) : 0);
              int stencil_iz = iz + (stage_dir == 2 ? (k_stencil - 2) : 0);
              const int point_index = m1_ghost_index(ig, stencil_ix, stencil_iy, stencil_iz,
                                                     grid_Nx, grid_Ny, grid_Nz);

              stencil_rE[stencil_index]  = rE_ghost[point_index];
              stencil_rF4[stencil_index] = rF4_ghost[point_index];
              stencil_rP4[stencil_index] = rP4_ghost[point_index];
              stencil_rN[stencil_index]  = rN_ghost[point_index];
            }
            // This might not work if NSPECIES > 1
            for(int k_flux = 0; k_flux < NSTENCIL_FLUX; ++k_flux) {
              const int offset = k_flux - 1;
              int stencil_ix = ix + (stage_dir == 0 ? offset : 0);
              int stencil_iy = iy + (stage_dir == 1 ? offset : 0);
              int stencil_iz = iz + (stage_dir == 2 ? offset : 0);
              const int point_index = m1_ghost_index(ig, stencil_ix, stencil_iy, stencil_iz,
                                                     grid_Nx, grid_Ny, grid_Nz);

              stencil_metric[k_flux]    = metric_ghost[point_index];
              stencil_adm_aux[k_flux]   = adm_aux_ghost[point_index];
              stencil_prims[k_flux]     = prims_ghost[point_index];
              stencil_abs_1[k_flux]     = abs_1_ghost[point_index];
              stencil_scat_1[k_flux]    = scat_1_ghost[point_index];
            }

            
            // printf("stencil Fx = %f, %f, %f, %f, %f\n", stencil_rF4[ig * NSTENCIL_CONS + 0].D[1], stencil_rF4[ig * NSTENCIL_CONS + 1].D[1], stencil_rF4[ig * NSTENCIL_CONS + 2].D[1], stencil_rF4[ig * NSTENCIL_CONS + 3].D[1], stencil_rF4[ig * NSTENCIL_CONS + 4].D[1]);
            ghl_m1_calc_flux(
                  timestep,
                  stencil_metric,
                  stencil_adm_aux,
                  stencil_prims,
                  thc_params,
                  stencil_rE,
                  stencil_rF4,
                  stencil_rP4,
                  stencil_rN,
                  delta,
                  stencil_abs_1,
                  stencil_scat_1,
                  stage_dir,
                  &rhs_flux);
          
                 // Total flux needs to be summed over all directions
            rhs_flux_grid[index].rN[ig] += rhs_flux.rN[ig];
            rhs_flux_grid[index].rF4[ig].D[1] += rhs_flux.rF4[ig].D[1];
            rhs_flux_grid[index].rF4[ig].D[2] += rhs_flux.rF4[ig].D[2];
            rhs_flux_grid[index].rF4[ig].D[3] += rhs_flux.rF4[ig].D[3];
            rhs_flux_grid[index].rE[ig] += rhs_flux.rE[ig];
          } // stage_dir
        }//ix
      }//iy
    }//iz
  }//ig

  free(rE_ghost);
  free(rF4_ghost);
  free(rP4_ghost);
  free(rN_ghost);
  free(metric_ghost);
  free(adm_aux_ghost);
  free(prims_ghost);
  free(abs_1_ghost);
  free(scat_1_ghost);

  /* Stage 2 ... */


  /* Stage 2: compute GR source terms directly from metric derivatives on the grid. */
#pragma omp parallel for
  for(int ig = 0; ig < NGROUPS * NSPECIES; ++ig) {
    for(int iz = kmin; iz < kmax + 1; ++iz) {
      for(int iy = jmin; iy < jmax + 1; ++iy) {
        for(int ix = imin; ix < imax + 1; ++ix) {
          const int index = M1_GRID_INDEX(ig, ix, iy, iz, grid_Nx, grid_Ny, grid_Nz);
          ghl_metric_quantities metric_derivs_x = {0};
          ghl_metric_quantities metric_derivs_y = {0};
          ghl_metric_quantities metric_derivs_z = {0};
          const ghl_extrinsic_curvature curv = curv_grid[index];

          const int idx0 = M1_GRID_INDEX(ig, ix, iy, iz, grid_Nx, grid_Ny, grid_Nz);

          const int idx_xm2 = (ix >= 2) ? M1_GRID_INDEX(ig, ix - 2, iy, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_xm1 = (ix >= 1) ? M1_GRID_INDEX(ig, ix - 1, iy, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_xp1 = (ix + 1 < grid_Nx) ? M1_GRID_INDEX(ig, ix + 1, iy, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_xp2 = (ix + 2 < grid_Nx) ? M1_GRID_INDEX(ig, ix + 2, iy, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_ym2 = (iy >= 2) ? M1_GRID_INDEX(ig, ix, iy - 2, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_ym1 = (iy >= 1) ? M1_GRID_INDEX(ig, ix, iy - 1, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_yp1 = (iy + 1 < grid_Ny) ? M1_GRID_INDEX(ig, ix, iy + 1, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_yp2 = (iy + 2 < grid_Ny) ? M1_GRID_INDEX(ig, ix, iy + 2, iz, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_zm2 = (iz >= 2) ? M1_GRID_INDEX(ig, ix, iy, iz - 2, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_zm1 = (iz >= 1) ? M1_GRID_INDEX(ig, ix, iy, iz - 1, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_zp1 = (iz + 1 < grid_Nz) ? M1_GRID_INDEX(ig, ix, iy, iz + 1, grid_Nx, grid_Ny, grid_Nz) : idx0;
          const int idx_zp2 = (iz + 2 < grid_Nz) ? M1_GRID_INDEX(ig, ix, iy, iz + 2, grid_Nx, grid_Ny, grid_Nz) : idx0;

          const double dfdx_c = (ix == 0)
              ? (-3.0 * metric_grid[idx0].lapse + 4.0 * metric_grid[idx_xp1].lapse - metric_grid[idx_xp2].lapse) / (2.0 * delta[0])
              : (ix == grid_Nx - 1)
                  ? (3.0 * metric_grid[idx0].lapse - 4.0 * metric_grid[idx_xm1].lapse + metric_grid[idx_xm2].lapse) / (2.0 * delta[0])
                  : (metric_grid[idx_xp1].lapse - metric_grid[idx_xm1].lapse) / (2.0 * delta[0]);
          const double dfdy_c = (iy == 0)
              ? (-3.0 * metric_grid[idx0].lapse + 4.0 * metric_grid[idx_yp1].lapse - metric_grid[idx_yp2].lapse) / (2.0 * delta[1])
              : (iy == grid_Ny - 1)
                  ? (3.0 * metric_grid[idx0].lapse - 4.0 * metric_grid[idx_ym1].lapse + metric_grid[idx_ym2].lapse) / (2.0 * delta[1])
                  : (metric_grid[idx_yp1].lapse - metric_grid[idx_ym1].lapse) / (2.0 * delta[1]);
          const double dfdz_c = (iz == 0)
              ? (-3.0 * metric_grid[idx0].lapse + 4.0 * metric_grid[idx_zp1].lapse - metric_grid[idx_zp2].lapse) / (2.0 * delta[2])
              : (iz == grid_Nz - 1)
                  ? (3.0 * metric_grid[idx0].lapse - 4.0 * metric_grid[idx_zm1].lapse + metric_grid[idx_zm2].lapse) / (2.0 * delta[2])
                  : (metric_grid[idx_zp1].lapse - metric_grid[idx_zm1].lapse) / (2.0 * delta[2]);

          metric_derivs_x.lapse = dfdx_c;
          metric_derivs_y.lapse = dfdy_c;
          metric_derivs_z.lapse = dfdz_c;

          for(int a = 0; a < 3; ++a) {
            metric_derivs_x.betaU[a] = (ix == 0)
                ? (-3.0 * metric_grid[idx0].betaU[a] + 4.0 * metric_grid[idx_xp1].betaU[a] - metric_grid[idx_xp2].betaU[a]) / (2.0 * delta[0])
                : (ix == grid_Nx - 1)
                    ? (3.0 * metric_grid[idx0].betaU[a] - 4.0 * metric_grid[idx_xm1].betaU[a] + metric_grid[idx_xm2].betaU[a]) / (2.0 * delta[0])
                    : (metric_grid[idx_xp1].betaU[a] - metric_grid[idx_xm1].betaU[a]) / (2.0 * delta[0]);
            metric_derivs_y.betaU[a] = (iy == 0)
                ? (-3.0 * metric_grid[idx0].betaU[a] + 4.0 * metric_grid[idx_yp1].betaU[a] - metric_grid[idx_yp2].betaU[a]) / (2.0 * delta[1])
                : (iy == grid_Ny - 1)
                    ? (3.0 * metric_grid[idx0].betaU[a] - 4.0 * metric_grid[idx_ym1].betaU[a] + metric_grid[idx_ym2].betaU[a]) / (2.0 * delta[1])
                    : (metric_grid[idx_yp1].betaU[a] - metric_grid[idx_ym1].betaU[a]) / (2.0 * delta[1]);
            metric_derivs_z.betaU[a] = (iz == 0)
                ? (-3.0 * metric_grid[idx0].betaU[a] + 4.0 * metric_grid[idx_zp1].betaU[a] - metric_grid[idx_zp2].betaU[a]) / (2.0 * delta[2])
                : (iz == grid_Nz - 1)
                    ? (3.0 * metric_grid[idx0].betaU[a] - 4.0 * metric_grid[idx_zm1].betaU[a] + metric_grid[idx_zm2].betaU[a]) / (2.0 * delta[2])
                    : (metric_grid[idx_zp1].betaU[a] - metric_grid[idx_zm1].betaU[a]) / (2.0 * delta[2]);
          }

          for(int a = 0; a < 3; ++a) {
            for(int b = 0; b < 3; ++b) {
              metric_derivs_x.gammaDD[a][b] = (ix == 0)
                  ? (-3.0 * metric_grid[idx0].gammaDD[a][b] + 4.0 * metric_grid[idx_xp1].gammaDD[a][b] - metric_grid[idx_xp2].gammaDD[a][b]) / (2.0 * delta[0])
                  : (ix == grid_Nx - 1)
                      ? (3.0 * metric_grid[idx0].gammaDD[a][b] - 4.0 * metric_grid[idx_xm1].gammaDD[a][b] + metric_grid[idx_xm2].gammaDD[a][b]) / (2.0 * delta[0])
                      : (metric_grid[idx_xp1].gammaDD[a][b] - metric_grid[idx_xm1].gammaDD[a][b]) / (2.0 * delta[0]);
              metric_derivs_y.gammaDD[a][b] = (iy == 0)
                  ? (-3.0 * metric_grid[idx0].gammaDD[a][b] + 4.0 * metric_grid[idx_yp1].gammaDD[a][b] - metric_grid[idx_yp2].gammaDD[a][b]) / (2.0 * delta[1])
                  : (iy == grid_Ny - 1)
                      ? (3.0 * metric_grid[idx0].gammaDD[a][b] - 4.0 * metric_grid[idx_ym1].gammaDD[a][b] + metric_grid[idx_ym2].gammaDD[a][b]) / (2.0 * delta[1])
                      : (metric_grid[idx_yp1].gammaDD[a][b] - metric_grid[idx_ym1].gammaDD[a][b]) / (2.0 * delta[1]);
              metric_derivs_z.gammaDD[a][b] = (iz == 0)
                  ? (-3.0 * metric_grid[idx0].gammaDD[a][b] + 4.0 * metric_grid[idx_zp1].gammaDD[a][b] - metric_grid[idx_zp2].gammaDD[a][b]) / (2.0 * delta[2])
                  : (iz == grid_Nz - 1)
                      ? (3.0 * metric_grid[idx0].gammaDD[a][b] - 4.0 * metric_grid[idx_zm1].gammaDD[a][b] + metric_grid[idx_zm2].gammaDD[a][b]) / (2.0 * delta[2])
                      : (metric_grid[idx_zp1].gammaDD[a][b] - metric_grid[idx_zm1].gammaDD[a][b]) / (2.0 * delta[2]);
            }
          }
          int gr_err = ghl_m1_calc_grsources(
                timestep,
                &metric_grid[index],
                &metric_derivs_x,
                &metric_derivs_y,
                &metric_derivs_z,
                *rU_p_grid[index].rE,
                &rU_p_grid[index].rF4,
                &rP4_grid[index],
                &curv_grid[index],
                &rhs_gr_grid[index]);
        }//ix
      }//iy
    }//iz
  }//ig
  

  /* Stage 3: pack RHS and advance the local state. */
  // printf("m1_update_grid Stage 3: pack RHS and advance the local state.\n");
#pragma omp parallel for
  ghl_m1_conservative_quantities *rhs_grid = calloc((size_t) npoints, sizeof(*rhs_grid));
  for(int ig = 0; ig < NGROUPS * NSPECIES; ++ig) {
    double abs_0_arr[NGROUPS * NSPECIES] = {0.0};
    double abs_1_arr[NGROUPS * NSPECIES] = {0.0};
    double eta_0_arr[NGROUPS * NSPECIES] = {0.0};
    double eta_1_arr[NGROUPS * NSPECIES] = {0.0};
    double scat_0_arr[NGROUPS * NSPECIES] = {0.0};
    double scat_1_arr[NGROUPS * NSPECIES] = {0.0};

    for(int iz = kmin; iz < kmax + 1; ++iz) {
      for(int iy = jmin; iy < jmax + 1; ++iy) {
        for(int ix = imin; ix < imax + 1; ++ix) {
          const int index = M1_GRID_INDEX(ig, ix, iy, iz, grid_Nx, grid_Ny, grid_Nz);
          
          // this is not the correct way to do all species, modify when handling multiple species
          abs_0_arr[ig]  = abs_0_grid[index];
          abs_1_arr[ig]  = abs_1_grid[index];
          scat_0_arr[ig] = scat_0_grid[index];
          scat_1_arr[ig] = scat_1_grid[index];
          eta_0_arr[ig]  = eta_0_grid[index];
          eta_1_arr[ig]  = eta_1_grid[index];

          ghl_m1_update(
                timestep,
                &metric_grid[index],
                &curv_grid[index],
                &adm_aux_grid[index],
                &prims_grid[index],
                thc_params,
                abs_0_arr,
                abs_1_arr,
                eta_0_arr,
                eta_1_arr,
                scat_0_arr, //scat_0 neutrinos number is conserved, so this is never used
                scat_1_arr,
                &rU_grid[index],
                &rU_p_grid[index],
                &rhs_flux_grid[index],
                &rhs_gr_grid[index]);
        }//ix
      }//iy
    }//iz
  }//ig

  free(rhs_flux_grid);
  free(rhs_gr_grid);

  return 0;
}