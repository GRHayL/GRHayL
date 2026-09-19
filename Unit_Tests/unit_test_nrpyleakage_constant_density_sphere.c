#include <float.h>

#include "ghl_unit_tests.h"

static inline bool sphere_pert_test_fail(
      const double trusted,
      const double computed,
      const double perturbed) {
  if(trusted == 0.0 && perturbed == 0.0) {
    return computed != 0.0;
  }
  return ghl_pert_test_fail_with_tolerance(
        trusted, computed, perturbed, 1024.0 * DBL_EPSILON, 0.0);
}

#define IDX3D(i0,i1,i2) ( (i0) + Nt0*( (i1) + Nt1*(i2) ) )

static void set_opacity_struct_from_gfs(
      const int index,
      const double *restrict kappa_0_nue,
      const double *restrict kappa_1_nue,
      const double *restrict kappa_0_anue,
      const double *restrict kappa_1_anue,
      const double *restrict kappa_0_nux,
      const double *restrict kappa_1_nux,
      ghl_neutrino_opacities *restrict kappa) {

  kappa->nue [0] = kappa_0_nue [index];
  kappa->nue [1] = kappa_1_nue [index];
  kappa->anue[0] = kappa_0_anue[index];
  kappa->anue[1] = kappa_1_anue[index];
  kappa->nux [0] = kappa_0_nux [index];
  kappa->nux [1] = kappa_1_nux [index];
}

static void set_optical_depths_struct_from_gfs(
      const int index,
      const double *restrict tau_0_nue,
      const double *restrict tau_1_nue,
      const double *restrict tau_0_anue,
      const double *restrict tau_1_anue,
      const double *restrict tau_0_nux,
      const double *restrict tau_1_nux,
      ghl_neutrino_optical_depths *restrict tau) {

  tau->nue [0] = tau_0_nue [index];
  tau->nue [1] = tau_1_nue [index];
  tau->anue[0] = tau_0_anue[index];
  tau->anue[1] = tau_1_anue[index];
  tau->nux [0] = tau_0_nux [index];
  tau->nux [1] = tau_1_nux [index];
}

/** Gather the sphere caller's stencil and update one cell. */
static void update_optical_depth_cell(
      const int i0,
      const int i1,
      const int i2,
      const int Nt0,
      const int Nt1,
      const double dxx[restrict static 3],
      const double stencil_gxx[restrict static 3],
      const double stencil_gyy[restrict static 3],
      const double stencil_gzz[restrict static 3],
      double *restrict const kappa_nue[restrict static 2],
      double *restrict const kappa_anue[restrict static 2],
      double *restrict const kappa_nux[restrict static 2],
      double *restrict const tau_nue[restrict static 2],
      double *restrict const tau_anue[restrict static 2],
      double *restrict const tau_nux[restrict static 2],
      ghl_neutrino_optical_depths *restrict tau_i_j_k) {
  const int i_j_k = IDX3D(i0, i1, i2);
  const int im1_j_k = IDX3D(i0 - 1, i1, i2);
  const int ip1_j_k = IDX3D(i0 + 1, i1, i2);
  const int i_jm1_k = IDX3D(i0, i1 - 1, i2);
  const int i_jp1_k = IDX3D(i0, i1 + 1, i2);
  const int i_j_km1 = IDX3D(i0, i1, i2 - 1);
  const int i_j_kp1 = IDX3D(i0, i1, i2 + 1);

  ghl_neutrino_opacities kappa_center, kappa_im1, kappa_ip1;
  ghl_neutrino_opacities kappa_jm1, kappa_jp1, kappa_km1, kappa_kp1;
  ghl_neutrino_optical_depths tau_im1, tau_ip1, tau_jm1, tau_jp1;
  ghl_neutrino_optical_depths tau_km1, tau_kp1;
#define LOAD_KAPPA(index, value)                                                       \
  set_opacity_struct_from_gfs(                                                         \
        index, kappa_nue[0], kappa_nue[1], kappa_anue[0], kappa_anue[1], kappa_nux[0], \
        kappa_nux[1], value)
#define LOAD_TAU(index, value)                                               \
  set_optical_depths_struct_from_gfs(                                        \
        index, tau_nue[0], tau_nue[1], tau_anue[0], tau_anue[1], tau_nux[0], \
        tau_nux[1], value)
  LOAD_KAPPA(i_j_k, &kappa_center);
  LOAD_KAPPA(im1_j_k, &kappa_im1);
  LOAD_KAPPA(ip1_j_k, &kappa_ip1);
  LOAD_KAPPA(i_jm1_k, &kappa_jm1);
  LOAD_KAPPA(i_jp1_k, &kappa_jp1);
  LOAD_KAPPA(i_j_km1, &kappa_km1);
  LOAD_KAPPA(i_j_kp1, &kappa_kp1);
  LOAD_TAU(im1_j_k, &tau_im1);
  LOAD_TAU(ip1_j_k, &tau_ip1);
  LOAD_TAU(i_jm1_k, &tau_jm1);
  LOAD_TAU(i_jp1_k, &tau_jp1);
  LOAD_TAU(i_j_km1, &tau_km1);
  LOAD_TAU(i_j_kp1, &tau_kp1);
#undef LOAD_KAPPA
#undef LOAD_TAU

  NRPyLeakage_optical_depths_PathOfLeastResistance(
        dxx, stencil_gxx, stencil_gyy, stencil_gzz, &kappa_im1, &kappa_ip1, &kappa_jm1,
        &kappa_jp1, &kappa_km1, &kappa_kp1, &tau_im1, &tau_ip1, &tau_jm1, &tau_jp1,
        &tau_km1, &tau_kp1, &kappa_center, tau_i_j_k);
}

/** Check the sphere gather/call path with an order-sensitive stencil. */
static void check_sphere_caller_neighbor_order(void) {
  enum { N = 27 };
  double kappa_storage[3][2][N];
  double tau_storage[3][2][N];
  double *kappa[3][2];
  double *tau[3][2];
  for(int species = 0; species < 3; species++) {
    for(int moment = 0; moment < 2; moment++) {
      kappa[species][moment] = kappa_storage[species][moment];
      tau[species][moment] = tau_storage[species][moment];
      for(int index = 0; index < N; index++) {
        kappa[species][moment][index] = 2.0;
        tau[species][moment][index] = 1000.0;
      }
    }
  }
  const int Nt0 = 3, Nt1 = 3;
  const int center = IDX3D(1, 1, 1);
  const int xminus = IDX3D(0, 1, 1);
  const int xplus = IDX3D(2, 1, 1);
  for(int species = 0; species < 3; species++) {
    for(int moment = 0; moment < 2; moment++) {
      kappa[species][moment][center] = 2.0;
      kappa[species][moment][xminus] = 4.0;
      kappa[species][moment][xplus] = 20.0;
      tau[species][moment][xminus] = 0.0;
      tau[species][moment][xplus] = 100.0;
    }
  }
  const double dxx[3] = { 1.0, 1.0, 1.0 };
  const double gxx[3] = { 1.0, 4.0, 25.0 };
  const double gyy[3] = { 1.0, 1.0, 1.0 };
  const double gzz[3] = { 1.0, 1.0, 1.0 };
  ghl_neutrino_optical_depths result;
  update_optical_depth_cell(
        1, 1, 1, Nt0, Nt1, dxx, gxx, gyy, gzz, kappa[0], kappa[1], kappa[2], tau[0],
        tau[1], tau[2], &result);
  const double expected = 3.0 * sqrt(2.5);
  const double values[] = {
    result.nue[0],  result.nue[1], result.anue[0],
    result.anue[1], result.nux[0], result.nux[1],
  };
  for(size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
    if(fabs(values[i] - expected) > 32.0 * DBL_EPSILON * expected) {
      ghl_error("Sphere caller reversed neighbor %zu\n", i);
    }
  }
}

void constantdensitysphere_test(
      const ghl_eos_parameters *restrict eos,
      const int test_key) {

  // Step 1: Set basic parameters
  const int N0     = 32;
  const int N1     = N0;
  const int N2     = N0;
  const int Ng0    =  2;
  const int Ng1    = Ng0;
  const int Ng2    = Ng0;
  const int Nt0    = N0 + 2*Ng0;
  const int Nt1    = N1 + 2*Ng1;
  const int Nt2    = N2 + 2*Ng2;
  const int Ntotal = Nt0 * Nt1 * Nt2;
  const double xmax  = +2;
  const double xmin  = -2;
  const double ymax  = xmax;
  const double ymin  = xmin;
  const double zmax  = xmax;
  const double zmin  = xmin;
  const double dx    = (xmax-xmin)/((double)N0);
  const double dy    = (xmax-xmin)/((double)N1);
  const double dz    = (xmax-xmin)/((double)N2);
  const double rSph  = 1;

  // Step 2: Allocate memory for the metric, opacities, and optical depths
  double **kappa_nue  = (double **)malloc(sizeof(double *)*2);
  double **kappa_anue = (double **)malloc(sizeof(double *)*2);
  double **kappa_nux  = (double **)malloc(sizeof(double *)*2);
  double **tau_nue    = (double **)malloc(sizeof(double *)*2);
  double **tau_anue   = (double **)malloc(sizeof(double *)*2);
  double **tau_nux    = (double **)malloc(sizeof(double *)*2);
  double **tau_nue_p  = (double **)malloc(sizeof(double *)*2);
  double **tau_anue_p = (double **)malloc(sizeof(double *)*2);
  double **tau_nux_p  = (double **)malloc(sizeof(double *)*2);
  for(int i=0;i<2;i++) {
   kappa_nue [i] = (double *)malloc(sizeof(double)*Ntotal);
   kappa_anue[i] = (double *)malloc(sizeof(double)*Ntotal);
   kappa_nux [i] = (double *)malloc(sizeof(double)*Ntotal);
   tau_nue   [i] = (double *)malloc(sizeof(double)*Ntotal);
   tau_anue  [i] = (double *)malloc(sizeof(double)*Ntotal);
   tau_nux   [i] = (double *)malloc(sizeof(double)*Ntotal);
   tau_nue_p [i] = (double *)malloc(sizeof(double)*Ntotal);
   tau_anue_p[i] = (double *)malloc(sizeof(double)*Ntotal);
   tau_nux_p [i] = (double *)malloc(sizeof(double)*Ntotal);
  }

  // Step 3: Define hydro quantities at the sphere interior and exterior
  //         The parameters below are extracted from https://arxiv.org/abs/2106.05356
  // Step 3.a: Interior of the sphere (sec. 4.2.1 of above referece)
  double rho_interior = 9.8e13 * NRPyLeakage_units_cgs_to_geom_D;
  double Y_e_interior = 0.1;
  double T_interior   = 8.0;
  // Step 3.b: Exterior of the sphere (sec. 4.2.1 of above referece)
  double rho_exterior = 6.0e7 * NRPyLeakage_units_cgs_to_geom_D;
  double Y_e_exterior = 0.5;
  const double perturbation_amplitude = 1e-14;
  const double minimum_safe_T
        = nextafter(eos->T_min / (1.0 - perturbation_amplitude), INFINITY);
  double T_exterior = fmax(0.01, minimum_safe_T);

  if( test_key == 1 ) {
    rho_interior *= (1 + randf(-1, 1) * perturbation_amplitude);
    Y_e_interior *= (1 + randf(-1, 1) * perturbation_amplitude);
    T_interior *= (1 + randf(-1, 1) * perturbation_amplitude);
    rho_exterior *= (1 + randf(-1, 1) * perturbation_amplitude);
    Y_e_exterior *= (1 + randf(-1, 1) * perturbation_amplitude);
    T_exterior *= (1 + randf(-1, 1) * perturbation_amplitude);
  }

  // Step 4: Compute opacities in the interior and exterior
  ghl_neutrino_optical_depths tau_in = {{0,0},{0,0},{0,0}};

  // Step 4.a: Interior
  ghl_neutrino_opacities kappa_interior;
  ghl_error_codes_t error = NRPyLeakage_compute_neutrino_opacities(eos,rho_interior,Y_e_interior,T_interior,
                                                                   &tau_in, &kappa_interior);
  ghl_abort_if_error(error);

  // Step 4.b: Exterior
  ghl_neutrino_opacities kappa_exterior;
  error = NRPyLeakage_compute_neutrino_opacities(eos,rho_exterior,Y_e_exterior,T_exterior,
                                                &tau_in, &kappa_exterior);
  ghl_abort_if_error(error);

  // Step 5: Print basic information
  ghl_info("Test information:\n");
  ghl_info("    Domain properties:\n");
  ghl_info("        - Sphere radius = %22.15e\n",rSph);
  ghl_info("        - xmin          = %22.15e\n",xmin);
  ghl_info("        - xmax          = %22.15e\n",xmax);
  ghl_info("        - ymin          = %22.15e\n",ymin);
  ghl_info("        - ymax          = %22.15e\n",ymax);
  ghl_info("        - zmin          = %22.15e\n",zmin);
  ghl_info("        - zmax          = %22.15e\n",zmax);
  ghl_info("        - Nx            = (%d) + 2x(%d)\n",N0,Ng0);
  ghl_info("        - Ny            = (%d) + 2x(%d)\n",N1,Ng1);
  ghl_info("        - Nz            = (%d) + 2x(%d)\n",N2,Ng2);
  ghl_info("        - dx            = %22.15e\n",dx);
  ghl_info("        - dy            = %22.15e\n",dy);
  ghl_info("        - dz            = %22.15e\n",dz);
  ghl_info("    Hydro quantities at sphere interior:\n");
  ghl_info("        - rho           = %22.15e\n", rho_interior);
  ghl_info("        - Y_e           = %22.15e\n", Y_e_interior);
  ghl_info("        -  T            = %22.15e\n", T_interior);
  ghl_info("        - kappa_0_nue   = %22.15e\n", kappa_interior.nue [0]);
  ghl_info("        - kappa_1_nue   = %22.15e\n", kappa_interior.nue [1]);
  ghl_info("        - kappa_0_anue  = %22.15e\n", kappa_interior.anue[0]);
  ghl_info("        - kappa_1_anue  = %22.15e\n", kappa_interior.anue[1]);
  ghl_info("        - kappa_0_nux   = %22.15e\n", kappa_interior.nux [0]);
  ghl_info("        - kappa_1_nux   = %22.15e\n", kappa_interior.nux [1]);
  ghl_info("    Hydro quantities at sphere exterior:\n");
  ghl_info("        - rho           = %22.15e\n", rho_exterior);
  ghl_info("        - Y_e           = %22.15e\n", Y_e_exterior);
  ghl_info("        -  T            = %22.15e\n", T_exterior);
  ghl_info("        - kappa_0_nue   = %22.15e\n", kappa_exterior.nue [0]);
  ghl_info("        - kappa_1_nue   = %22.15e\n", kappa_exterior.nue [1]);
  ghl_info("        - kappa_0_anue  = %22.15e\n", kappa_exterior.anue[0]);
  ghl_info("        - kappa_1_anue  = %22.15e\n", kappa_exterior.anue[1]);
  ghl_info("        - kappa_0_nux   = %22.15e\n", kappa_exterior.nux [0]);
  ghl_info("        - kappa_1_nux   = %22.15e\n", kappa_exterior.nux [1]);

  // Step 6: Loop over the grid, set opacities
  //         and initialize optical depth to zero
#pragma omp parallel for
  for(int i2=0;i2<Nt2;i2++) {
    const double z = zmin + (i2-Ng2+0.5)*dz;
    for(int i1=0;i1<Nt1;i1++) {
      const double y = ymin + (i1-Ng1+0.5)*dy;
      for(int i0=0;i0<Nt0;i0++) {
        const double x = xmin + (i0-Ng0+0.5)*dx;

        // Step 3.a: Set local index
        const int index = IDX3D(i0,i1,i2);

        // Step 3.c: Set local opacities
        const double r = sqrt( x*x + y*y + z*z );
        if( r > rSph ) {
          // Exterior
          kappa_nue [0][index] = kappa_exterior.nue [0];
          kappa_nue [1][index] = kappa_exterior.nue [1];
          kappa_anue[0][index] = kappa_exterior.anue[0];
          kappa_anue[1][index] = kappa_exterior.anue[1];
          kappa_nux [0][index] = kappa_exterior.nux [0];
          kappa_nux [1][index] = kappa_exterior.nux [1];
        }
        else {
          // Interior
          kappa_nue [0][index] = kappa_interior.nue [0];
          kappa_nue [1][index] = kappa_interior.nue [1];
          kappa_anue[0][index] = kappa_interior.anue[0];
          kappa_anue[1][index] = kappa_interior.anue[1];
          kappa_nux [0][index] = kappa_interior.nux [0];
          kappa_nux [1][index] = kappa_interior.nux [1];
        }

        // Step 3.d: Initialize optical depths to zero
        tau_nue [0][index] = tau_nue_p [0][index] = 0.0;
        tau_nue [1][index] = tau_nue_p [1][index] = 0.0;
        tau_anue[0][index] = tau_anue_p[0][index] = 0.0;
        tau_anue[1][index] = tau_anue_p[1][index] = 0.0;
        tau_nux [0][index] = tau_nux_p [0][index] = 0.0;
        tau_nux [1][index] = tau_nux_p [1][index] = 0.0;
      }
    }
  }

  // Step 4: Now update the optical depth
  int n;
  for(n=0;n<20*N0;n++) {

    // Step 4.a: Copy optical depth to previous level
    double l2norm = 0.0;
#pragma omp parallel for reduction(+:l2norm)
    for(int i2=0;i2<Nt2;i2++) {
      for(int i1=0;i1<Nt1;i1++) {
        for(int i0=0;i0<Nt0;i0++) {
          const int index = IDX3D(i0,i1,i2);
          // Read from main memory
          const double tau_0_nue    = tau_nue   [0][index];
          const double tau_1_nue    = tau_nue   [1][index];
          const double tau_0_anue   = tau_anue  [0][index];
          const double tau_1_anue   = tau_anue  [1][index];
          const double tau_0_nux    = tau_nux   [0][index];
          const double tau_1_nux    = tau_nux   [1][index];
          const double tau_0_nue_p  = tau_nue_p [0][index];
          const double tau_1_nue_p  = tau_nue_p [1][index];
          const double tau_0_anue_p = tau_anue_p[0][index];
          const double tau_1_anue_p = tau_anue_p[1][index];
          const double tau_0_nux_p  = tau_nux_p [0][index];
          const double tau_1_nux_p  = tau_nux_p [1][index];

          // Now compute the l2 norm
          const double diff_tau_0_nue  = tau_0_nue -tau_0_nue_p;
          const double diff_tau_1_nue  = tau_1_nue -tau_1_nue_p;
          const double diff_tau_0_anue = tau_0_anue-tau_0_anue_p;
          const double diff_tau_1_anue = tau_1_anue-tau_1_anue_p;
          const double diff_tau_0_nux  = tau_0_nux -tau_0_nux_p;
          const double diff_tau_1_nux  = tau_1_nux -tau_1_nux_p;

          l2norm += diff_tau_0_nue*diff_tau_0_nue;
          l2norm += diff_tau_1_nue*diff_tau_1_nue;
          l2norm += diff_tau_0_anue*diff_tau_0_anue;
          l2norm += diff_tau_1_anue*diff_tau_1_anue;
          l2norm += diff_tau_0_nux*diff_tau_0_nux;
          l2norm += diff_tau_1_nux*diff_tau_1_nux;

          tau_nue_p [0][index] = tau_nue [0][index];
          tau_nue_p [1][index] = tau_nue [1][index];
          tau_anue_p[0][index] = tau_anue[0][index];
          tau_anue_p[1][index] = tau_anue[1][index];
          tau_nux_p [0][index] = tau_nux [0][index];
          tau_nux_p [1][index] = tau_nux [1][index];
        }
      }
    }

    l2norm = sqrt(l2norm);
    printf("it. %02d: l2norm = %.15e\n", n, l2norm);
    if(n>0 && l2norm<1e-16)
      break;

    // Step 4.b: Update optical depth using path of least resistance algorithm
    const double dxx[3] = {dx, dy, dz};
    const double stencil_gxx[3] = {1, 1, 1}; // Flat spatial metric
    const double stencil_gyy[3] = {1, 1, 1}; // Flat spatial metric
    const double stencil_gzz[3] = {1, 1, 1}; // Flat spatial metric
#pragma omp parallel for
    for(int i2=Ng2;i2<Nt2-Ng2;i2++) {
      for(int i1=Ng1;i1<Nt1-Ng1;i1++) {
        for(int i0=Ng0;i0<Nt0-Ng0;i0++) {

          ghl_neutrino_optical_depths tau_i_j_k;
          update_optical_depth_cell(
                i0, i1, i2, Nt0, Nt1, dxx, stencil_gxx, stencil_gyy, stencil_gzz,
                kappa_nue, kappa_anue, kappa_nux, tau_nue_p, tau_anue_p, tau_nux_p,
                &tau_i_j_k);

          const int i_j_k = IDX3D(i0, i1, i2);

          tau_nue [0][i_j_k] = tau_i_j_k.nue [0];
          tau_nue [1][i_j_k] = tau_i_j_k.nue [1];
          tau_anue[0][i_j_k] = tau_i_j_k.anue[0];
          tau_anue[1][i_j_k] = tau_i_j_k.anue[1];
          tau_nux [0][i_j_k] = tau_i_j_k.nux [0];
          tau_nux [1][i_j_k] = tau_i_j_k.nux [1];
        }
      }
    }
  }

  // Step 5: Dump the data
  if( test_key != 2 ) {
    FILE *fp;
    if( test_key )
      fp = fopen_with_check("nrpyleakage_constant_density_sphere_perturbed.bin", "wb");
    else
      fp = fopen_with_check("nrpyleakage_constant_density_sphere_unperturbed.bin", "wb");

    fwrite(kappa_nue [0], sizeof(double), Ntotal, fp);
    fwrite(kappa_nue [1], sizeof(double), Ntotal, fp);
    fwrite(kappa_anue[0], sizeof(double), Ntotal, fp);
    fwrite(kappa_anue[1], sizeof(double), Ntotal, fp);
    fwrite(kappa_nux [0], sizeof(double), Ntotal, fp);
    fwrite(kappa_nux [1], sizeof(double), Ntotal, fp);
    fwrite(tau_nue   [0], sizeof(double), Ntotal, fp);
    fwrite(tau_nue   [1], sizeof(double), Ntotal, fp);
    fwrite(tau_anue  [0], sizeof(double), Ntotal, fp);
    fwrite(tau_anue  [1], sizeof(double), Ntotal, fp);
    fwrite(tau_nux   [0], sizeof(double), Ntotal, fp);
    fwrite(tau_nux   [1], sizeof(double), Ntotal, fp);
    fclose(fp);
  }
  else {
    double **kappa_nue_unpert  = (double **)malloc(sizeof(double *)*2);
    double **kappa_anue_unpert = (double **)malloc(sizeof(double *)*2);
    double **kappa_nux_unpert  = (double **)malloc(sizeof(double *)*2);
    double **tau_nue_unpert    = (double **)malloc(sizeof(double *)*2);
    double **tau_anue_unpert   = (double **)malloc(sizeof(double *)*2);
    double **tau_nux_unpert    = (double **)malloc(sizeof(double *)*2);
    double **kappa_nue_pert    = (double **)malloc(sizeof(double *)*2);
    double **kappa_anue_pert   = (double **)malloc(sizeof(double *)*2);
    double **kappa_nux_pert    = (double **)malloc(sizeof(double *)*2);
    double **tau_nue_pert      = (double **)malloc(sizeof(double *)*2);
    double **tau_anue_pert     = (double **)malloc(sizeof(double *)*2);
    double **tau_nux_pert      = (double **)malloc(sizeof(double *)*2);
    for(int i=0;i<2;i++) {
      kappa_nue_unpert [i] = (double *)malloc(sizeof(double)*Ntotal);
      kappa_anue_unpert[i] = (double *)malloc(sizeof(double)*Ntotal);
      kappa_nux_unpert [i] = (double *)malloc(sizeof(double)*Ntotal);
      tau_nue_unpert   [i] = (double *)malloc(sizeof(double)*Ntotal);
      tau_anue_unpert  [i] = (double *)malloc(sizeof(double)*Ntotal);
      tau_nux_unpert   [i] = (double *)malloc(sizeof(double)*Ntotal);
      kappa_nue_pert   [i] = (double *)malloc(sizeof(double)*Ntotal);
      kappa_anue_pert  [i] = (double *)malloc(sizeof(double)*Ntotal);
      kappa_nux_pert   [i] = (double *)malloc(sizeof(double)*Ntotal);
      tau_nue_pert     [i] = (double *)malloc(sizeof(double)*Ntotal);
      tau_anue_pert    [i] = (double *)malloc(sizeof(double)*Ntotal);
      tau_nux_pert     [i] = (double *)malloc(sizeof(double)*Ntotal);
    }
    // Read in the unperturbed data
    FILE *fp_unpert = fopen_with_check("nrpyleakage_constant_density_sphere_unperturbed.bin", "rb");
    int err = 0;
    err += fread(kappa_nue_unpert [0], sizeof(double), Ntotal, fp_unpert);
    err += fread(kappa_nue_unpert [1], sizeof(double), Ntotal, fp_unpert);
    err += fread(kappa_anue_unpert[0], sizeof(double), Ntotal, fp_unpert);
    err += fread(kappa_anue_unpert[1], sizeof(double), Ntotal, fp_unpert);
    err += fread(kappa_nux_unpert [0], sizeof(double), Ntotal, fp_unpert);
    err += fread(kappa_nux_unpert [1], sizeof(double), Ntotal, fp_unpert);
    err += fread(tau_nue_unpert   [0], sizeof(double), Ntotal, fp_unpert);
    err += fread(tau_nue_unpert   [1], sizeof(double), Ntotal, fp_unpert);
    err += fread(tau_anue_unpert  [0], sizeof(double), Ntotal, fp_unpert);
    err += fread(tau_anue_unpert  [1], sizeof(double), Ntotal, fp_unpert);
    err += fread(tau_nux_unpert   [0], sizeof(double), Ntotal, fp_unpert);
    err += fread(tau_nux_unpert   [1], sizeof(double), Ntotal, fp_unpert);
    fclose(fp_unpert);
    if( err != 12*Ntotal )
      ghl_error("Error reading file nrpyleakage_constant_density_sphere_unperturbed.bin\n");

    // Read in the perturbed data
    FILE *fp_pert   = fopen_with_check("nrpyleakage_constant_density_sphere_perturbed.bin", "rb");
    err = 0;
    err += fread(kappa_nue_pert [0], sizeof(double), Ntotal, fp_pert);
    err += fread(kappa_nue_pert [1], sizeof(double), Ntotal, fp_pert);
    err += fread(kappa_anue_pert[0], sizeof(double), Ntotal, fp_pert);
    err += fread(kappa_anue_pert[1], sizeof(double), Ntotal, fp_pert);
    err += fread(kappa_nux_pert [0], sizeof(double), Ntotal, fp_pert);
    err += fread(kappa_nux_pert [1], sizeof(double), Ntotal, fp_pert);
    err += fread(tau_nue_pert   [0], sizeof(double), Ntotal, fp_pert);
    err += fread(tau_nue_pert   [1], sizeof(double), Ntotal, fp_pert);
    err += fread(tau_anue_pert  [0], sizeof(double), Ntotal, fp_pert);
    err += fread(tau_anue_pert  [1], sizeof(double), Ntotal, fp_pert);
    err += fread(tau_nux_pert   [0], sizeof(double), Ntotal, fp_pert);
    err += fread(tau_nux_pert   [1], sizeof(double), Ntotal, fp_pert);
    if( err != 12*Ntotal )
      ghl_error("Error reading file nrpyleakage_constant_density_sphere_perturbed.bin\n");
    fclose(fp_pert);

    // Perform the validation
    for(int i2=0;i2<Nt2;i2++) {
      for(int i1=0;i1<Nt1;i1++) {
        for(int i0=0;i0<Nt0;i0++) {
          const int index = IDX3D(i0,i1,i2);
          if(sphere_pert_test_fail(
                   kappa_nue_unpert[0][index], kappa_nue[0][index],
                   kappa_nue_pert[0][index])) {
            ghl_error(
                  "Validation failed for kappa_nue[0] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, kappa_nue_unpert[0][index], kappa_nue[0][index],
                  kappa_nue_pert[0][index]);
          }
          if(sphere_pert_test_fail(
                   kappa_nue_unpert[1][index], kappa_nue[1][index],
                   kappa_nue_pert[1][index])) {
            ghl_error(
                  "Validation failed for kappa_nue[1] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, kappa_nue_unpert[1][index], kappa_nue[1][index],
                  kappa_nue_pert[1][index]);
          }
          if(sphere_pert_test_fail(
                   kappa_anue_unpert[0][index], kappa_anue[0][index],
                   kappa_anue_pert[0][index])) {
            ghl_error(
                  "Validation failed for kappa_anue[0] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, kappa_anue_unpert[0][index], kappa_anue[0][index],
                  kappa_anue_pert[0][index]);
          }
          if(sphere_pert_test_fail(
                   kappa_anue_unpert[1][index], kappa_anue[1][index],
                   kappa_anue_pert[1][index])) {
            ghl_error(
                  "Validation failed for kappa_anue[1] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, kappa_anue_unpert[1][index], kappa_anue[1][index],
                  kappa_anue_pert[1][index]);
          }
          if(sphere_pert_test_fail(
                   kappa_nux_unpert[0][index], kappa_nux[0][index],
                   kappa_nux_pert[0][index])) {
            ghl_error(
                  "Validation failed for kappa_nux[0] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, kappa_nux_unpert[0][index], kappa_nux[0][index],
                  kappa_nux_pert[0][index]);
          }
          if(sphere_pert_test_fail(
                   kappa_nux_unpert[1][index], kappa_nux[1][index],
                   kappa_nux_pert[1][index])) {
            ghl_error(
                  "Validation failed for kappa_nux[1] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, kappa_nux_unpert[1][index], kappa_nux[1][index],
                  kappa_nux_pert[1][index]);
          }
          if(sphere_pert_test_fail(
                   tau_nue_unpert[0][index], tau_nue[0][index],
                   tau_nue_pert[0][index])) {
            ghl_error(
                  "Validation failed for tau_nue[0] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, tau_nue_unpert[0][index], tau_nue[0][index],
                  tau_nue_pert[0][index]);
          }
          if(sphere_pert_test_fail(
                   tau_nue_unpert[1][index], tau_nue[1][index],
                   tau_nue_pert[1][index])) {
            ghl_error(
                  "Validation failed for tau_nue[1] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, tau_nue_unpert[1][index], tau_nue[1][index],
                  tau_nue_pert[1][index]);
          }
          if(sphere_pert_test_fail(
                   tau_anue_unpert[0][index], tau_anue[0][index],
                   tau_anue_pert[0][index])) {
            ghl_error(
                  "Validation failed for tau_anue[0] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, tau_anue_unpert[0][index], tau_anue[0][index],
                  tau_anue_pert[0][index]);
          }
          if(sphere_pert_test_fail(
                   tau_anue_unpert[1][index], tau_anue[1][index],
                   tau_anue_pert[1][index])) {
            ghl_error(
                  "Validation failed for tau_anue[1] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, tau_anue_unpert[1][index], tau_anue[1][index],
                  tau_anue_pert[1][index]);
          }
          if(sphere_pert_test_fail(
                   tau_nux_unpert[0][index], tau_nux[0][index],
                   tau_nux_pert[0][index])) {
            ghl_error(
                  "Validation failed for tau_nux[0] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, tau_nux_unpert[0][index], tau_nux[0][index],
                  tau_nux_pert[0][index]);
          }
          if(sphere_pert_test_fail(
                   tau_nux_unpert[1][index], tau_nux[1][index],
                   tau_nux_pert[1][index])) {
            ghl_error(
                  "Validation failed for tau_nux[1] at (i0, i1, i2) = (%d, %d, %d), "
                  "index %d: "
                  "trusted %.17e, computed %.17e, perturbed %.17e\n",
                  i0, i1, i2, index, tau_nux_unpert[1][index], tau_nux[1][index],
                  tau_nux_pert[1][index]);
          }
        }
      }
    }

    for(int i=0;i<2;i++) {
      free(kappa_nue_unpert [i]);
      free(kappa_anue_unpert[i]);
      free(kappa_nux_unpert [i]);
      free(tau_nue_unpert   [i]);
      free(tau_anue_unpert  [i]);
      free(tau_nux_unpert   [i]);
      free(kappa_nue_pert   [i]);
      free(kappa_anue_pert  [i]);
      free(kappa_nux_pert   [i]);
      free(tau_nue_pert     [i]);
      free(tau_anue_pert    [i]);
      free(tau_nux_pert     [i]);
    }
    free(kappa_nue_unpert );
    free(kappa_anue_unpert);
    free(kappa_nux_unpert );
    free(tau_nue_unpert   );
    free(tau_anue_unpert  );
    free(tau_nux_unpert   );
    free(kappa_nue_pert   );
    free(kappa_anue_pert  );
    free(kappa_nux_pert   );
    free(tau_nue_pert     );
    free(tau_anue_pert    );
    free(tau_nux_pert     );
  }

  // Step 6: Free memory
  for(int i=0;i<2;i++) {
    free(kappa_nue [i]);
    free(kappa_anue[i]);
    free(kappa_nux [i]);
    free(tau_nue   [i]);
    free(tau_anue  [i]);
    free(tau_nux   [i]);
    free(tau_nue_p [i]);
    free(tau_anue_p[i]);
    free(tau_nux_p [i]);
  }
  free(kappa_nue);
  free(kappa_anue);
  free(kappa_nux);
  free(tau_nue);
  free(tau_anue);
  free(tau_nux);
  free(tau_nue_p);
  free(tau_anue_p);
  free(tau_nux_p);
}

void
generate_test_data(const ghl_eos_parameters *restrict eos) {
  for(int perturb=0;perturb<=1;perturb++)
    constantdensitysphere_test(eos, perturb);
}

void
run_unit_test(const ghl_eos_parameters *restrict eos) {
  check_sphere_caller_neighbor_order();
  if(!sphere_pert_test_fail(0.0, DBL_MIN, 0.0) || sphere_pert_test_fail(0.0, 0.0, 0.0)) {
    ghl_error("Sphere comparison did not enforce exact-zero references\n");
  }
  if(!sphere_pert_test_fail(1.0, NAN, NAN) || !sphere_pert_test_fail(NAN, NAN, 1.0)) {
    ghl_error("Sphere comparison accepted nonfinite values\n");
  }
  constantdensitysphere_test(eos, 2);
}

#include "nrpyleakage_main.h"
