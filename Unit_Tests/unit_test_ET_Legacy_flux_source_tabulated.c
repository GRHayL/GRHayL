#include "unit_tests.h"

#define AM2 -0.0625
#define AM1 0.5625
#define A0 0.5625
#define A1 -0.0625
#define COMPUTE_FCVAL(METRICm2, METRICm1, METRIC, METRICp1)                                       \
  (AM2 * (METRICm2) + AM1 * (METRICm1) + A0 * (METRIC) + A1 * (METRICp1))

static inline void calculate_face_value(
    const int flux_dirn,
    const int dirlength,
    const int ghostzone,
    const double *restrict cell_var,
    double *restrict face_var) {

  const int xdir = (flux_dirn == 0);
  const int ydir = (flux_dirn == 1);
  const int zdir = (flux_dirn == 2);

  for(int k = ghostzone - 1; k < dirlength - (ghostzone - 2); k++) {
    for(int j = ghostzone - 1; j < dirlength - (ghostzone - 2); j++) {
      for(int i = ghostzone - 1; i < dirlength - (ghostzone - 2); i++) {
        const int indm2 = indexf(dirlength, i - 2 * xdir, j - 2 * ydir, k - 2 * zdir);
        const int indm1 = indexf(dirlength, i - xdir, j - ydir, k - zdir);
        const int index = indexf(dirlength, i, j, k);
        const int indp1 = indexf(dirlength, i + xdir, j + ydir, k + zdir);

        face_var[index]
            = COMPUTE_FCVAL(cell_var[indm2], cell_var[indm1], cell_var[index], cell_var[indp1]);
      }
    }
  }
}

int main(int argc, char **argv) {

  if(argc != 2) {
    ghl_info("Usage: %s <EOS table path>\n", argv[0]);
    exit(1);
  }

  const char *tablepath = argv[1];
  const double rho_b_atm = 1e-12;
  const double rho_b_min = -1;
  const double rho_b_max = -1;
  const double Y_e_atm = 0.5;
  const double Y_e_min = -1;
  const double Y_e_max = -1;
  const double T_atm = 1e-2;
  const double T_min = -1;
  const double T_max = -1;

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
      tablepath, rho_b_atm, rho_b_min, rho_b_max, Y_e_atm, Y_e_min, Y_e_max, T_atm, T_min, T_max,
      &eos);

  FILE *infile = fopen_with_check("ET_Legacy_flux_source_input.bin", "rb");

  int dirlength;
  int key = fread(&dirlength, sizeof(int), 1, infile);
  const int ghostzone = 3;
  const int arraylength = dirlength * dirlength * dirlength;
  const double invdx = 1.0 / 0.1;
  const double invdy = 1.0 / 0.1;
  const double invdz = 1.0 / 0.1;
  const double poison = 0.0 / 0.0;

  // Allocate memory for metric
  double *lapse = (double *)malloc(sizeof(double) * arraylength);
  double *betax = (double *)malloc(sizeof(double) * arraylength);
  double *betay = (double *)malloc(sizeof(double) * arraylength);
  double *betaz = (double *)malloc(sizeof(double) * arraylength);
  double *gxx = (double *)malloc(sizeof(double) * arraylength);
  double *gxy = (double *)malloc(sizeof(double) * arraylength);
  double *gxz = (double *)malloc(sizeof(double) * arraylength);
  double *gyy = (double *)malloc(sizeof(double) * arraylength);
  double *gyz = (double *)malloc(sizeof(double) * arraylength);
  double *gzz = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for face-interpolated metric
  double *face_lapse = (double *)malloc(sizeof(double) * arraylength);
  double *face_betax = (double *)malloc(sizeof(double) * arraylength);
  double *face_betay = (double *)malloc(sizeof(double) * arraylength);
  double *face_betaz = (double *)malloc(sizeof(double) * arraylength);
  double *face_gxx = (double *)malloc(sizeof(double) * arraylength);
  double *face_gxy = (double *)malloc(sizeof(double) * arraylength);
  double *face_gxz = (double *)malloc(sizeof(double) * arraylength);
  double *face_gyy = (double *)malloc(sizeof(double) * arraylength);
  double *face_gyz = (double *)malloc(sizeof(double) * arraylength);
  double *face_gzz = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for extrinsic curvature
  double *kxx = (double *)malloc(sizeof(double) * arraylength);
  double *kxy = (double *)malloc(sizeof(double) * arraylength);
  double *kxz = (double *)malloc(sizeof(double) * arraylength);
  double *kyy = (double *)malloc(sizeof(double) * arraylength);
  double *kyz = (double *)malloc(sizeof(double) * arraylength);
  double *kzz = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for primitives
  double *rho = (double *)malloc(sizeof(double) * arraylength);
  double *press = (double *)malloc(sizeof(double) * arraylength);
  double *vx = (double *)malloc(sizeof(double) * arraylength);
  double *vy = (double *)malloc(sizeof(double) * arraylength);
  double *vz = (double *)malloc(sizeof(double) * arraylength);
  double *Bx = (double *)malloc(sizeof(double) * arraylength);
  double *By = (double *)malloc(sizeof(double) * arraylength);
  double *Bz = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e = (double *)malloc(sizeof(double) * arraylength);
  double *entropy = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for right face
  double *rho_r = (double *)malloc(sizeof(double) * arraylength);
  double *press_r = (double *)malloc(sizeof(double) * arraylength);
  double *vx_r = (double *)malloc(sizeof(double) * arraylength);
  double *vy_r = (double *)malloc(sizeof(double) * arraylength);
  double *vz_r = (double *)malloc(sizeof(double) * arraylength);
  double *Bx_r = (double *)malloc(sizeof(double) * arraylength);
  double *By_r = (double *)malloc(sizeof(double) * arraylength);
  double *Bz_r = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e_r = (double *)malloc(sizeof(double) * arraylength);
  double *entropy_r = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for left face
  double *rho_l = (double *)malloc(sizeof(double) * arraylength);
  double *press_l = (double *)malloc(sizeof(double) * arraylength);
  double *vx_l = (double *)malloc(sizeof(double) * arraylength);
  double *vy_l = (double *)malloc(sizeof(double) * arraylength);
  double *vz_l = (double *)malloc(sizeof(double) * arraylength);
  double *Bx_l = (double *)malloc(sizeof(double) * arraylength);
  double *By_l = (double *)malloc(sizeof(double) * arraylength);
  double *Bz_l = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e_l = (double *)malloc(sizeof(double) * arraylength);
  double *entropy_l = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for fluxes
  double *tau_flux = (double *)malloc(sizeof(double) * arraylength);
  double *rho_star_flux = (double *)malloc(sizeof(double) * arraylength);
  double *S_x_flux = (double *)malloc(sizeof(double) * arraylength);
  double *S_y_flux = (double *)malloc(sizeof(double) * arraylength);
  double *S_z_flux = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e_star_flux = (double *)malloc(sizeof(double) * arraylength);
  double *S_star_flux = (double *)malloc(sizeof(double) * arraylength);

  // Allocate memory for right-hand sides
  double *rho_star_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *tau_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *S_x_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *S_y_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *S_z_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *Y_e_star_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *S_star_rhs = (double *)malloc(sizeof(double) * arraylength);

  // clang-format off
  int count=0;
  key  = fread(gxx,     sizeof(double), arraylength, infile); count++;
  key += fread(gxy,     sizeof(double), arraylength, infile); count++;
  key += fread(gxz,     sizeof(double), arraylength, infile); count++;
  key += fread(gyy,     sizeof(double), arraylength, infile); count++;
  key += fread(gyz,     sizeof(double), arraylength, infile); count++;
  key += fread(gzz,     sizeof(double), arraylength, infile); count++;
  key += fread(lapse,   sizeof(double), arraylength, infile); count++;
  key += fread(betax,   sizeof(double), arraylength, infile); count++;
  key += fread(betay,   sizeof(double), arraylength, infile); count++;
  key += fread(betaz,   sizeof(double), arraylength, infile); count++;
  key += fread(kxx,     sizeof(double), arraylength, infile); count++;
  key += fread(kxy,     sizeof(double), arraylength, infile); count++;
  key += fread(kxz,     sizeof(double), arraylength, infile); count++;
  key += fread(kyy,     sizeof(double), arraylength, infile); count++;
  key += fread(kyz,     sizeof(double), arraylength, infile); count++;
  key += fread(kzz,     sizeof(double), arraylength, infile); count++;
  key += fread(rho,     sizeof(double), arraylength, infile); count++;
  key += fread(press,   sizeof(double), arraylength, infile); count++;
  key += fread(vx,      sizeof(double), arraylength, infile); count++;
  key += fread(vy,      sizeof(double), arraylength, infile); count++;
  key += fread(vz,      sizeof(double), arraylength, infile); count++;
  key += fread(Bx,      sizeof(double), arraylength, infile); count++;
  key += fread(By,      sizeof(double), arraylength, infile); count++;
  key += fread(Bz,      sizeof(double), arraylength, infile); count++;
	key += fread(Y_e,     sizeof(double), arraylength, infile); count++;
	key += fread(entropy, sizeof(double), arraylength, infile); count++;
  // clang-format on

  if(key != arraylength * count) {
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
              "is up-to-date with current test version.\n");
  }

  for(int k = ghostzone; k < dirlength - ghostzone; k++) {
    for(int j = ghostzone; j < dirlength - ghostzone; j++) {
      for(int i = ghostzone; i < dirlength - ghostzone; i++) {
        const int index = indexf(dirlength, i, j, k);
        rho_star_rhs[index] = 0.0;
        tau_rhs[index] = 0.0;
        S_x_rhs[index] = 0.0;
        S_y_rhs[index] = 0.0;
        S_z_rhs[index] = 0.0;
        Y_e_star_rhs[index] = 0.0;
        S_star_rhs[index] = 0.0;
      }
    }
  }

  // Function pointer to allow for loop over fluxes
  void (*calculate_HLLE_fluxes)(
      ghl_primitive_quantities *restrict, ghl_primitive_quantities *restrict,
      const ghl_eos_parameters *restrict, const ghl_metric_quantities *restrict, const double,
      const double, ghl_conservative_quantities *restrict);

  void (*calculate_characteristic_speed)(
      ghl_primitive_quantities *restrict, ghl_primitive_quantities *restrict,
      const ghl_eos_parameters *restrict, const ghl_metric_quantities *restrict, double *restrict,
      double *restrict);

  // Loop over flux directions (x,y,z)
  for(int flux_dirn = 0; flux_dirn < 3; flux_dirn++) {
    const int xdir = (flux_dirn == 0);
    const int ydir = (flux_dirn == 1);
    const int zdir = (flux_dirn == 2);

    // Set function pointer to specific function for a given direction
    switch(flux_dirn) {
      case 0:
        calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy;
        calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
        break;
      case 1:
        calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy;
        calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
        break;
      case 2:
        calculate_HLLE_fluxes = &ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy;
        calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
        break;
    }

    // clang-format off
    count = 0;
    key  = fread(rho_r,    sizeof(double), arraylength, infile); count++;
    key += fread(press_r,  sizeof(double), arraylength, infile); count++;
    key += fread(vx_r,     sizeof(double), arraylength, infile); count++;
    key += fread(vy_r,     sizeof(double), arraylength, infile); count++;
    key += fread(vz_r,     sizeof(double), arraylength, infile); count++;
    key += fread(Bx_r,     sizeof(double), arraylength, infile); count++;
    key += fread(By_r,     sizeof(double), arraylength, infile); count++;
    key += fread(Bz_r,     sizeof(double), arraylength, infile); count++;
    key += fread(Y_e_r,    sizeof(double), arraylength, infile); count++;
    key += fread(entropy_r,sizeof(double), arraylength, infile); count++;

    key += fread(rho_l,    sizeof(double), arraylength, infile); count++;
    key += fread(press_l,  sizeof(double), arraylength, infile); count++;
    key += fread(vx_l,     sizeof(double), arraylength, infile); count++;
    key += fread(vy_l,     sizeof(double), arraylength, infile); count++;
    key += fread(vz_l,     sizeof(double), arraylength, infile); count++;
    key += fread(Bx_l,     sizeof(double), arraylength, infile); count++;
    key += fread(By_l,     sizeof(double), arraylength, infile); count++;
    key += fread(Bz_l,     sizeof(double), arraylength, infile); count++;
		key += fread(Y_e_l,    sizeof(double), arraylength, infile); count++;
    key += fread(entropy_l,sizeof(double), arraylength, infile); count++;
    // clang-format on

    if(key != arraylength * count) {
      ghl_error("An error has occured with reading in initial data. Please check that data\n"
                "is up-to-date with current test version.\n");
    }

    // Calculate the values of the metric quantities at the faces by interpolating
    // from the cell-centered quantities
    calculate_face_value(flux_dirn, dirlength, ghostzone, lapse, face_lapse);
    calculate_face_value(flux_dirn, dirlength, ghostzone, betax, face_betax);
    calculate_face_value(flux_dirn, dirlength, ghostzone, betay, face_betay);
    calculate_face_value(flux_dirn, dirlength, ghostzone, betaz, face_betaz);
    calculate_face_value(flux_dirn, dirlength, ghostzone, gxx, face_gxx);
    calculate_face_value(flux_dirn, dirlength, ghostzone, gxy, face_gxy);
    calculate_face_value(flux_dirn, dirlength, ghostzone, gxz, face_gxz);
    calculate_face_value(flux_dirn, dirlength, ghostzone, gyy, face_gyy);
    calculate_face_value(flux_dirn, dirlength, ghostzone, gyz, face_gyz);
    calculate_face_value(flux_dirn, dirlength, ghostzone, gzz, face_gzz);

    // Upper bound includes 1 ghostzone for RHS calculation in following
    // loop.
    for(int k = ghostzone; k < dirlength - (ghostzone - 1); k++) {
      for(int j = ghostzone; j < dirlength - (ghostzone - 1); j++) {
        for(int i = ghostzone; i < dirlength - (ghostzone - 1); i++) {
          const int index = indexf(dirlength, i, j, k);

          ghl_metric_quantities metric_face;
          ghl_initialize_metric(
              face_lapse[index], face_betax[index], face_betay[index], face_betaz[index],
              face_gxx[index], face_gxy[index], face_gxz[index], face_gyy[index], face_gyz[index],
              face_gzz[index], &metric_face);

          // B's get rescaled to match IGM's definition of B
          ghl_primitive_quantities prims_r, prims_l;
          ghl_initialize_primitives(
              rho_r[index], press_r[index], poison, vx_r[index], vy_r[index], vz_r[index],
              ONE_OVER_SQRT_4PI * Bx_r[index], ONE_OVER_SQRT_4PI * By_r[index],
              ONE_OVER_SQRT_4PI * Bz_r[index], entropy_r[index], Y_e_r[index], eos.table_T_min,
              &prims_r);

          ghl_initialize_primitives(
              rho_l[index], press_l[index], poison, vx_l[index], vy_l[index], vz_l[index],
              ONE_OVER_SQRT_4PI * Bx_l[index], ONE_OVER_SQRT_4PI * By_l[index],
              ONE_OVER_SQRT_4PI * Bz_l[index], entropy_l[index], Y_e_l[index], eos.table_T_min,
              &prims_l);

          // Generate randomized u^0
          prims_r.u0 = rho_r[index] * Bx_r[index] / vy_r[index];
          prims_l.u0 = rho_l[index] * Bx_l[index] / vy_l[index];

          double cmin, cmax;
          calculate_characteristic_speed(&prims_r, &prims_l, &eos, &metric_face, &cmin, &cmax);

          ghl_conservative_quantities cons_fluxes;
          calculate_HLLE_fluxes(&prims_r, &prims_l, &eos, &metric_face, cmin, cmax, &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux[index] = cons_fluxes.tau;
          S_x_flux[index] = cons_fluxes.SD[0];
          S_y_flux[index] = cons_fluxes.SD[1];
          S_z_flux[index] = cons_fluxes.SD[2];
          Y_e_star_flux[index] = cons_fluxes.Y_e;
          S_star_flux[index] = cons_fluxes.entropy;
        }
      }
    }

    for(int k = ghostzone; k < dirlength - ghostzone; k++) {
      for(int j = ghostzone; j < dirlength - ghostzone; j++) {
        for(int i = ghostzone; i < dirlength - ghostzone; i++) {
          const int index = indexf(dirlength, i, j, k);
          const int indp1 = indexf(dirlength, i + xdir, j + ydir, k + zdir);

          rho_star_rhs[index] += invdx * (rho_star_flux[index] - rho_star_flux[indp1]);
          tau_rhs[index] += invdx * (tau_flux[index] - tau_flux[indp1]);
          S_x_rhs[index] += invdx * (S_x_flux[index] - S_x_flux[indp1]);
          S_y_rhs[index] += invdx * (S_y_flux[index] - S_y_flux[indp1]);
          S_z_rhs[index] += invdx * (S_z_flux[index] - S_z_flux[indp1]);
          Y_e_star_rhs[index] += invdx * (Y_e_star_flux[index] - Y_e_star_flux[indp1]);
          S_star_rhs[index] += invdx * (S_star_flux[index] - S_star_flux[indp1]);
        }
      }
    }
  }

  for(int k = ghostzone; k < dirlength - ghostzone; k++) {
    for(int j = ghostzone; j < dirlength - ghostzone; j++) {
      for(int i = ghostzone; i < dirlength - ghostzone; i++) {
        const int index = indexf(dirlength, i, j, k);
        const int indip1 = indexf(dirlength, i + 1, j, k);
        const int indjp1 = indexf(dirlength, i, j + 1, k);
        const int indkp1 = indexf(dirlength, i, j, k + 1);

        ghl_metric_quantities metric;
        ghl_initialize_metric(
            lapse[index], betax[index], betay[index], betaz[index], gxx[index], gxy[index],
            gxz[index], gyy[index], gyz[index], gzz[index], &metric);

        ghl_metric_quantities metric_derivs_x;
        metric_derivs_x.lapse = invdx * (face_lapse[indip1] - face_lapse[index]);
        metric_derivs_x.betaU[0] = invdx * (face_betax[indip1] - face_betax[index]);
        metric_derivs_x.betaU[1] = invdx * (face_betay[indip1] - face_betay[index]);
        metric_derivs_x.betaU[2] = invdx * (face_betaz[indip1] - face_betaz[index]);
        metric_derivs_x.gammaDD[0][0] = invdx * (face_gxx[indip1] - face_gxx[index]);
        metric_derivs_x.gammaDD[0][1] = invdx * (face_gxy[indip1] - face_gxy[index]);
        metric_derivs_x.gammaDD[0][2] = invdx * (face_gxz[indip1] - face_gxz[index]);
        metric_derivs_x.gammaDD[1][1] = invdx * (face_gyy[indip1] - face_gyy[index]);
        metric_derivs_x.gammaDD[1][2] = invdx * (face_gyz[indip1] - face_gyz[index]);
        metric_derivs_x.gammaDD[2][2] = invdx * (face_gzz[indip1] - face_gzz[index]);

        ghl_metric_quantities metric_derivs_y;
        metric_derivs_y.lapse = invdy * (face_lapse[indjp1] - face_lapse[index]);
        metric_derivs_y.betaU[0] = invdy * (face_betax[indjp1] - face_betax[index]);
        metric_derivs_y.betaU[1] = invdy * (face_betay[indjp1] - face_betay[index]);
        metric_derivs_y.betaU[2] = invdy * (face_betaz[indjp1] - face_betaz[index]);
        metric_derivs_y.gammaDD[0][0] = invdy * (face_gxx[indjp1] - face_gxx[index]);
        metric_derivs_y.gammaDD[0][1] = invdy * (face_gxy[indjp1] - face_gxy[index]);
        metric_derivs_y.gammaDD[0][2] = invdy * (face_gxz[indjp1] - face_gxz[index]);
        metric_derivs_y.gammaDD[1][1] = invdy * (face_gyy[indjp1] - face_gyy[index]);
        metric_derivs_y.gammaDD[1][2] = invdy * (face_gyz[indjp1] - face_gyz[index]);
        metric_derivs_y.gammaDD[2][2] = invdy * (face_gzz[indjp1] - face_gzz[index]);

        ghl_metric_quantities metric_derivs_z;
        metric_derivs_z.lapse = invdz * (face_lapse[indkp1] - face_lapse[index]);
        metric_derivs_z.betaU[0] = invdz * (face_betax[indkp1] - face_betax[index]);
        metric_derivs_z.betaU[1] = invdz * (face_betay[indkp1] - face_betay[index]);
        metric_derivs_z.betaU[2] = invdz * (face_betaz[indkp1] - face_betaz[index]);
        metric_derivs_z.gammaDD[0][0] = invdz * (face_gxx[indkp1] - face_gxx[index]);
        metric_derivs_z.gammaDD[0][1] = invdz * (face_gxy[indkp1] - face_gxy[index]);
        metric_derivs_z.gammaDD[0][2] = invdz * (face_gxz[indkp1] - face_gxz[index]);
        metric_derivs_z.gammaDD[1][1] = invdz * (face_gyy[indkp1] - face_gyy[index]);
        metric_derivs_z.gammaDD[1][2] = invdz * (face_gyz[indkp1] - face_gyz[index]);
        metric_derivs_z.gammaDD[2][2] = invdz * (face_gzz[indkp1] - face_gzz[index]);

        ghl_extrinsic_curvature curv;
        ghl_initialize_extrinsic_curvature(
            kxx[index], kxy[index], kxz[index], kyy[index], kyz[index], kzz[index], &curv);

        // B's get rescaled to match IGM's definition of B
        ghl_primitive_quantities prims;
        ghl_initialize_primitives(
            rho[index], press[index], poison, vx[index], vy[index], vz[index],
            ONE_OVER_SQRT_4PI * Bx[index], ONE_OVER_SQRT_4PI * By[index],
            ONE_OVER_SQRT_4PI * Bz[index], entropy[index], Y_e[index], eos.table_T_min, &prims);
        prims.u0 = rho[index] * Bx[index] / vy[index];

        ghl_conservative_quantities cons_sources;
        ghl_calculate_source_terms(
            &eos, &prims, &metric, &metric_derivs_x, &metric_derivs_y, &metric_derivs_z, &curv,
            &cons_sources);

        tau_rhs[index] += cons_sources.tau;
        S_x_rhs[index] += cons_sources.SD[0];
        S_y_rhs[index] += cons_sources.SD[1];
        S_z_rhs[index] += cons_sources.SD[2];
      }
    }
  }

  fclose(infile);

  // Allocate memory for comparison data
  double *trusted_rho_star_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *trusted_tau_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *trusted_S_x_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *trusted_S_y_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *trusted_S_z_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *trusted_Y_e_star_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *trusted_S_star_rhs = (double *)malloc(sizeof(double) * arraylength);

  double *pert_rho_star_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *pert_tau_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *pert_S_x_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *pert_S_y_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *pert_S_z_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *pert_Y_e_star_rhs = (double *)malloc(sizeof(double) * arraylength);
  double *pert_S_star_rhs = (double *)malloc(sizeof(double) * arraylength);

  FILE *outfile = fopen_with_check("ET_Legacy_flux_source_output.bin", "rb");

  // clang-format off
	count = 0;
  key  = fread(trusted_rho_star_rhs, sizeof(double), arraylength, outfile); count++;
  key += fread(trusted_tau_rhs,      sizeof(double), arraylength, outfile); count++;
  key += fread(trusted_S_x_rhs,      sizeof(double), arraylength, outfile); count++;
  key += fread(trusted_S_y_rhs,      sizeof(double), arraylength, outfile); count++;
  key += fread(trusted_S_z_rhs,      sizeof(double), arraylength, outfile); count++;
	key += fread(trusted_Y_e_star_rhs, sizeof(double), arraylength, outfile); count++;
	key += fread(trusted_S_star_rhs,   sizeof(double), arraylength, outfile); count++;
  // clang-format on

  if(key != arraylength * count) {
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
              "is up-to-date with current test version.\n");
  }

  fclose(outfile);

  outfile = fopen_with_check("ET_Legacy_flux_source_output_pert.bin", "rb");

  // clang-format off
	count = 0;
  key  = fread(pert_rho_star_rhs, sizeof(double), arraylength, outfile); count++;
  key += fread(pert_tau_rhs,      sizeof(double), arraylength, outfile); count++;
  key += fread(pert_S_x_rhs,      sizeof(double), arraylength, outfile); count++;
  key += fread(pert_S_y_rhs,      sizeof(double), arraylength, outfile); count++;
  key += fread(pert_S_z_rhs,      sizeof(double), arraylength, outfile); count++;
	key += fread(pert_Y_e_star_rhs, sizeof(double), arraylength, outfile); count++;
	key += fread(pert_S_star_rhs,   sizeof(double), arraylength, outfile); count++;
  // clang-format on

  if(key != arraylength * count) {
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
              "is up-to-date with current test version.\n");
  }

  fclose(outfile);

  for(int k = ghostzone; k < dirlength - ghostzone; k++) {
    for(int j = ghostzone; j < dirlength - ghostzone; j++) {
      for(int i = ghostzone; i < dirlength - ghostzone; i++) {
        const int index = indexf(dirlength, i, j, k);

        if(ghl_pert_test_fail(
               trusted_rho_star_rhs[index], rho_star_rhs[index], pert_rho_star_rhs[index])) {
          ghl_error(
              "Test unit_test_ET_Legacy_flux_source has failed for variable rho_star_rhs.\n"
              "  rho_star_rhs trusted %.14e computed %.14e perturbed %.14e\n"
              "  rel.err. %.14e %.14e\n",
              trusted_rho_star_rhs[index], rho_star_rhs[index], pert_rho_star_rhs[index],
              relative_error(trusted_rho_star_rhs[index], rho_star_rhs[index]),
              relative_error(trusted_rho_star_rhs[index], pert_rho_star_rhs[index]));
        }
        if(ghl_pert_test_fail(trusted_tau_rhs[index], tau_rhs[index], pert_tau_rhs[index])) {
          ghl_error(
              "Test unit_test_ET_Legacy_flux_source has failed for variable tau_rhs.\n"
              "  tau_rhs trusted %.14e computed %.14e perturbed %.14e\n"
              "  rel.err. %.14e %.14e\n",
              trusted_tau_rhs[index], tau_rhs[index], pert_tau_rhs[index],
              relative_error(trusted_tau_rhs[index], tau_rhs[index]),
              relative_error(trusted_tau_rhs[index], pert_tau_rhs[index]));
        }
        if(ghl_pert_test_fail(trusted_S_x_rhs[index], S_x_rhs[index], pert_S_x_rhs[index])) {
          ghl_error(
              "Test unit_test_ET_Legacy_flux_source has failed for variable S_x_rhs.\n"
              "  S_x_rhs trusted %.14e computed %.14e perturbed %.14e\n"
              "  rel.err. %.14e %.14e\n",
              trusted_S_x_rhs[index], S_x_rhs[index], pert_S_x_rhs[index],
              relative_error(trusted_S_x_rhs[index], S_x_rhs[index]),
              relative_error(trusted_S_x_rhs[index], pert_S_x_rhs[index]));
        }
        if(ghl_pert_test_fail(trusted_S_y_rhs[index], S_y_rhs[index], pert_S_y_rhs[index])) {
          ghl_error(
              "Test unit_test_ET_Legacy_flux_source has failed for variable S_y_rhs.\n"
              "  S_y_rhs trusted %.14e computed %.14e perturbed %.14e\n"
              "  rel.err. %.14e %.14e\n",
              trusted_S_y_rhs[index], S_y_rhs[index], pert_S_y_rhs[index],
              relative_error(trusted_S_y_rhs[index], S_y_rhs[index]),
              relative_error(trusted_S_y_rhs[index], pert_S_y_rhs[index]));
        }
        if(ghl_pert_test_fail(trusted_S_z_rhs[index], S_z_rhs[index], pert_S_z_rhs[index])) {
          ghl_error(
              "Test unit_test_ET_Legacy_flux_source has failed for variable S_z_rhs.\n"
              "  S_z_rhs trusted %.14e computed %.14e perturbed %.14e\n"
              "  rel.err. %.14e %.14e\n",
              trusted_S_z_rhs[index], S_z_rhs[index], pert_S_z_rhs[index],
              relative_error(trusted_S_z_rhs[index], S_z_rhs[index]),
              relative_error(trusted_S_z_rhs[index], pert_S_z_rhs[index]));
        }
        if(ghl_pert_test_fail(
               trusted_Y_e_star_rhs[index], Y_e_star_rhs[index], pert_Y_e_star_rhs[index])) {
          ghl_error(
              "Test unit_test_ET_Legacy_flux_source has failed for variable Y_e_star_rhs.\n"
              "  Y_e_star_rhs trusted %.14e computed %.14e perturbed %.14e\n"
              "  rel.err. %.14e %.14e\n",
              trusted_Y_e_star_rhs[index], Y_e_star_rhs[index], pert_Y_e_star_rhs[index],
              relative_error(trusted_Y_e_star_rhs[index], Y_e_star_rhs[index]),
              relative_error(trusted_Y_e_star_rhs[index], pert_Y_e_star_rhs[index]));
        }
        if(ghl_pert_test_fail(
               trusted_S_star_rhs[index], S_star_rhs[index], pert_S_star_rhs[index])) {
          ghl_error(
              "Test unit_test_ET_Legacy_flux_source has failed for variable S_star_rhs.\n"
              "  S_star_rhs trusted %.14e computed %.14e perturbed %.14e\n"
              "  rel.err. %.14e %.14e\n",
              trusted_S_star_rhs[index], S_star_rhs[index], pert_S_star_rhs[index],
              relative_error(trusted_S_star_rhs[index], S_star_rhs[index]),
              relative_error(trusted_S_star_rhs[index], pert_S_star_rhs[index]));
        }
      }
    }
    ghl_info("ET_Legacy flux/source test has passed!\n");

    free(gxx);
    free(gxy);
    free(gxz);
    free(gyy);
    free(gyz);
    free(gzz);
    free(lapse);
    free(betax);
    free(betay);
    free(betaz);
    free(face_gxx);
    free(face_gxy);
    free(face_gxz);
    free(face_gyy);
    free(face_gyz);
    free(face_gzz);
    free(face_lapse);
    free(face_betax);
    free(face_betay);
    free(face_betaz);
    free(kxx);
    free(kxy);
    free(kxz);
    free(kyy);
    free(kyz);
    free(kzz);
    free(rho);
    free(press);
    free(vx);
    free(vy);
    free(vz);
    free(Bx);
    free(By);
    free(Bz);
    free(Y_e);
    free(entropy);
    free(rho_r);
    free(press_r);
    free(vx_r);
    free(vy_r);
    free(vz_r);
    free(Bx_r);
    free(By_r);
    free(Bz_r);
    free(Y_e_r);
    free(entropy_r);
    free(rho_l);
    free(press_l);
    free(vx_l);
    free(vy_l);
    free(vz_l);
    free(Bx_l);
    free(By_l);
    free(Bz_l);
    free(Y_e_l);
    free(entropy_l);
    free(rho_star_flux);
    free(tau_flux);
    free(S_x_flux);
    free(S_y_flux);
    free(S_z_flux);
    free(Y_e_star_flux);
    free(S_star_flux);
    free(rho_star_rhs);
    free(tau_rhs);
    free(S_x_rhs);
    free(S_y_rhs);
    free(S_z_rhs);
    free(Y_e_star_rhs);
    free(S_star_rhs);
    free(trusted_rho_star_rhs);
    free(trusted_tau_rhs);
    free(trusted_S_x_rhs);
    free(trusted_S_y_rhs);
    free(trusted_S_z_rhs);
    free(trusted_Y_e_star_rhs);
    free(trusted_S_star_rhs);
    free(pert_rho_star_rhs);
    free(pert_tau_rhs);
    free(pert_S_x_rhs);
    free(pert_S_y_rhs);
    free(pert_S_z_rhs);
    free(pert_Y_e_star_rhs);
    free(pert_S_star_rhs);
  }

  return 0;
}
