#include "unit_tests.h"
#include "flux_source_unit_test.h"

static inline void compute_h_and_cs2(struct eos_parameters const *restrict eos,
                                     primitive_quantities const *restrict prims,
                                     double *restrict h,
                                     double *restrict cs2) {


// CCTK_REAL h_ = U[PRESSURE]*U[VX]/U[VZ];
*h = prims->press*prims->vx / prims->vz;

// CCTK_REAL c_s_squared  = U[PRESSURE]*U[VZ]/(h);
*cs2 = prims->rho*prims->vz*(*h)/1e4;
// printf("works!!\n");
}

#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

static inline void calculate_face_value(
      const int flux_dirn,
      const int dirlength,
      const int ghostzone,
      const double *restrict cell_var,
      double *restrict face_var) {

  const int xdir = (flux_dirn == 0);
  const int ydir = (flux_dirn == 1);
  const int zdir = (flux_dirn == 2);

#pragma omp parallel for
  for(int k=ghostzone-1; k<dirlength-(ghostzone-2); k++)
    for(int j=ghostzone-1; j<dirlength-(ghostzone-2); j++)
      for(int i=ghostzone-1; i<dirlength-(ghostzone-2); i++) {
        const int indm2  = indexf(dirlength, i-2*xdir, j-2*ydir, k-2*zdir);
        const int indm1  = indexf(dirlength, i-xdir,   j-ydir,   k-zdir);
        const int index  = indexf(dirlength, i,        j ,       k);
        const int indp1  = indexf(dirlength, i+xdir,   j+ydir,   k+zdir);

        face_var[index] = COMPUTE_FCVAL(cell_var[indm2],
                                        cell_var[indm1],
                                        cell_var[index],
                                        cell_var[indp1]);
  }
}

int main(int argc, char **argv) {

  // Set up test data
  FILE* infile = fopen("ET_Legacy_flux_source_input.bin", "rb");
  check_file_was_successfully_open(infile, "ET_Legacy_flux_source_input.bin");

  int dirlength;
  int key = fread(&dirlength, sizeof(int), 1, infile);
  const int ghostzone = 3;
  const int arraylength = dirlength*dirlength*dirlength;
  double invdx = 1.0/0.1;

  const double poison = 0.0/0.0;

  eos_parameters eos;
  eos.compute_h_and_cs2 = &compute_h_and_cs2;  

  // Allocate memory for metric
  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);
  double *gxx   = (double*) malloc(sizeof(double)*arraylength);
  double *gxy   = (double*) malloc(sizeof(double)*arraylength);
  double *gxz   = (double*) malloc(sizeof(double)*arraylength);
  double *gyy   = (double*) malloc(sizeof(double)*arraylength);
  double *gyz   = (double*) malloc(sizeof(double)*arraylength);
  double *gzz   = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for face-interpolated metric
  double *face_lapse = (double*) malloc(sizeof(double)*arraylength);
  double *face_betax = (double*) malloc(sizeof(double)*arraylength);
  double *face_betay = (double*) malloc(sizeof(double)*arraylength);
  double *face_betaz = (double*) malloc(sizeof(double)*arraylength);
  double *face_gxx   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gxy   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gxz   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gyy   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gyz   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gzz   = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for metric quantity derivatives
  double *D_lapse[3];
  double *D_betax[3];
  double *D_betay[3];
  double *D_betaz[3];
  double *D_gxx[3];
  double *D_gxy[3];
  double *D_gxz[3];
  double *D_gyy[3];
  double *D_gyz[3];
  double *D_gzz[3];
  for(int i=0; i<3; i++) {
    D_lapse[i] = (double*) malloc(sizeof(double)*arraylength);
    D_betax[i] = (double*) malloc(sizeof(double)*arraylength);
    D_betay[i] = (double*) malloc(sizeof(double)*arraylength);
    D_betaz[i] = (double*) malloc(sizeof(double)*arraylength);
    D_gxx[i]   = (double*) malloc(sizeof(double)*arraylength);
    D_gxy[i]   = (double*) malloc(sizeof(double)*arraylength);
    D_gxz[i]   = (double*) malloc(sizeof(double)*arraylength);
    D_gyy[i]   = (double*) malloc(sizeof(double)*arraylength);
    D_gyz[i]   = (double*) malloc(sizeof(double)*arraylength);
    D_gzz[i]   = (double*) malloc(sizeof(double)*arraylength);
  }

  // Allocate memory for extrinsic curvature
  double *kxx   = (double*) malloc(sizeof(double)*arraylength);
  double *kxy   = (double*) malloc(sizeof(double)*arraylength);
  double *kxz   = (double*) malloc(sizeof(double)*arraylength);
  double *kyy   = (double*) malloc(sizeof(double)*arraylength);
  double *kyz   = (double*) malloc(sizeof(double)*arraylength);
  double *kzz   = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for primitives
  double *rho   = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx    = (double*) malloc(sizeof(double)*arraylength);
  double *vy    = (double*) malloc(sizeof(double)*arraylength);
  double *vz    = (double*) malloc(sizeof(double)*arraylength);
  double *Bx    = (double*) malloc(sizeof(double)*arraylength);
  double *By    = (double*) malloc(sizeof(double)*arraylength);
  double *Bz    = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for right face
  double *rho_r   = (double*) malloc(sizeof(double)*arraylength);
  double *press_r = (double*) malloc(sizeof(double)*arraylength);
  double *vx_r    = (double*) malloc(sizeof(double)*arraylength);
  double *vy_r    = (double*) malloc(sizeof(double)*arraylength);
  double *vz_r    = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_r    = (double*) malloc(sizeof(double)*arraylength);
  double *By_r    = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_r    = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for left face
  double *rho_l   = (double*) malloc(sizeof(double)*arraylength);
  double *press_l = (double*) malloc(sizeof(double)*arraylength);
  double *vx_l    = (double*) malloc(sizeof(double)*arraylength);
  double *vy_l    = (double*) malloc(sizeof(double)*arraylength);
  double *vz_l    = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_l    = (double*) malloc(sizeof(double)*arraylength);
  double *By_l    = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_l    = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for fluxes
  double *tau_flux = (double*) malloc(sizeof(double)*arraylength);
  double *rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_flux = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_flux = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_flux = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for right-hand sides
  double *rho_star_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *tau_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_rhs = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(gxx,     sizeof(double), arraylength, infile);
  key += fread(gxy,     sizeof(double), arraylength, infile);
  key += fread(gxz,     sizeof(double), arraylength, infile);
  key += fread(gyy,     sizeof(double), arraylength, infile);
  key += fread(gyz,     sizeof(double), arraylength, infile);
  key += fread(gzz,     sizeof(double), arraylength, infile);
  key += fread(lapse,   sizeof(double), arraylength, infile);
  key += fread(betax,   sizeof(double), arraylength, infile);
  key += fread(betay,   sizeof(double), arraylength, infile);
  key += fread(betaz,   sizeof(double), arraylength, infile);
  key += fread(kxx,     sizeof(double), arraylength, infile);
  key += fread(kxy,     sizeof(double), arraylength, infile);
  key += fread(kxz,     sizeof(double), arraylength, infile);
  key += fread(kyy,     sizeof(double), arraylength, infile);
  key += fread(kyz,     sizeof(double), arraylength, infile);
  key += fread(kzz,     sizeof(double), arraylength, infile);

  key += fread(rho,     sizeof(double), arraylength, infile);
  key += fread(press,   sizeof(double), arraylength, infile);
  key += fread(vx,      sizeof(double), arraylength, infile);
  key += fread(vy,      sizeof(double), arraylength, infile);
  key += fread(vz,      sizeof(double), arraylength, infile);
  key += fread(Bx,      sizeof(double), arraylength, infile);
  key += fread(By,      sizeof(double), arraylength, infile);
  key += fread(Bz,      sizeof(double), arraylength, infile);

  if(key != arraylength*24)
    grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Initialize rhs variables to 0 so we can safely use += operator
#pragma omp parallel for
  for(int k=ghostzone; k<dirlength-ghostzone; k++)
    for(int j=ghostzone; j<dirlength-ghostzone; j++)
      for(int i=ghostzone; i<dirlength-ghostzone; i++) {
        const int index  = indexf(dirlength, i, j ,k);
        rho_star_rhs[index] = 0.0;
        tau_rhs[index] = 0.0;
        S_x_rhs[index] = 0.0;
        S_y_rhs[index] = 0.0;
        S_z_rhs[index] = 0.0;
  }

  // Function pointer to allow for loop over fluxes
  void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict, const primitive_quantities *restrict,
                              const eos_parameters *restrict, const metric_quantities *restrict, conservative_quantities *restrict);

  // Loop over flux directions (x,y,z)
  for(int flux_dirn=0; flux_dirn<3; flux_dirn++) {
    const int xdir = (flux_dirn == 0);
    const int ydir = (flux_dirn == 1);
    const int zdir = (flux_dirn == 2);

    // Set function pointer to specific function for a given direction
    switch(flux_dirn) {
      case 0:
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn0;
        break;
      case 1:
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn1;
        break;
      case 2:
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn2;
        break;
    }

    key  = fread(rho_r,   sizeof(double), arraylength, infile);
    key += fread(press_r, sizeof(double), arraylength, infile);
    key += fread(vx_r,    sizeof(double), arraylength, infile);
    key += fread(vy_r,    sizeof(double), arraylength, infile);
    key += fread(vz_r,    sizeof(double), arraylength, infile);
    key += fread(Bx_r,    sizeof(double), arraylength, infile);
    key += fread(By_r,    sizeof(double), arraylength, infile);
    key += fread(Bz_r,    sizeof(double), arraylength, infile);

    key += fread(rho_l,   sizeof(double), arraylength, infile);
    key += fread(press_l, sizeof(double), arraylength, infile);
    key += fread(vx_l,    sizeof(double), arraylength, infile);
    key += fread(vy_l,    sizeof(double), arraylength, infile);
    key += fread(vz_l,    sizeof(double), arraylength, infile);
    key += fread(Bx_l,    sizeof(double), arraylength, infile);
    key += fread(By_l,    sizeof(double), arraylength, infile);
    key += fread(Bz_l,    sizeof(double), arraylength, infile);

    if(key != arraylength*16)
      grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                   "is up-to-date with current test version.\n");

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
#pragma omp parallel for
    for(int k=ghostzone; k<dirlength-(ghostzone-1); k++)
      for(int j=ghostzone; j<dirlength-(ghostzone-1); j++)
        for(int i=ghostzone; i<dirlength-(ghostzone-1); i++) {
          const int index  = indexf(dirlength, i, j ,k);

          metric_quantities metric_face;
          initialize_metric(face_lapse[index],
                            face_gxx[index], face_gxy[index], face_gxz[index],
                            face_gyy[index], face_gyz[index], face_gzz[index],
                            face_betax[index], face_betay[index], face_betaz[index],
                            &metric_face);

          primitive_quantities prims_r, prims_l;
          initialize_primitives(rho_r[index], press_r[index], poison,
                                vx_r[index], vy_r[index], vz_r[index],
                                Bx_r[index], By_r[index], Bz_r[index],
                                poison, poison, poison, // entropy, Y_e, temp
                                &prims_r);

          initialize_primitives(rho_l[index], press_l[index], poison,
                                vx_l[index], vy_l[index], vz_l[index],
                                Bx_l[index], By_l[index], Bz_l[index],
                                poison, poison, poison, // entropy, Y_e, temp
                                &prims_l);

          // Generate randomized u^0
          prims_r.u0 = rho_r[index]*Bx_r[index]/vy_r[index];
          prims_l.u0 = rho_l[index]*Bx_l[index]/vy_l[index];

          conservative_quantities cons_fluxes;
          calculate_HLLE_fluxes(&prims_r, 
                                &prims_l,
                                &eos,
                                &metric_face, 
                                &cons_fluxes);

          rho_star_flux[index]  = cons_fluxes.rho;
          tau_flux[index]       = cons_fluxes.tau;
          S_x_flux[index]       = cons_fluxes.S_x;
          S_y_flux[index]       = cons_fluxes.S_y;
          S_z_flux[index]       = cons_fluxes.S_z;
    }

#pragma omp parallel for
    for(int k=ghostzone; k<dirlength-ghostzone; k++)
      for(int j=ghostzone; j<dirlength-ghostzone; j++)
        for(int i=ghostzone; i<dirlength-ghostzone; i++) {
          const int index  = indexf(dirlength, i, j ,k);
          const int indp1  = indexf(dirlength, i+xdir, j+ydir, k+zdir);

          rho_star_rhs[index] += invdx*(rho_star_flux[index] - rho_star_flux[indp1]);
          tau_rhs[index] += invdx*(tau_flux[index] - tau_flux[indp1]);
          S_x_rhs[index] += invdx*(S_x_flux[index] - S_x_flux[indp1]);
          S_y_rhs[index] += invdx*(S_y_flux[index] - S_y_flux[indp1]);
          S_z_rhs[index] += invdx*(S_z_flux[index] - S_z_flux[indp1]);

          D_lapse[flux_dirn][index] = invdx*(face_lapse[indp1] - face_lapse[index]);
          D_betax[flux_dirn][index] = invdx*(face_betax[indp1] - face_betax[index]);
          D_betay[flux_dirn][index] = invdx*(face_betay[indp1] - face_betay[index]);
          D_betaz[flux_dirn][index] = invdx*(face_betaz[indp1] - face_betaz[index]);
          D_gxx[flux_dirn][index] = invdx*(face_gxx[indp1] - face_gxx[index]);
          D_gxy[flux_dirn][index] = invdx*(face_gxy[indp1] - face_gxy[index]);
          D_gxz[flux_dirn][index] = invdx*(face_gxz[indp1] - face_gxz[index]);
          D_gyy[flux_dirn][index] = invdx*(face_gyy[indp1] - face_gyy[index]);
          D_gyz[flux_dirn][index] = invdx*(face_gyz[indp1] - face_gyz[index]);
          D_gzz[flux_dirn][index] = invdx*(face_gzz[indp1] - face_gzz[index]);
    }
  }

#pragma omp parallel for
  for(int k=ghostzone; k<dirlength-ghostzone; k++)
    for(int j=ghostzone; j<dirlength-ghostzone; j++)
      for(int i=ghostzone; i<dirlength-ghostzone; i++) {
        const int index  = indexf(dirlength, i, j ,k);

        metric_quantities metric;
        initialize_metric(lapse[index],
                          gxx[index], gxy[index], gxz[index],
                          gyy[index], gyz[index], gzz[index],
                          betax[index], betay[index], betaz[index],
                          &metric);

        extrinsic_curvature curv;
        initialize_extrinsic_curvature(
                          kxx[index], kxy[index], kxz[index],
                          kyy[index], kyz[index], kzz[index],
                          &curv);

        metric_derivatives metric_derivs;
        for(int dirn=0; dirn<3; dirn++) {
          metric_derivs.lapse[dirn] = D_lapse[dirn][index];
          metric_derivs.betax[dirn] = D_betax[dirn][index];
          metric_derivs.betay[dirn] = D_betay[dirn][index];
          metric_derivs.betaz[dirn] = D_betaz[dirn][index];
          metric_derivs.adm_gxx[dirn] = D_gxx[dirn][index];
          metric_derivs.adm_gxy[dirn] = D_gxy[dirn][index];
          metric_derivs.adm_gxz[dirn] = D_gxz[dirn][index];
          metric_derivs.adm_gyy[dirn] = D_gyy[dirn][index];
          metric_derivs.adm_gyz[dirn] = D_gyz[dirn][index];
          metric_derivs.adm_gzz[dirn] = D_gzz[dirn][index];
        }

        primitive_quantities prims;
        initialize_primitives(rho[index], press[index], poison,
                              vx[index], vy[index], vz[index],
                              Bx[index], By[index], Bz[index],
                              poison, poison, poison, // entropy, Y_e, temp
                              &prims);
        prims.u0  = rho[index]*Bx[index] / vy[index];


        conservative_quantities cons_sources;
        calculate_all_source_terms(&prims,
                                   &eos,
                                   &metric,
                                   &curv,
                                   &metric_derivs,
                                   &cons_sources);
                                    
        tau_rhs[index] += cons_sources.tau;
        S_x_rhs[index] += cons_sources.S_x;
        S_y_rhs[index] += cons_sources.S_y;
        S_z_rhs[index] += cons_sources.S_z;
  }

  fclose(infile);

  // Allocate memory for comparison data
  double *trusted_rho_star_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_tau_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_x_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_y_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *trusted_S_z_rhs = (double*) malloc(sizeof(double)*arraylength);

  double *pert_rho_star_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *pert_tau_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_x_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_y_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *pert_S_z_rhs = (double*) malloc(sizeof(double)*arraylength);

  FILE *outfile = fopen("ET_Legacy_flux_source_output.bin", "rb");
  check_file_was_successfully_open(outfile, "ET_Legacy_flux_source_output.bin")

  key  = fread(trusted_rho_star_rhs, sizeof(double), arraylength, outfile);
  key += fread(trusted_tau_rhs,      sizeof(double), arraylength, outfile);
  key += fread(trusted_S_x_rhs,      sizeof(double), arraylength, outfile);
  key += fread(trusted_S_y_rhs,      sizeof(double), arraylength, outfile);
  key += fread(trusted_S_z_rhs,      sizeof(double), arraylength, outfile);

  if(key != arraylength*5)
    grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  fclose(outfile);

  outfile = fopen("ET_Legacy_flux_source_output_pert.bin", "rb");
  check_file_was_successfully_open(outfile, "ET_Legacy_flux_source_output_pert.bin");

  key  = fread(pert_rho_star_rhs, sizeof(double), arraylength, outfile);
  key += fread(pert_tau_rhs,      sizeof(double), arraylength, outfile);
  key += fread(pert_S_x_rhs,      sizeof(double), arraylength, outfile);
  key += fread(pert_S_y_rhs,      sizeof(double), arraylength, outfile);
  key += fread(pert_S_z_rhs,      sizeof(double), arraylength, outfile);

  if(key != arraylength*5)
    grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  fclose(outfile);

#pragma omp parallel for
  for(int k=ghostzone; k<dirlength-ghostzone; k++)
    for(int j=ghostzone; j<dirlength-ghostzone; j++)
      for(int i=ghostzone; i<dirlength-ghostzone; i++) {
        const int index  = indexf(dirlength, i, j ,k);

        if( validate(trusted_rho_star_rhs[index], rho_star_rhs[index], pert_rho_star_rhs[index]) )
          grhayl_error("Test unit_test_ET_Legacy_flux_source has failed for variable rho_star_rhs.\n"
                       "  rho_star_rhs trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_rho_star_rhs[index], rho_star_rhs[index], pert_rho_star_rhs[index],
                                                   relative_error(trusted_rho_star_rhs[index], rho_star_rhs[index]),
                                                   relative_error(trusted_rho_star_rhs[index], pert_rho_star_rhs[index]));
        if( validate(trusted_tau_rhs[index], tau_rhs[index], pert_tau_rhs[index]) )
          grhayl_error("Test unit_test_ET_Legacy_flux_source has failed for variable tau_rhs.\n"
                       "  tau_rhs trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_tau_rhs[index], tau_rhs[index], pert_tau_rhs[index],
                                                   relative_error(trusted_tau_rhs[index], tau_rhs[index]),
                                                   relative_error(trusted_tau_rhs[index], pert_tau_rhs[index]));
        if( validate(trusted_S_x_rhs[index], S_x_rhs[index], pert_S_x_rhs[index]) )
          grhayl_error("Test unit_test_ET_Legacy_flux_source has failed for variable S_x_rhs.\n"
                       "  S_x_rhs trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_x_rhs[index], S_x_rhs[index], pert_S_x_rhs[index],
                                                   relative_error(trusted_S_x_rhs[index], S_x_rhs[index]),
                                                   relative_error(trusted_S_x_rhs[index], pert_S_x_rhs[index]));
        if( validate(trusted_S_y_rhs[index], S_y_rhs[index], pert_S_y_rhs[index]) )
          grhayl_error("Test unit_test_ET_Legacy_flux_source has failed for variable S_y_rhs.\n"
                       "  S_y_rhs trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_y_rhs[index], S_y_rhs[index], pert_S_y_rhs[index],
                                                   relative_error(trusted_S_y_rhs[index], S_y_rhs[index]),
                                                   relative_error(trusted_S_y_rhs[index], pert_S_y_rhs[index]));
        if( validate(trusted_S_z_rhs[index], S_z_rhs[index], pert_S_z_rhs[index]) )
          grhayl_error("Test unit_test_ET_Legacy_flux_source has failed for variable S_z_rhs.\n"
                       "  S_z_rhs trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", trusted_S_z_rhs[index], S_z_rhs[index], pert_S_z_rhs[index],
                                                   relative_error(trusted_S_z_rhs[index], S_z_rhs[index]),
                                                   relative_error(trusted_S_z_rhs[index], pert_S_z_rhs[index]));
  }
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(face_gxx); free(face_gxy); free(face_gxz);
  free(face_gyy); free(face_gyz); free(face_gzz);
  free(face_lapse);
  free(face_betax); free(face_betay); free(face_betaz);
  for(int i=0; i<3; i++) {
    free(D_gxx[i]); free(D_gxy[i]); free(D_gxz[i]);
    free(D_gyy[i]); free(D_gyz[i]); free(D_gzz[i]);
    free(D_lapse[i]);
    free(D_betax[i]); free(D_betay[i]); free(D_betaz[i]);
  }
  grhayl_info("ET_Legacy flux/source test has passed!\n");
  free(kxx); free(kxy); free(kxz);
  free(kyy); free(kyz); free(kzz);
  free(rho); free(press);
  free(vx); free(vy); free(vz);
  free(Bx); free(By); free(Bz);
  free(rho_r); free(press_r);
  free(vx_r); free(vy_r); free(vz_r);
  free(Bx_r); free(By_r); free(Bz_r);
  free(rho_l); free(press_l);
  free(vx_l); free(vy_l); free(vz_l);
  free(Bx_l); free(By_l); free(Bz_l);
  free(rho_star_flux); free(tau_flux);
  free(S_x_flux); free(S_y_flux); free(S_z_flux);
  free(rho_star_rhs); free(tau_rhs);
  free(S_x_rhs); free(S_y_rhs); free(S_z_rhs);
  free(trusted_rho_star_rhs); free(trusted_tau_rhs);
  free(trusted_S_x_rhs); free(trusted_S_y_rhs); free(trusted_S_z_rhs);
  free(pert_rho_star_rhs); free(pert_tau_rhs);
  free(pert_S_x_rhs); free(pert_S_y_rhs); free(pert_S_z_rhs);
}
