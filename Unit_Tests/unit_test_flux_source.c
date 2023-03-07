#include "unit_tests.h"
#include "flux_source_unit_test.h"

static inline void compute_h_and_cs2(struct eos_parameters const *restrict eos,
                                     primitive_quantities const *restrict prims,
                                     double *restrict h,
                                     double *restrict cs2) {


// CCTK_REAL h_ = U[PRESSURE]*U[VX]/U[VZ];
*h = prims->press*prims->vx / prims->vz;

// CCTK_REAL c_s_squared  = U[RHOB]*U[VZ]*(h)/1e4;
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

  for(int k=ghostzone-1; k<dirlength-(ghostzone-2); k++)
    for(int j=ghostzone-1; j<dirlength-(ghostzone-2); j++)
      for(int i=ghostzone-1; i<dirlength-(ghostzone-2); i++) {
        const int indm2  = indexf(dirlength, i-2*xdir, j-2*ydir, k-2*zdir);
        const int indm1  = indexf(dirlength, i-xdir,   j-ydir,   k-zdir);
        const int index  = indexf(dirlength, i,        j ,       k);
        const int indp1  = indexf(dirlength, i+xdir,   j+ydir,   k+zdir);

    face_var[index] = COMPUTE_FCVAL(cell_var[indm2], cell_var[indm1], cell_var[index], cell_var[indp1]);
  }
}

/*
 * // main() function:
 * // Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
 * // Step 1: Write test data to gridfunctions
 * // Step 2: Overwrite all data in ghost zones with NaNs
 * // Step 3: Apply curvilinear boundary conditions
 * // Step 4: Print gridfunction data after curvilinear boundary conditions have been applied
 * // Step 5: Free all allocated memory
 */
int main(int argc, char **argv) {

  eos_parameters eos;
  eos.compute_h_and_cs2 = &compute_h_and_cs2;  
  const double poison = 0.0/0.0;

  const int ghostzone = 3; 
  const int dirlength = Nxx_plus_2NGHOSTS0;
  const int arraylength = dirlength*dirlength*dirlength;
 
  // Step 0.m: Allocate memory for non y_n_gfs. We do this here to free up
  //         memory for setting up initial data (for cases in which initial
  //         data setup is memory intensive.)
  double *restrict     evol_gfs    = (double *restrict)malloc(sizeof(double) * NUM_EVOL_GFS * arraylength);
  double *restrict etk_evol_gfs    = (double *restrict)malloc(sizeof(double) * NUM_EVOL_GFS * arraylength);
  double *restrict auxevol_gfs = (double *restrict)malloc(sizeof(double) * NUM_AUXEVOL_GFS * arraylength);

  // Allocate memory for metric
  double *gxx   = (double*) malloc(sizeof(double)*arraylength);
  double *gxy   = (double*) malloc(sizeof(double)*arraylength);
  double *gxz   = (double*) malloc(sizeof(double)*arraylength);
  double *gyy   = (double*) malloc(sizeof(double)*arraylength);
  double *gyz   = (double*) malloc(sizeof(double)*arraylength);
  double *gzz   = (double*) malloc(sizeof(double)*arraylength);
  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for face-interpolated metric
  double *face_gxx   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gxy   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gxz   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gyy   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gyz   = (double*) malloc(sizeof(double)*arraylength);
  double *face_gzz   = (double*) malloc(sizeof(double)*arraylength);
  double *face_lapse = (double*) malloc(sizeof(double)*arraylength);
  double *face_betax = (double*) malloc(sizeof(double)*arraylength);
  double *face_betay = (double*) malloc(sizeof(double)*arraylength);
  double *face_betaz = (double*) malloc(sizeof(double)*arraylength);

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
  double *tau_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *rho_star_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_rhs = (double*) malloc(sizeof(double)*arraylength);

  // We begin by generating random inital data.

  for(int which_gf=0; which_gf<NUM_AUXEVOL_GFS; which_gf++) for(int idx=0; idx<arraylength; idx++) {
    auxevol_gfs[IDX4ptS(which_gf, idx)] = poison;
        }

  for(int which_gf=0; which_gf<NUM_EVOL_GFS; which_gf++) for(int idx=0; idx<arraylength; idx++) {
    evol_gfs[IDX4ptS(which_gf, idx)] = poison;
        }

  for(int i=0; i<arraylength; i++) {
    rho_star_rhs[i] = 0.0;
    tau_rhs[i] = 0.0;
    S_x_rhs[i] = 0.0;
    S_y_rhs[i] = 0.0;
    S_z_rhs[i] = 0.0;
  }

  /*
  After reading in the same random data used in the original IGM code, now we
  calculate face cell-face values and derivatives at cell centers
  */

  void (*calculate_HLLE_fluxes)(const primitive_quantities *restrict, const primitive_quantities *restrict,
                              const eos_parameters *restrict, const metric_quantities *restrict, conservative_quantities *restrict);

  for(int flux_dirn=0; flux_dirn<3; flux_dirn++) {
    const int xdir = (flux_dirn == 0);
    const int ydir = (flux_dirn == 1);
    const int zdir = (flux_dirn == 2);

    switch(flux_dirn) {
      case 0:
        read_from_binary_file_all("flux_source_dirn1.bin", auxevol_gfs);
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn0;
        for(int i=0; i<arraylength; i++) {
          gxx[i]   = auxevol_gfs[IDX4ptS(GAMMADD00GF, i)];
          gxy[i]   = auxevol_gfs[IDX4ptS(GAMMADD01GF, i)];
          gxz[i]   = auxevol_gfs[IDX4ptS(GAMMADD02GF, i)];
          gyy[i]   = auxevol_gfs[IDX4ptS(GAMMADD11GF, i)];
          gyz[i]   = auxevol_gfs[IDX4ptS(GAMMADD12GF, i)];
          gzz[i]   = auxevol_gfs[IDX4ptS(GAMMADD22GF, i)];
          lapse[i] = auxevol_gfs[IDX4ptS(ALPHAGF, i)];
          betax[i] = auxevol_gfs[IDX4ptS(BETAU0GF, i)];
          betay[i] = auxevol_gfs[IDX4ptS(BETAU1GF, i)];
          betaz[i] = auxevol_gfs[IDX4ptS(BETAU2GF, i)];

          kxx[i] = auxevol_gfs[IDX4ptS(KDD00GF, i)];
          kxy[i] = auxevol_gfs[IDX4ptS(KDD01GF, i)];
          kxz[i] = auxevol_gfs[IDX4ptS(KDD02GF, i)];
          kyy[i] = auxevol_gfs[IDX4ptS(KDD11GF, i)];
          kyz[i] = auxevol_gfs[IDX4ptS(KDD12GF, i)];
          kzz[i] = auxevol_gfs[IDX4ptS(KDD22GF, i)];

          rho[i] = auxevol_gfs[IDX4ptS(RHOBGF, i)];
          press[i] = auxevol_gfs[IDX4ptS(PGF, i)];
          vx[i] = auxevol_gfs[IDX4ptS(VU0GF, i)];
          vy[i] = auxevol_gfs[IDX4ptS(VU1GF, i)];
          vz[i] = auxevol_gfs[IDX4ptS(VU2GF, i)];
          Bx[i] = auxevol_gfs[IDX4ptS(BU0GF, i)];
          By[i] = auxevol_gfs[IDX4ptS(BU1GF, i)];
          Bz[i] = auxevol_gfs[IDX4ptS(BU2GF, i)];
        }
        break;
      case 1:
        read_from_binary_file_recons("flux_source_dirn2.bin", auxevol_gfs);
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn1;
        break;
      case 2:
        read_from_binary_file_recons("flux_source_dirn3.bin", auxevol_gfs);
        calculate_HLLE_fluxes = &calculate_HLLE_fluxes_dirn2;
        break;
    }

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

    for(int k=ghostzone-1; k<dirlength-(ghostzone-1); k++)
      for(int j=ghostzone-1; j<dirlength-(ghostzone-1); j++)
        for(int i=ghostzone-1; i<dirlength-(ghostzone-1); i++) {
          const int index  = indexf(dirlength, i, j ,k);
          int idx  = IDX3S(i, j, k);
          rho_r[index] = auxevol_gfs[IDX4ptS(RHOB_RGF, idx)];
          press_r[index] = auxevol_gfs[IDX4ptS(P_RGF, idx)];
          vx_r[index] = auxevol_gfs[IDX4ptS(VRU0GF, idx)];
          vy_r[index] = auxevol_gfs[IDX4ptS(VRU1GF, idx)];
          vz_r[index] = auxevol_gfs[IDX4ptS(VRU2GF, idx)];
          Bx_r[index] = auxevol_gfs[IDX4ptS(BRU0GF, idx)];
          By_r[index] = auxevol_gfs[IDX4ptS(BRU1GF, idx)];
          Bz_r[index] = auxevol_gfs[IDX4ptS(BRU2GF, idx)];

          rho_l[index] = auxevol_gfs[IDX4ptS(RHOB_LGF, idx)];
          press_l[index] = auxevol_gfs[IDX4ptS(P_LGF, idx)];
          vx_l[index] = auxevol_gfs[IDX4ptS(VLU0GF, idx)];
          vy_l[index] = auxevol_gfs[IDX4ptS(VLU1GF, idx)];
          vz_l[index] = auxevol_gfs[IDX4ptS(VLU2GF, idx)];
          Bx_l[index] = auxevol_gfs[IDX4ptS(BLU0GF, idx)];
          By_l[index] = auxevol_gfs[IDX4ptS(BLU1GF, idx)];
          Bz_l[index] = auxevol_gfs[IDX4ptS(BLU2GF, idx)];
    }

    for(int k=ghostzone-1; k<dirlength-(ghostzone-1); k++)
      for(int j=ghostzone-1; j<dirlength-(ghostzone-1); j++)
        for(int i=ghostzone-1; i<dirlength-(ghostzone-1); i++) {
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

          prims_r.u0 = rho_r[index]*Bx_r[index] / vy_r[index];
          prims_l.u0 = rho_l[index]*Bx_l[index] / vy_l[index];

          conservative_quantities cons_fluxes;
          calculate_HLLE_fluxes(&prims_r, 
                                &prims_l,
                                &eos,
                                &metric_face, 
                                &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux[index] = cons_fluxes.tau;
          S_x_flux[index] = cons_fluxes.S_x;
          S_y_flux[index] = cons_fluxes.S_y;
          S_z_flux[index] = cons_fluxes.S_z;
    }

    for(int k=ghostzone-1; k<dirlength-(ghostzone-1); k++)
      for(int j=ghostzone-1; j<dirlength-(ghostzone-1); j++)
        for(int i=ghostzone-1; i<dirlength-(ghostzone-1); i++) {
          const int index  = indexf(dirlength, i, j ,k);
          const int indp1  = indexf(dirlength, i+xdir, j+ydir, k+zdir);

          rho_star_rhs[index]  += invdx*(rho_star_flux[index] - rho_star_flux[indp1]);
          tau_rhs[index]  += invdx*(tau_flux[index] - tau_flux[indp1]);
          S_x_rhs[index]  += invdx*(S_x_flux[index] - S_x_flux[indp1]);
          S_y_rhs[index]  += invdx*(S_y_flux[index] - S_y_flux[indp1]);
          S_z_rhs[index]  += invdx*(S_z_flux[index] - S_z_flux[indp1]);

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

    for(int k=ghostzone-1; k<dirlength-(ghostzone-1); k++)
      for(int j=ghostzone-1; j<dirlength-(ghostzone-1); j++)
        for(int i=ghostzone-1; i<dirlength-(ghostzone-1); i++) {
          const int index  = indexf(dirlength, i, j ,k);

          metric_quantities metric;
          initialize_metric(lapse[index],
                            gxx[index], gxy[index], gxz[index],
                            gyy[index], gyz[index], gzz[index],
                            betax[index], betay[index], betaz[index],
                            &metric);

          extrinsic_curvature curv;
          curv.Kxy = kxy[index];
          curv.Kxz = kxz[index];
          curv.Kyz = kyz[index];
          curv.Kxx = kxx[index];
          curv.Kyy = kyy[index];
          curv.Kzz = kzz[index];

          /* Consider alternative code setup
          metric_quantities metric_derivs[3];
          for(int dirn=0; dirn<3; dirn++) {
            metric_derivs[dirn].lapse = D_lapse[dirn][index];
            metric_derivs[dirn].betax = D_betax[dirn][index];
            metric_derivs[dirn].betay = D_betay[dirn][index];
            metric_derivs[dirn].betaz = D_betaz[dirn][index];
            metric_derivs[dirn].adm_gxx = D_gxx[dirn][index];
            metric_derivs[dirn].adm_gxy = D_gxy[dirn][index];
            metric_derivs[dirn].adm_gxz = D_gxz[dirn][index];
            metric_derivs[dirn].adm_gyy = D_gyy[dirn][index];
            metric_derivs[dirn].adm_gyz = D_gyz[dirn][index];
            metric_derivs[dirn].adm_gzz = D_gzz[dirn][index];
          }
          */
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
                                      
          S_x_rhs[index]  += cons_sources.S_x;
          S_y_rhs[index]  += cons_sources.S_y;
          S_z_rhs[index]  += cons_sources.S_z;
          tau_rhs[index] += cons_sources.tau;
    }

    FILE *infile = fopen("flux_source_output_rhs_data.bin", "rb");
    double rhs_correct_magic_number = 9.524300707856655e-3;
    double rhs_magic_number1, rhs_magic_number2, rhs_magic_number3;
    
    int key;
    
    key  = fread(&rhs_magic_number1, sizeof(double), 1, infile);
    key += fread(etk_evol_gfs + arraylength*RHO_STARGF, sizeof(double), arraylength, infile);
    key += fread(etk_evol_gfs + arraylength*TAU_TILDEGF, sizeof(double), arraylength, infile);
    key += fread(etk_evol_gfs + arraylength*STILDED0GF, sizeof(double), arraylength, infile);
    key += fread(&rhs_magic_number2, sizeof(double), 1, infile);
    key += fread(etk_evol_gfs + arraylength*STILDED1GF, sizeof(double), arraylength, infile);
    key += fread(etk_evol_gfs + arraylength*STILDED2GF, sizeof(double), arraylength, infile);
    key += fread(&rhs_magic_number3, sizeof(double), 1, infile);
    
    fclose(infile);
    if(rhs_magic_number1!=rhs_correct_magic_number){ printf("ERROR: magic_number1 does not match"); exit(1);}
    if(rhs_magic_number2!=rhs_correct_magic_number){ printf("ERROR: magic_number2 does not match"); exit(1);}
    if(rhs_magic_number3!=rhs_correct_magic_number){ printf("ERROR: magic_number3 does not match"); exit(1);}

    for(int k=ghostzone; k<dirlength-ghostzone-1; k++)
      for(int j=ghostzone; j<dirlength-ghostzone-1; j++)
        for(int i=ghostzone; i<dirlength-ghostzone-1; i++) {
          const int index  = indexf(dirlength, i, j ,k);
          int idx  = IDX3S(i, j, k);

      const double rho_rel_diff = log10(fabs(0.5*(rho_star_rhs[index] - etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)]) / (rho_star_rhs[index] + etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)])));
      const double tau_rel_diff = log10(fabs(0.5*(tau_rhs[index] - etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)]) / (tau_rhs[index] + etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)])));
      const double S_x_rel_diff = log10(fabs(0.5*(S_x_rhs[index] - etk_evol_gfs[IDX4ptS(STILDED0GF, idx)]) / (S_x_rhs[index] + etk_evol_gfs[IDX4ptS(STILDED0GF, idx)])));
      const double S_y_rel_diff = log10(fabs(0.5*(S_y_rhs[index] - etk_evol_gfs[IDX4ptS(STILDED1GF, idx)]) / (S_y_rhs[index] + etk_evol_gfs[IDX4ptS(STILDED1GF, idx)])));
      const double S_z_rel_diff = log10(fabs(0.5*(S_z_rhs[index] - etk_evol_gfs[IDX4ptS(STILDED2GF, idx)]) / (S_z_rhs[index] + etk_evol_gfs[IDX4ptS(STILDED2GF, idx)])));

      if(isnan(rho_star_rhs[index]) && isfinite(etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)])) {
        printf("ERROR: rho_star evaluated to nan!\n");
        exit(1);
      } else if(rho_rel_diff > -9) {
        printf("ERROR: S_x_rel_diff is too high\n");
        exit(1);
      }
      if(isnan(tau_rhs[index]) && isfinite(etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)])) {
        printf("ERROR: tau evaluated to nan!\n");
        exit(1);
      } else if(tau_rel_diff > -9) {
        printf("ERROR: tau_rel_diff is too high\n");
        exit(1);
      }
      if(isnan(S_x_rhs[index]) && isfinite(etk_evol_gfs[IDX4ptS(STILDED0GF, idx)])) {
        printf("ERROR: S_x evaluated to nan!\n");
        exit(1);
      } else if(S_x_rel_diff > -9) {
        printf("ERROR: S_x_rel_diff is too high\n");
        exit(1);
      }
      if(isnan(S_y_rhs[index]) && isfinite(etk_evol_gfs[IDX4ptS(STILDED1GF, idx)])) {
        printf("ERROR: S_y evaluated to nan!\n");
        exit(1);
      } else if(S_y_rel_diff > -9) {
        printf("ERROR: S_y_rel_diff is too high\n");
        exit(1);
      }
      if(isnan(S_z_rhs[index]) && isfinite(etk_evol_gfs[IDX4ptS(STILDED2GF, idx)])) {
        printf("ERROR: S_z evaluated to nan!\n");
        exit(1);
      } else if(S_z_rel_diff > -9) {
        printf("ERROR: S_z_rel_diff is too high\n");
        exit(1);
      }
  }

  // Step 4: Free all allocated memory
  free(evol_gfs);
  free(etk_evol_gfs);
  free(auxevol_gfs);
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
  free(kxx); free(kxy); free(kxz);
  free(kyy); free(kyz); free(kzz);
  free(rho);
  free(press);
  free(vx); free(vy); free(vz);
  free(Bx); free(By); free(Bz);
  free(rho_r);
  free(press_r);
  free(vx_r); free(vy_r); free(vz_r);
  free(Bx_r); free(By_r); free(Bz_r);
  free(rho_l);
  free(press_l);
  free(vx_l); free(vy_l); free(vz_l);
  free(Bx_l); free(By_l); free(Bz_l);
  free(rho_star_flux);
  free(tau_flux);
  free(S_x_flux); free(S_y_flux); free(S_z_flux);
  free(rho_star_rhs);
  free(tau_rhs);
  free(S_x_rhs); free(S_y_rhs); free(S_z_rhs);
  return 0;
}
