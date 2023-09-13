#include "unit_tests.h"

int main(int argc, char **argv) {

  const double poison = 1e300;

  const char tablepath[] = "LS220_234r_136t_50y_analmu_20091212_SVNr26.h5";
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Ye_min    = 5e-2;
  const double Ye_max    = 5e-1;
  const double T_min     = 1e-2;
  const double T_max     = 2e2;

  ghl_parameters params;
  params.max_lorenz_factor = W_max;
  params.inv_sq_max_lorenz_factor = 1.0/SQR(W_max);

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
        tablepath,
        rho_b_min, rho_b_min, rho_b_max,
        Ye_min, Ye_min, Ye_max,
        T_min, T_min, T_max,
        &eos);

  ghl_compute_h_and_cs2 = &ghl_test_compute_h_and_cs2;

  const int arraylength = 10000;

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

  // Allocate memory for right face
  double *rho_r = (double*) malloc(sizeof(double)*arraylength);
  double *Ye_r  = (double*) malloc(sizeof(double)*arraylength);
  double *T_r   = (double*) malloc(sizeof(double)*arraylength);
  double *vx_r  = (double*) malloc(sizeof(double)*arraylength);
  double *vy_r  = (double*) malloc(sizeof(double)*arraylength);
  double *vz_r  = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_r  = (double*) malloc(sizeof(double)*arraylength);
  double *By_r  = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_r  = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for left face
  double *rho_l = (double*) malloc(sizeof(double)*arraylength);
  double *Ye_l  = (double*) malloc(sizeof(double)*arraylength);
  double *T_l   = (double*) malloc(sizeof(double)*arraylength);
  double *vx_l  = (double*) malloc(sizeof(double)*arraylength);
  double *vy_l  = (double*) malloc(sizeof(double)*arraylength);
  double *vz_l  = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_l  = (double*) malloc(sizeof(double)*arraylength);
  double *By_l  = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_l  = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for characteristic speeds
  double *cxmin = (double*) malloc(sizeof(double)*arraylength);
  double *cxmax = (double*) malloc(sizeof(double)*arraylength);
  double *cymin = (double*) malloc(sizeof(double)*arraylength);
  double *cymax = (double*) malloc(sizeof(double)*arraylength);
  double *czmin = (double*) malloc(sizeof(double)*arraylength);
  double *czmax = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for fluxes
  double *rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *Y_e_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *tau_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *ent_flux      = (double*) malloc(sizeof(double)*arraylength);

  // Initialize random data
  for(int index=0; index<arraylength; index++) {
    ghl_randomize_metric(
          &lapse[index], &betax[index], &betay[index], &betaz[index],
          &gxx[index], &gxy[index], &gxz[index],
          &gyy[index], &gyz[index], &gzz[index]);

    lapse[index] = randf(1e-8, 1);
    betax[index] = randf(0, 1e-5);
    betay[index] = randf(0, 1e-5);
    betaz[index] = randf(0, 1e-5);

    /* Table limits:
    rho 1.618871e-15 1.618871e-02
    Y_e 3.500000e-02 5.500000e-01
    T   1.000000e-02 2.511886e+02
    */
    const double rhomin = log(1e-12);
    const double rhomax = log(1e-2);
    rho_r[index] = exp(randf(rhomin, rhomax));
    rho_l[index] = exp(randf(rhomin, rhomax));

    Ye_r[index]  = randf(5e-2, 5e-1);
    Ye_l[index]  = randf(5e-2, 5e-1);

    const double Tmin = log(2e-2);
    const double Tmax = log(1e2);
    T_r[index]   = exp(randf(Tmin, Tmax));
    T_l[index]   = exp(randf(Tmin, Tmax));

    double press;
    ghl_tabulated_compute_P_from_T(
          &eos, rho_r[index], Ye_r[index], T_r[index],
          &press);

    ghl_randomize_primitives(
          &eos, rho_r[index], press,
          &vx_r[index], &vy_r[index], &vz_r[index],
          &Bx_r[index], &By_r[index], &Bz_r[index]);

    ghl_tabulated_compute_P_from_T(
          &eos, rho_l[index], Ye_l[index], T_l[index],
          &press);

    ghl_randomize_primitives(
          &eos, rho_l[index], press,
          &vx_l[index], &vy_l[index], &vz_l[index],
          &Bx_l[index], &By_l[index], &Bz_l[index]);

    ghl_metric_quantities ADM_metric;
    ghl_initialize_metric(
          lapse[index],
          betax[index], betay[index], betaz[index],
          gxx[index], gxy[index], gxz[index],
          gyy[index], gyz[index], gzz[index],
          &ADM_metric);

    ghl_primitive_quantities prims_r, prims_l;
    ghl_initialize_primitives(
          rho_r[index], poison, poison,
          vx_r[index], vy_r[index], vz_r[index],
          Bx_r[index], By_r[index], Bz_r[index],
          poison, Ye_r[index], T_r[index],
          &prims_r);
    ghl_tabulated_compute_P_eps_S_from_T(
          &eos,
          prims_r.rho, prims_r.Y_e, prims_r.temperature,
          &prims_r.press, &prims_r.eps, &prims_r.entropy);

    ghl_initialize_primitives(
          rho_l[index], poison, poison,
          vx_l[index], vy_l[index], vz_l[index],
          Bx_l[index], By_l[index], Bz_l[index],
          poison, Ye_l[index], T_l[index],
          &prims_l);
    ghl_tabulated_compute_P_eps_S_from_T(
          &eos,
          prims_l.rho, prims_l.Y_e, prims_l.temperature,
          &prims_l.press, &prims_l.eps, &prims_l.entropy);

    int speed_limit __attribute__((unused)) = ghl_limit_v_and_compute_u0(
          &params, &ADM_metric, &prims_r);
    speed_limit = ghl_limit_v_and_compute_u0(
          &params, &ADM_metric, &prims_l);

    ghl_calculate_characteristic_speed_dirn0(
          &prims_r, &prims_l, &eos,
          &ADM_metric, &cxmin[index], &cxmax[index]);

    ghl_calculate_characteristic_speed_dirn1(
          &prims_r, &prims_l, &eos,
          &ADM_metric, &cymin[index], &cymax[index]);

    ghl_calculate_characteristic_speed_dirn2(
          &prims_r, &prims_l, &eos,
          &ADM_metric, &czmin[index], &czmax[index]);
  }

  char filename[100];
  for(int perturb=0; perturb<2; perturb++) {
    if(perturb) {
      sprintf(filename,"tabulated_flux_input_pert.bin");
      for(int index=0; index<arraylength; index++) {
        lapse[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        betax[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        betay[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        betaz[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        gxx[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;
        gxy[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;
        gxz[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;
        gyy[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;
        gyz[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;
        gzz[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;

        rho_r[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        Ye_r[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        T_r[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;
        vx_r[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        vy_r[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        vz_r[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        Bx_r[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        By_r[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        Bz_r[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;

        rho_l[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        Ye_l[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        T_l[index]   *= 1 + randf(-1.0,1.0)*1.0e-12;
        vx_l[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        vy_l[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        vz_l[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        Bx_l[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        By_l[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;
        Bz_l[index]  *= 1 + randf(-1.0,1.0)*1.0e-12;

        cxmin[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        cxmax[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        cymin[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        cymax[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        czmin[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
        czmax[index] *= 1 + randf(-1.0,1.0)*1.0e-12;
      }
    } else {
      sprintf(filename,"tabulated_flux_input.bin");
      FILE* infile = fopen_with_check(filename, "wb");

      fwrite(&arraylength, sizeof(int), 1, infile);
      fwrite(lapse, sizeof(double), arraylength, infile);
      fwrite(betax, sizeof(double), arraylength, infile);
      fwrite(betay, sizeof(double), arraylength, infile);
      fwrite(betaz, sizeof(double), arraylength, infile);
      fwrite(gxx,   sizeof(double), arraylength, infile);
      fwrite(gxy,   sizeof(double), arraylength, infile);
      fwrite(gxz,   sizeof(double), arraylength, infile);
      fwrite(gyy,   sizeof(double), arraylength, infile);
      fwrite(gyz,   sizeof(double), arraylength, infile);
      fwrite(gzz,   sizeof(double), arraylength, infile);

      fwrite(rho_r, sizeof(double), arraylength, infile);
      fwrite(Ye_r,  sizeof(double), arraylength, infile);
      fwrite(T_r,   sizeof(double), arraylength, infile);
      fwrite(vx_r,  sizeof(double), arraylength, infile);
      fwrite(vy_r,  sizeof(double), arraylength, infile);
      fwrite(vz_r,  sizeof(double), arraylength, infile);
      fwrite(Bx_r,  sizeof(double), arraylength, infile);
      fwrite(By_r,  sizeof(double), arraylength, infile);
      fwrite(Bz_r,  sizeof(double), arraylength, infile);

      fwrite(rho_l, sizeof(double), arraylength, infile);
      fwrite(Ye_l,  sizeof(double), arraylength, infile);
      fwrite(T_l,   sizeof(double), arraylength, infile);
      fwrite(vx_l,  sizeof(double), arraylength, infile);
      fwrite(vy_l,  sizeof(double), arraylength, infile);
      fwrite(vz_l,  sizeof(double), arraylength, infile);
      fwrite(Bx_l,  sizeof(double), arraylength, infile);
      fwrite(By_l,  sizeof(double), arraylength, infile);
      fwrite(Bz_l,  sizeof(double), arraylength, infile);

      fwrite(cxmin, sizeof(double), arraylength, infile);
      fwrite(cxmax, sizeof(double), arraylength, infile);
      fwrite(cymin, sizeof(double), arraylength, infile);
      fwrite(cymax, sizeof(double), arraylength, infile);
      fwrite(czmin, sizeof(double), arraylength, infile);
      fwrite(czmax, sizeof(double), arraylength, infile);
      fclose(infile);
    }

    sprintf(filename,"tabulated_flux_output.bin");
    if(perturb)
      sprintf(filename,"tabulated_flux_output_pert.bin");

    FILE *outfile = fopen_with_check(filename, "wb");

    double *cmin;
    double *cmax;
    for(int entropy=0; entropy<2; entropy++) {
      for(int flux_dir=0; flux_dir<3; flux_dir++) {
        void (*calculate_HLLE_fluxes)(
              ghl_primitive_quantities *restrict,
              ghl_primitive_quantities *restrict,
              const ghl_eos_parameters *restrict,
              const ghl_metric_quantities *restrict,
              const double,
              const double,
              ghl_conservative_quantities *restrict);

        switch(flux_dir) {
          case 0:
            cmin = cxmin;
            cmax = cxmax;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy : &ghl_calculate_HLLE_fluxes_dirn0_tabulated;
            break;
          case 1:
            cmin = cymin;
            cmax = cymax;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy : &ghl_calculate_HLLE_fluxes_dirn1_tabulated;
            break;
          case 2:
            cmin = czmin;
            cmax = czmax;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy : &ghl_calculate_HLLE_fluxes_dirn2_tabulated;
            break;
        }

        for(int index=0; index<arraylength; index++) {
          ghl_metric_quantities ADM_metric;
          ghl_initialize_metric(
                lapse[index],
                betax[index], betay[index], betaz[index],
                gxx[index], gxy[index], gxz[index],
                gyy[index], gyz[index], gzz[index],
                &ADM_metric);

          ghl_primitive_quantities prims_r, prims_l;
          ghl_initialize_primitives(
                rho_r[index], poison, poison,
                vx_r[index], vy_r[index], vz_r[index],
                Bx_r[index], By_r[index], Bz_r[index],
                poison, Ye_r[index], T_r[index],
                &prims_r);
          ghl_tabulated_compute_P_eps_S_from_T(
                &eos,
                prims_r.rho, prims_r.Y_e, prims_r.temperature,
                &prims_r.press, &prims_r.eps, &prims_r.entropy);

          ghl_initialize_primitives(
                rho_l[index], poison, poison,
                vx_l[index], vy_l[index], vz_l[index],
                Bx_l[index], By_l[index], Bz_l[index],
                poison, Ye_l[index], T_l[index],
                &prims_l);
          ghl_tabulated_compute_P_eps_S_from_T(
                &eos,
                prims_l.rho, prims_l.Y_e, prims_l.temperature,
                &prims_l.press, &prims_l.eps, &prims_l.entropy);

          int speed_limit __attribute__((unused)) = ghl_limit_v_and_compute_u0(
                &params, &ADM_metric, &prims_r);
          speed_limit = ghl_limit_v_and_compute_u0(
                &params, &ADM_metric, &prims_l);

          ghl_conservative_quantities cons_fluxes;
          calculate_HLLE_fluxes(
                &prims_r, &prims_l, &eos,
                &ADM_metric, cmin[index], cmax[index],
                &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          Y_e_flux[index]      = cons_fluxes.Y_e;
          tau_flux[index]      = cons_fluxes.tau;
          S_x_flux[index]      = cons_fluxes.SD[0];
          S_y_flux[index]      = cons_fluxes.SD[1];
          S_z_flux[index]      = cons_fluxes.SD[2];
          ent_flux[index]      = cons_fluxes.entropy;
        } // arraylength
        fwrite(rho_star_flux, sizeof(double), arraylength, outfile);
        fwrite(Y_e_flux,      sizeof(double), arraylength, outfile);
        fwrite(tau_flux,      sizeof(double), arraylength, outfile);
        fwrite(S_x_flux,      sizeof(double), arraylength, outfile);
        fwrite(S_y_flux,      sizeof(double), arraylength, outfile);
        fwrite(S_z_flux,      sizeof(double), arraylength, outfile);
        if(entropy)
          fwrite(ent_flux,    sizeof(double), arraylength, outfile);
      } //flux_dir
    } //entropy
    fclose(outfile);
  } //perturb
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(rho_r); free(Ye_r); free(T_r);
  free(vx_r); free(vy_r); free(vz_r);
  free(Bx_r); free(By_r); free(Bz_r);
  free(rho_l); free(Ye_l); free(T_l);
  free(vx_l); free(vy_l); free(vz_l);
  free(Bx_l); free(By_l); free(Bz_l);
  free(cxmin); free(cymin); free(czmin);
  free(cxmax); free(cymax); free(czmax);
  free(rho_star_flux); free(tau_flux); free(Y_e_flux);
  free(S_x_flux); free(S_y_flux); free(S_z_flux);
  free(ent_flux);
}
