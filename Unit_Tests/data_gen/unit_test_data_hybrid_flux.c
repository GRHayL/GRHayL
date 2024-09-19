#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  const int arraylength = 10000;

  const double poison = 1e300;

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  ghl_parameters params;
  params.max_Lorentz_factor = W_max;
  params.inv_sq_max_Lorentz_factor = 1.0/SQR(W_max);

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);

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

  // Allocate memory for characteristic speeds
  double *cxmin = (double*) malloc(sizeof(double)*arraylength);
  double *cxmax = (double*) malloc(sizeof(double)*arraylength);
  double *cymin = (double*) malloc(sizeof(double)*arraylength);
  double *cymax = (double*) malloc(sizeof(double)*arraylength);
  double *czmin = (double*) malloc(sizeof(double)*arraylength);
  double *czmax = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for fluxes
  double *rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *tau_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *ent_flux      = (double*) malloc(sizeof(double)*arraylength);

  // Initialize random data. Note that for this test,
  // we needn't worry too much with physical reasonableness.
  for(int index=0; index<arraylength; index++) {
    ghl_randomize_metric(
          &lapse[index], &betax[index], &betay[index], &betaz[index],
          &gxx[index], &gxy[index], &gxz[index],
          &gyy[index], &gyz[index], &gzz[index]);

    lapse[index] = randf(1e-8, 1);
    betax[index] = randf(0, 1e-5);
    betay[index] = randf(0, 1e-5);
    betaz[index] = randf(0, 1e-5);

    rho_r[index] = pow(10, randf(-12, -2));
    double P_cold = 0.0;
    ghl_hybrid_compute_P_cold(&eos, rho_r[index], &P_cold);
    double P_min = log(0.9*P_cold);
    double P_max = log(1000*P_cold);
    press_r[index] = exp(randf(P_min, P_max));

    rho_l[index] = pow(10, randf(-12, -2));
    ghl_hybrid_compute_P_cold(&eos, rho_l[index], &P_cold);
    P_min = log(0.9*P_cold);
    P_max = log(1000*P_cold);
    press_l[index] = exp(randf(P_min, P_max));

    ghl_randomize_primitives(
          &eos, rho_r[index], press_r[index],
          &vx_r[index], &vy_r[index], &vz_r[index],
          &Bx_r[index], &By_r[index], &Bz_r[index]);

    ghl_randomize_primitives(
          &eos, rho_l[index], press_l[index],
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
          rho_r[index], press_r[index], poison,
          vx_r[index], vy_r[index], vz_r[index],
          Bx_r[index], By_r[index], Bz_r[index],
          poison, poison, poison, // entropy, Y_e, temp
          &prims_r);

    ghl_initialize_primitives(
          rho_l[index], press_l[index], poison,
          vx_l[index], vy_l[index], vz_l[index],
          Bx_l[index], By_l[index], Bz_l[index],
          poison, poison, poison, // entropy, Y_e, temp
          &prims_l);

    bool speed_limit;
    ghl_error_codes_t error = ghl_limit_v_and_compute_u0(
          &params, &ADM_metric, &prims_r, &speed_limit);
    if(error)
      ghl_read_error_codes(error);
    error = ghl_limit_v_and_compute_u0(
          &params, &ADM_metric, &prims_l, &speed_limit);
    if(error)
      ghl_read_error_codes(error);

    prims_r.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_r.rho, prims_r.press);
    prims_l.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_l.rho, prims_l.press);

    ghl_calculate_characteristic_speed(
          0, &eos, &ADM_metric,
          &prims_r, &prims_l, &cxmin[index], &cxmax[index]);
    ghl_calculate_characteristic_speed(
          1, &eos, &ADM_metric,
          &prims_r, &prims_l, &cymin[index], &cymax[index]);
    ghl_calculate_characteristic_speed(
          1,&eos, &ADM_metric,
          &prims_r, &prims_l, &czmin[index], &czmax[index]);
  }

  double *cmin;
  double *cmax;
  char filename[100];
  for(int perturb=0; perturb<2; perturb++) {
    if(perturb) {
      for(int index=0; index<arraylength; index++) {
        lapse[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        betax[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        betay[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        betaz[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        gxx[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        gxy[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        gxz[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        gyy[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        gyz[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        gzz[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;

        rho_r[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        press_r[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        vx_r[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        vy_r[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        vz_r[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        Bx_r[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        By_r[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        Bz_r[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;

        rho_l[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        press_l[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        vx_l[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        vy_l[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        vz_l[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        Bx_l[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        By_l[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;
        Bz_l[index]    *= 1 + randf(-1.0,1.0)*1.0e-14;

        cxmin[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        cxmax[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        cymin[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        cymax[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        czmin[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
        czmax[index] *= 1 + randf(-1.0,1.0)*1.0e-14;
      }
    } else {
      sprintf(filename,"hybrid_flux_input.bin");
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

      fwrite(rho_r,   sizeof(double), arraylength, infile);
      fwrite(press_r, sizeof(double), arraylength, infile);
      fwrite(vx_r,    sizeof(double), arraylength, infile);
      fwrite(vy_r,    sizeof(double), arraylength, infile);
      fwrite(vz_r,    sizeof(double), arraylength, infile);
      fwrite(Bx_r,    sizeof(double), arraylength, infile);
      fwrite(By_r,    sizeof(double), arraylength, infile);
      fwrite(Bz_r,    sizeof(double), arraylength, infile);

      fwrite(rho_l,   sizeof(double), arraylength, infile);
      fwrite(press_l, sizeof(double), arraylength, infile);
      fwrite(vx_l,    sizeof(double), arraylength, infile);
      fwrite(vy_l,    sizeof(double), arraylength, infile);
      fwrite(vz_l,    sizeof(double), arraylength, infile);
      fwrite(Bx_l,    sizeof(double), arraylength, infile);
      fwrite(By_l,    sizeof(double), arraylength, infile);
      fwrite(Bz_l,    sizeof(double), arraylength, infile);

      fwrite(cxmin, sizeof(double), arraylength, infile);
      fwrite(cxmax, sizeof(double), arraylength, infile);
      fwrite(cymin, sizeof(double), arraylength, infile);
      fwrite(cymax, sizeof(double), arraylength, infile);
      fwrite(czmin, sizeof(double), arraylength, infile);
      fwrite(czmax, sizeof(double), arraylength, infile);
      fclose(infile);
    }

    sprintf(filename,"hybrid_flux_output.bin");
    if(perturb)
      sprintf(filename,"hybrid_flux_output_pert.bin");
    FILE *outfile = fopen_with_check(filename, "wb");

    void (*calculate_HLLE_fluxes)(
          ghl_primitive_quantities *restrict,
          ghl_primitive_quantities *restrict,
          const ghl_eos_parameters *restrict,
          const ghl_metric_quantities *restrict,
          const double,
          const double,
          ghl_conservative_quantities *restrict);

    for(int entropy=0; entropy<2; entropy++) {
      for(int flux_dir=0; flux_dir<3; flux_dir++) {
        switch(flux_dir) {
          case 0:
            cmin = cxmin;
            cmax = cxmax;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn0_hybrid;
            break;
          case 1:
            cmin = cymin;
            cmax = cymax;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn1_hybrid;
            break;
          case 2:
            cmin = czmin;
            cmax = czmax;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn2_hybrid;
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
                rho_r[index], press_r[index], poison,
                vx_r[index], vy_r[index], vz_r[index],
                Bx_r[index], By_r[index], Bz_r[index],
                poison, poison, poison, // entropy, Y_e, temp
                &prims_r);

          ghl_initialize_primitives(
                rho_l[index], press_l[index], poison,
                vx_l[index], vy_l[index], vz_l[index],
                Bx_l[index], By_l[index], Bz_l[index],
                poison, poison, poison, // entropy, Y_e, temp
                &prims_l);

          bool speed_limit;
          ghl_error_codes_t error = ghl_limit_v_and_compute_u0(
                &params, &ADM_metric, &prims_r, &speed_limit);
          if(error)
            ghl_read_error_codes(error);
          error = ghl_limit_v_and_compute_u0(
                &params, &ADM_metric, &prims_l, &speed_limit);
          if(error)
            ghl_read_error_codes(error);

          prims_r.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_r.rho, prims_r.press);
          prims_l.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_l.rho, prims_l.press);

          ghl_conservative_quantities cons_fluxes;
          calculate_HLLE_fluxes(
                &prims_r, &prims_l, &eos,
                &ADM_metric, cmin[index], cmax[index],
                &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux[index]      = cons_fluxes.tau;
          S_x_flux[index]      = cons_fluxes.SD[0];
          S_y_flux[index]      = cons_fluxes.SD[1];
          S_z_flux[index]      = cons_fluxes.SD[2];
          ent_flux[index]      = cons_fluxes.entropy;
        } // arraylength
        fwrite(rho_star_flux, sizeof(double), arraylength, outfile);
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
  free(rho_r); free(press_r);
  free(vx_r); free(vy_r); free(vz_r);
  free(Bx_r); free(By_r); free(Bz_r);
  free(rho_l); free(press_l);
  free(vx_l); free(vy_l); free(vz_l);
  free(Bx_l); free(By_l); free(Bz_l);
  free(cxmin); free(cymin); free(czmin);
  free(cxmax); free(cymax); free(czmax);
  free(rho_star_flux); free(tau_flux); free(ent_flux);
  free(S_x_flux); free(S_y_flux); free(S_z_flux);
}
