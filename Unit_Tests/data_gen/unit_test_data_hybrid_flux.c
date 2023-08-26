#include "unit_tests.h"

int main(int argc, char **argv) {

  const double poison = 1e300;

  int neos = 1;
  double W_max = 10.0; //IGM default
  double rho_b_min = 1e-12;
  double rho_b_max = 1e300; //IGM default
  double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  double rho_ppoly[1] = {0.0};
  double Gamma_ppoly[1] = {2.0};
  double k_ppoly0 = 1.0;

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

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

  double *rho_star_flux = (double*) malloc(sizeof(double)*arraylength);
  double *tau_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_flux      = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_flux      = (double*) malloc(sizeof(double)*arraylength);

  double dummy;
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

    rho_r[index] = randf(1e-12, 1e-2);
    double P_cold = 0.0;
    ghl_hybrid_compute_P_cold(&eos, rho_r[index], &P_cold);
    press_r[index] = randf(0.9*P_cold, 1000*P_cold);

    rho_l[index] = randf(1e-12,1.0);
    ghl_hybrid_compute_P_cold(&eos, rho_l[index], &P_cold);
    press_l[index] = randf(0.9*P_cold, 1000*P_cold);

    ghl_randomize_primitives(
          &eos, rho_r[index], press_r[index], &dummy, // no need for eps
          &vx_r[index], &vy_r[index], &vz_r[index],
          &Bx_r[index], &By_r[index], &Bz_r[index]);

    ghl_randomize_primitives(
          &eos, rho_l[index], press_l[index], &dummy, // no need for eps
          &vx_l[index], &vy_l[index], &vz_l[index],
          &Bx_l[index], &By_l[index], &Bz_l[index]);

    ghl_metric_quantities metric;
    ghl_initialize_metric(
          lapse[index],
          betax[index], betay[index], betaz[index],
          gxx[index], gxy[index], gxz[index],
          gyy[index], gyz[index], gzz[index],
          &metric);

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
  }

  for(int perturb=0; perturb<2; perturb++) {
    char filename[100];
    sprintf(filename,"hybrid_flux_input.bin");
    if(perturb) {
      sprintf(filename,"hybrid_flux_input_pert.bin");
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
      }
    }
    FILE* outfile = fopen_with_check(filename, "wb");

    fwrite(&arraylength, sizeof(int), 1, outfile);
    fwrite(lapse,   sizeof(double), arraylength, outfile);
    fwrite(betax,   sizeof(double), arraylength, outfile);
    fwrite(betay,   sizeof(double), arraylength, outfile);
    fwrite(betaz,   sizeof(double), arraylength, outfile);
    fwrite(gxx,     sizeof(double), arraylength, outfile);
    fwrite(gxy,     sizeof(double), arraylength, outfile);
    fwrite(gxz,     sizeof(double), arraylength, outfile);
    fwrite(gyy,     sizeof(double), arraylength, outfile);
    fwrite(gyz,     sizeof(double), arraylength, outfile);
    fwrite(gzz,     sizeof(double), arraylength, outfile);

    fwrite(rho_r,   sizeof(double), arraylength, outfile);
    fwrite(press_r, sizeof(double), arraylength, outfile);
    fwrite(vx_r,    sizeof(double), arraylength, outfile);
    fwrite(vy_r,    sizeof(double), arraylength, outfile);
    fwrite(vz_r,    sizeof(double), arraylength, outfile);
    fwrite(Bx_r,    sizeof(double), arraylength, outfile);
    fwrite(By_r,    sizeof(double), arraylength, outfile);
    fwrite(Bz_r,    sizeof(double), arraylength, outfile);

    fwrite(rho_l,   sizeof(double), arraylength, outfile);
    fwrite(press_l, sizeof(double), arraylength, outfile);
    fwrite(vx_l,    sizeof(double), arraylength, outfile);
    fwrite(vy_l,    sizeof(double), arraylength, outfile);
    fwrite(vz_l,    sizeof(double), arraylength, outfile);
    fwrite(Bx_l,    sizeof(double), arraylength, outfile);
    fwrite(By_l,    sizeof(double), arraylength, outfile);
    fwrite(Bz_l,    sizeof(double), arraylength, outfile);
    fclose(outfile);

    sprintf(filename,"hybrid_flux_output.bin");
    if(perturb)
      sprintf(filename,"hybrid_flux_output_pert.bin");
    outfile = fopen_with_check(filename, "wb");

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

        void (*calculate_characteristic_speed)(
              ghl_primitive_quantities *restrict,
              ghl_primitive_quantities *restrict,
              const ghl_eos_parameters *restrict,
              const ghl_metric_quantities *restrict,
              double *restrict,
              double *restrict);
        switch(flux_dir) {
          case 0:
            calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn0_hybrid;
            break;
          case 1:
            calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
            calculate_HLLE_fluxes = (entropy) ? &ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy : &ghl_calculate_HLLE_fluxes_dirn1_hybrid;
            break;
          case 2:
            calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
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

          int speed_limit __attribute__((unused)) = ghl_limit_v_and_compute_u0(
                &eos, &ADM_metric, &prims_r);
          speed_limit = ghl_limit_v_and_compute_u0(
                &eos, &ADM_metric, &prims_l);

          prims_r.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_r.rho, prims_r.press);
          prims_l.entropy = ghl_hybrid_compute_entropy_function(&eos, prims_l.rho, prims_l.press);

          double cmin, cmax;
          calculate_characteristic_speed(
                &prims_r, &prims_l, &eos,
                &ADM_metric, &cmin, &cmax);

          ghl_conservative_quantities cons_fluxes;
          calculate_HLLE_fluxes(
                &prims_r, &prims_l, &eos,
                &ADM_metric, cmin, cmax,
                &cons_fluxes);

          rho_star_flux[index] = cons_fluxes.rho;
          tau_flux[index]      = cons_fluxes.tau;
          S_x_flux[index]      = cons_fluxes.SD[0];
          S_y_flux[index]      = cons_fluxes.SD[1];
          S_z_flux[index]      = cons_fluxes.SD[2];
        } // arraylength
        fwrite(rho_star_flux, sizeof(double), arraylength, outfile);
        fwrite(tau_flux,      sizeof(double), arraylength, outfile);
        fwrite(S_x_flux,      sizeof(double), arraylength, outfile);
        fwrite(S_y_flux,      sizeof(double), arraylength, outfile);
        fwrite(S_z_flux,      sizeof(double), arraylength, outfile);
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
  free(rho_star_flux); free(tau_flux);
  free(S_x_flux); free(S_y_flux); free(S_z_flux);
}
