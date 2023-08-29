#include "unit_tests.h"

int main(int argc, char **argv) {

  int neos = 1;
  double W_max = 10.0; //IGM default
  double rho_b_min = 1e-12;
  double rho_b_max = 1e300; //IGM default
  double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  double rho_ppoly[1] = {0.0};
  double Gamma_ppoly[1] = {2.0};
  double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  const int dirlength = 20;
  const int arraylength = dirlength*dirlength*dirlength;

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

  double dummy;
  // Initialize random data. Note that for this test,
  // we needn't worry too much with physical reasonableness.
  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        ghl_randomize_metric(
              &lapse[index], &betax[index], &betay[index], &betaz[index],
              &gxx[index], &gxy[index], &gxz[index],
              &gyy[index], &gyz[index], &gzz[index]);

        gxx[index] = 1.23; /*randf(scale);*/ \
        gxy[index] = 0.03; /*randf(scale);*/ \
        gxz[index] = 0.07; /*randf(scale);*/ \
        gyy[index] = 2.33; /*randf(scale);*/ \
        gyz[index] = 0.041; /*randf(scale);*/ \
        gzz[index] = 3.11; /*randf(scale);*/ \

        kxx[index] = randf(0,1.0);
        kxy[index] = randf(0,1.0);
        kxz[index] = randf(0,1.0);
        kyy[index] = randf(0,1.0);
        kyz[index] = randf(0,1.0);
        kzz[index] = randf(0,1.0);

        rho[index] = randf(0,1.0);
        press[index] = randf(0,1.0);

        rho_r[index] = randf(0,1.0);
        press_r[index] = randf(0,1.0);

        rho_l[index] = randf(0,1.0);
        press_l[index] = randf(0,1.0);

        ghl_randomize_primitives(
              &eos, rho[index], press[index],
              &vx[index], &vy[index], &vz[index],
              &Bx[index], &By[index], &Bz[index]);

        ghl_randomize_primitives(
              &eos, rho_r[index], press_r[index],
              &vx_r[index], &vy_r[index], &vz_r[index],
              &Bx_r[index], &By_r[index], &Bz_r[index]);

        ghl_randomize_primitives(
              &eos, rho_l[index], press_l[index],
              &vx_l[index], &vy_l[index], &vz_l[index],
              &Bx_l[index], &By_l[index], &Bz_l[index]);
  }

  FILE* outfile = fopen_with_check("ET_Legacy_flux_source_input.bin", "wb");
  fwrite(&dirlength, sizeof(int), 1, outfile);
  fwrite(gxx,     sizeof(double), arraylength, outfile);
  fwrite(gxy,     sizeof(double), arraylength, outfile);
  fwrite(gxz,     sizeof(double), arraylength, outfile);
  fwrite(gyy,     sizeof(double), arraylength, outfile);
  fwrite(gyz,     sizeof(double), arraylength, outfile);
  fwrite(gzz,     sizeof(double), arraylength, outfile);
  fwrite(lapse,   sizeof(double), arraylength, outfile);
  fwrite(betax,   sizeof(double), arraylength, outfile);
  fwrite(betay,   sizeof(double), arraylength, outfile);
  fwrite(betaz,   sizeof(double), arraylength, outfile);
  fwrite(kxx,     sizeof(double), arraylength, outfile);
  fwrite(kxy,     sizeof(double), arraylength, outfile);
  fwrite(kxz,     sizeof(double), arraylength, outfile);
  fwrite(kyy,     sizeof(double), arraylength, outfile);
  fwrite(kyz,     sizeof(double), arraylength, outfile);
  fwrite(kzz,     sizeof(double), arraylength, outfile);

  fwrite(rho,     sizeof(double), arraylength, outfile);
  fwrite(press,   sizeof(double), arraylength, outfile);
  fwrite(vx,      sizeof(double), arraylength, outfile);
  fwrite(vy,      sizeof(double), arraylength, outfile);
  fwrite(vz,      sizeof(double), arraylength, outfile);
  fwrite(Bx,      sizeof(double), arraylength, outfile);
  fwrite(By,      sizeof(double), arraylength, outfile);
  fwrite(Bz,      sizeof(double), arraylength, outfile);

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

  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        rho[index]     *= 1 + randf(-1.0,1.0)*1.0e-14;
        press[index]   *= 1 + randf(-1.0,1.0)*1.0e-14;
        vx[index]      *= 1 + randf(-1.0,1.0)*1.0e-14;
        vy[index]      *= 1 + randf(-1.0,1.0)*1.0e-14;
        vz[index]      *= 1 + randf(-1.0,1.0)*1.0e-14;
        Bx[index]      *= 1 + randf(-1.0,1.0)*1.0e-14;
        By[index]      *= 1 + randf(-1.0,1.0)*1.0e-14;
        Bz[index]      *= 1 + randf(-1.0,1.0)*1.0e-14;

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

  FILE* outpert = fopen_with_check("ET_Legacy_flux_source_input_pert.bin", "wb");
  fwrite(gxx,     sizeof(double), arraylength, outpert);
  fwrite(gxy,     sizeof(double), arraylength, outpert);
  fwrite(gxz,     sizeof(double), arraylength, outpert);
  fwrite(gyy,     sizeof(double), arraylength, outpert);
  fwrite(gyz,     sizeof(double), arraylength, outpert);
  fwrite(gzz,     sizeof(double), arraylength, outpert);
  fwrite(lapse,   sizeof(double), arraylength, outpert);
  fwrite(betax,   sizeof(double), arraylength, outpert);
  fwrite(betay,   sizeof(double), arraylength, outpert);
  fwrite(betaz,   sizeof(double), arraylength, outpert);
  fwrite(kxx,     sizeof(double), arraylength, outpert);
  fwrite(kxy,     sizeof(double), arraylength, outpert);
  fwrite(kxz,     sizeof(double), arraylength, outpert);
  fwrite(kyy,     sizeof(double), arraylength, outpert);
  fwrite(kyz,     sizeof(double), arraylength, outpert);
  fwrite(kzz,     sizeof(double), arraylength, outpert);

  fwrite(rho,     sizeof(double), arraylength, outpert);
  fwrite(press,   sizeof(double), arraylength, outpert);
  fwrite(vx,      sizeof(double), arraylength, outpert);
  fwrite(vy,      sizeof(double), arraylength, outpert);
  fwrite(vz,      sizeof(double), arraylength, outpert);
  fwrite(Bx,      sizeof(double), arraylength, outpert);
  fwrite(By,      sizeof(double), arraylength, outpert);
  fwrite(Bz,      sizeof(double), arraylength, outpert);

  fwrite(rho_r,   sizeof(double), arraylength, outpert);
  fwrite(press_r, sizeof(double), arraylength, outpert);
  fwrite(vx_r,    sizeof(double), arraylength, outpert);
  fwrite(vy_r,    sizeof(double), arraylength, outpert);
  fwrite(vz_r,    sizeof(double), arraylength, outpert);
  fwrite(Bx_r,    sizeof(double), arraylength, outpert);
  fwrite(By_r,    sizeof(double), arraylength, outpert);
  fwrite(Bz_r,    sizeof(double), arraylength, outpert);

  fwrite(rho_l,   sizeof(double), arraylength, outpert);
  fwrite(press_l, sizeof(double), arraylength, outpert);
  fwrite(vx_l,    sizeof(double), arraylength, outpert);
  fwrite(vy_l,    sizeof(double), arraylength, outpert);
  fwrite(vz_l,    sizeof(double), arraylength, outpert);
  fwrite(Bx_l,    sizeof(double), arraylength, outpert);
  fwrite(By_l,    sizeof(double), arraylength, outpert);
  fwrite(Bz_l,    sizeof(double), arraylength, outpert);

  // Set r/l faces for next direction
  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        rho_r[index] = randf(0,1.0);
        press_r[index] = randf(0,1.0);

        rho_l[index] = randf(0,1.0);
        press_l[index] = randf(0,1.0);

        ghl_randomize_primitives(
              &eos, rho_r[index], press_r[index],
              &vx_r[index], &vy_r[index], &vz_r[index],
              &Bx_r[index], &By_r[index], &Bz_r[index]);

        ghl_randomize_primitives(
              &eos, rho_l[index], press_l[index],
              &vx_l[index], &vy_l[index], &vz_l[index],
              &Bx_l[index], &By_l[index], &Bz_l[index]);
  }

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

  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

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

  fwrite(rho_r,   sizeof(double), arraylength, outpert);
  fwrite(press_r, sizeof(double), arraylength, outpert);
  fwrite(vx_r,    sizeof(double), arraylength, outpert);
  fwrite(vy_r,    sizeof(double), arraylength, outpert);
  fwrite(vz_r,    sizeof(double), arraylength, outpert);
  fwrite(Bx_r,    sizeof(double), arraylength, outpert);
  fwrite(By_r,    sizeof(double), arraylength, outpert);
  fwrite(Bz_r,    sizeof(double), arraylength, outpert);

  fwrite(rho_l,   sizeof(double), arraylength, outpert);
  fwrite(press_l, sizeof(double), arraylength, outpert);
  fwrite(vx_l,    sizeof(double), arraylength, outpert);
  fwrite(vy_l,    sizeof(double), arraylength, outpert);
  fwrite(vz_l,    sizeof(double), arraylength, outpert);
  fwrite(Bx_l,    sizeof(double), arraylength, outpert);
  fwrite(By_l,    sizeof(double), arraylength, outpert);
  fwrite(Bz_l,    sizeof(double), arraylength, outpert);

  // Set r/l faces for next direction
  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        rho_r[index] = randf(0,1.0);
        press_r[index] = randf(0,1.0);

        rho_l[index] = randf(0,1.0);
        press_l[index] = randf(0,1.0);

        ghl_randomize_primitives(
              &eos, rho_r[index], press_r[index],
              &vx_r[index], &vy_r[index], &vz_r[index],
              &Bx_r[index], &By_r[index], &Bz_r[index]);

        ghl_randomize_primitives(
              &eos, rho_l[index], press_l[index],
              &vx_l[index], &vy_l[index], &vz_l[index],
              &Bx_l[index], &By_l[index], &Bz_l[index]);
  }

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

  for(int k=0; k<dirlength; k++)
    for(int j=0; j<dirlength; j++)
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

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

  fwrite(rho_r,   sizeof(double), arraylength, outpert);
  fwrite(press_r, sizeof(double), arraylength, outpert);
  fwrite(vx_r,    sizeof(double), arraylength, outpert);
  fwrite(vy_r,    sizeof(double), arraylength, outpert);
  fwrite(vz_r,    sizeof(double), arraylength, outpert);
  fwrite(Bx_r,    sizeof(double), arraylength, outpert);
  fwrite(By_r,    sizeof(double), arraylength, outpert);
  fwrite(Bz_r,    sizeof(double), arraylength, outpert);

  fwrite(rho_l,   sizeof(double), arraylength, outpert);
  fwrite(press_l, sizeof(double), arraylength, outpert);
  fwrite(vx_l,    sizeof(double), arraylength, outpert);
  fwrite(vy_l,    sizeof(double), arraylength, outpert);
  fwrite(vz_l,    sizeof(double), arraylength, outpert);
  fwrite(Bx_l,    sizeof(double), arraylength, outpert);
  fwrite(By_l,    sizeof(double), arraylength, outpert);
  fwrite(Bz_l,    sizeof(double), arraylength, outpert);

  fclose(outfile);
  fclose(outpert);
}
