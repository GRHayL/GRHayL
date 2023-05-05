#include "Neutrinos.h"
#include "unit_tests.h"

#define IDX3D(i0,i1,i2) ( (i0) + Nt0*( (i1) + Nt1*(i2) ) )

void update_optical_depths_path_of_least_resistance(
      const int Nt0,
      const int Nt1,
      const int Nt2,
      const int Ng0,
      const int Ng1,
      const int Ng2,
      const double dxx0,
      const double dxx1,
      const double dxx2,
      const double *restrict gammaDD00,
      const double *restrict gammaDD11,
      const double *restrict gammaDD22,
      const double *restrict kappa_0_nue,
      const double *restrict kappa_1_nue,
      const double *restrict kappa_0_anue,
      const double *restrict kappa_1_anue,
      const double *restrict kappa_0_nux,
      const double *restrict kappa_1_nux,
      const double *restrict tau_0_nue_p,
      const double *restrict tau_1_nue_p,
      const double *restrict tau_0_anue_p,
      const double *restrict tau_1_anue_p,
      const double *restrict tau_0_nux_p,
      const double *restrict tau_1_nux_p,
      double *restrict tau_0_nue,
      double *restrict tau_1_nue,
      double *restrict tau_0_anue,
      double *restrict tau_1_anue,
      double *restrict tau_0_nux,
      double *restrict tau_1_nux ) {

#pragma omp parallel for
  for(int i2=Ng2;i2<Nt2-Ng2;i2++) {
    for(int i1=Ng1;i1<Nt1-Ng1;i1++) {
      for(int i0=Ng0;i0<Nt0-Ng0;i0++) {

        // Step 1: Set gridpoint indices
        const int i0_i1_i2   = IDX3D(i0  ,i1,i2  );
        const int i0p1_i1_i2 = IDX3D(i0+1,i1,i2  );
        const int i0m1_i1_i2 = IDX3D(i0-1,i1,i2  );
        const int i0_i1p1_i2 = IDX3D(i0,i1+1,i2  );
        const int i0_i1m1_i2 = IDX3D(i0,i1-1,i2  );
        const int i0_i1_i2p1 = IDX3D(i0,i1  ,i2+1);
        const int i0_i1_i2m1 = IDX3D(i0,i1  ,i2-1);

        // Step 2: Read in metric gfs from main memory
        const double gammaDD00_i0_i1_i2 = gammaDD00[i0_i1_i2];
        const double gammaDD00_i0p1_i1_i2 = gammaDD00[i0p1_i1_i2];
        const double gammaDD00_i0m1_i1_i2 = gammaDD00[i0m1_i1_i2];

        const double gammaDD11_i0_i1_i2 = gammaDD11[i0_i1_i2];
        const double gammaDD11_i0_i1p1_i2 = gammaDD11[i0_i1p1_i2];
        const double gammaDD11_i0_i1m1_i2 = gammaDD11[i0_i1m1_i2];

        const double gammaDD22_i0_i1_i2 = gammaDD22[i0_i1_i2];
        const double gammaDD22_i0_i1_i2p1 = gammaDD22[i0_i1_i2p1];
        const double gammaDD22_i0_i1_i2m1 = gammaDD22[i0_i1_i2m1];


        // Step 3: Read in opacity gfs from main memory
        const double kappa_0_nue_i0_i1_i2 = kappa_0_nue[i0_i1_i2];
        const double kappa_0_nue_i0p1_i1_i2 = kappa_0_nue[i0p1_i1_i2];
        const double kappa_0_nue_i0m1_i1_i2 = kappa_0_nue[i0m1_i1_i2];
        const double kappa_0_nue_i0_i1p1_i2 = kappa_0_nue[i0_i1p1_i2];
        const double kappa_0_nue_i0_i1m1_i2 = kappa_0_nue[i0_i1m1_i2];
        const double kappa_0_nue_i0_i1_i2p1 = kappa_0_nue[i0_i1_i2p1];
        const double kappa_0_nue_i0_i1_i2m1 = kappa_0_nue[i0_i1_i2m1];

        const double kappa_1_nue_i0_i1_i2 = kappa_1_nue[i0_i1_i2];
        const double kappa_1_nue_i0p1_i1_i2 = kappa_1_nue[i0p1_i1_i2];
        const double kappa_1_nue_i0m1_i1_i2 = kappa_1_nue[i0m1_i1_i2];
        const double kappa_1_nue_i0_i1p1_i2 = kappa_1_nue[i0_i1p1_i2];
        const double kappa_1_nue_i0_i1m1_i2 = kappa_1_nue[i0_i1m1_i2];
        const double kappa_1_nue_i0_i1_i2p1 = kappa_1_nue[i0_i1_i2p1];
        const double kappa_1_nue_i0_i1_i2m1 = kappa_1_nue[i0_i1_i2m1];

        const double kappa_0_anue_i0_i1_i2 = kappa_0_anue[i0_i1_i2];
        const double kappa_0_anue_i0p1_i1_i2 = kappa_0_anue[i0p1_i1_i2];
        const double kappa_0_anue_i0m1_i1_i2 = kappa_0_anue[i0m1_i1_i2];
        const double kappa_0_anue_i0_i1p1_i2 = kappa_0_anue[i0_i1p1_i2];
        const double kappa_0_anue_i0_i1m1_i2 = kappa_0_anue[i0_i1m1_i2];
        const double kappa_0_anue_i0_i1_i2p1 = kappa_0_anue[i0_i1_i2p1];
        const double kappa_0_anue_i0_i1_i2m1 = kappa_0_anue[i0_i1_i2m1];

        const double kappa_1_anue_i0_i1_i2 = kappa_1_anue[i0_i1_i2];
        const double kappa_1_anue_i0p1_i1_i2 = kappa_1_anue[i0p1_i1_i2];
        const double kappa_1_anue_i0m1_i1_i2 = kappa_1_anue[i0m1_i1_i2];
        const double kappa_1_anue_i0_i1p1_i2 = kappa_1_anue[i0_i1p1_i2];
        const double kappa_1_anue_i0_i1m1_i2 = kappa_1_anue[i0_i1m1_i2];
        const double kappa_1_anue_i0_i1_i2p1 = kappa_1_anue[i0_i1_i2p1];
        const double kappa_1_anue_i0_i1_i2m1 = kappa_1_anue[i0_i1_i2m1];

        const double kappa_0_nux_i0_i1_i2 = kappa_0_nux[i0_i1_i2];
        const double kappa_0_nux_i0p1_i1_i2 = kappa_0_nux[i0p1_i1_i2];
        const double kappa_0_nux_i0m1_i1_i2 = kappa_0_nux[i0m1_i1_i2];
        const double kappa_0_nux_i0_i1p1_i2 = kappa_0_nux[i0_i1p1_i2];
        const double kappa_0_nux_i0_i1m1_i2 = kappa_0_nux[i0_i1m1_i2];
        const double kappa_0_nux_i0_i1_i2p1 = kappa_0_nux[i0_i1_i2p1];
        const double kappa_0_nux_i0_i1_i2m1 = kappa_0_nux[i0_i1_i2m1];

        const double kappa_1_nux_i0_i1_i2 = kappa_1_nux[i0_i1_i2];
        const double kappa_1_nux_i0p1_i1_i2 = kappa_1_nux[i0p1_i1_i2];
        const double kappa_1_nux_i0m1_i1_i2 = kappa_1_nux[i0m1_i1_i2];
        const double kappa_1_nux_i0_i1p1_i2 = kappa_1_nux[i0_i1p1_i2];
        const double kappa_1_nux_i0_i1m1_i2 = kappa_1_nux[i0_i1m1_i2];
        const double kappa_1_nux_i0_i1_i2p1 = kappa_1_nux[i0_i1_i2p1];
        const double kappa_1_nux_i0_i1_i2m1 = kappa_1_nux[i0_i1_i2m1];


        // Step 4: Compute metric at cell faces
        const double gammaDD00_i0phalf_i1_i2 = 0.5*(gammaDD00_i0_i1_i2 + gammaDD00_i0p1_i1_i2);
        const double gammaDD00_i0mhalf_i1_i2 = 0.5*(gammaDD00_i0_i1_i2 + gammaDD00_i0m1_i1_i2);

        const double gammaDD11_i0_i1phalf_i2 = 0.5*(gammaDD11_i0_i1_i2 + gammaDD11_i0_i1p1_i2);
        const double gammaDD11_i0_i1mhalf_i2 = 0.5*(gammaDD11_i0_i1_i2 + gammaDD11_i0_i1m1_i2);

        const double gammaDD22_i0_i1_i2phalf = 0.5*(gammaDD22_i0_i1_i2 + gammaDD22_i0_i1_i2p1);
        const double gammaDD22_i0_i1_i2mhalf = 0.5*(gammaDD22_i0_i1_i2 + gammaDD22_i0_i1_i2m1);


        // Step 5: Compute ds^{i} = sqrt(gamma_{ii}dx^{i}dx^{i})
        const double ds_i0phalf_i1_i2 = sqrt(dxx0*dxx0*gammaDD00_i0phalf_i1_i2);
        const double ds_i0mhalf_i1_i2 = sqrt(dxx0*dxx0*gammaDD00_i0mhalf_i1_i2);
        const double ds_i0_i1phalf_i2 = sqrt(dxx1*dxx1*gammaDD11_i0_i1phalf_i2);
        const double ds_i0_i1mhalf_i2 = sqrt(dxx1*dxx1*gammaDD11_i0_i1mhalf_i2);
        const double ds_i0_i1_i2phalf = sqrt(dxx2*dxx2*gammaDD22_i0_i1_i2phalf);
        const double ds_i0_i1_i2mhalf = sqrt(dxx2*dxx2*gammaDD22_i0_i1_i2mhalf);

        // Step 6: Compute opacities at cell faces
        const double kappa_0_nue_i0phalf_i1_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0p1_i1_i2);
        const double kappa_0_nue_i0mhalf_i1_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0m1_i1_i2);
        const double kappa_0_nue_i0_i1phalf_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1p1_i2);
        const double kappa_0_nue_i0_i1mhalf_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1m1_i2);
        const double kappa_0_nue_i0_i1_i2phalf = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1_i2p1);
        const double kappa_0_nue_i0_i1_i2mhalf = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1_i2m1);

        const double kappa_1_nue_i0phalf_i1_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0p1_i1_i2);
        const double kappa_1_nue_i0mhalf_i1_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0m1_i1_i2);
        const double kappa_1_nue_i0_i1phalf_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1p1_i2);
        const double kappa_1_nue_i0_i1mhalf_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1m1_i2);
        const double kappa_1_nue_i0_i1_i2phalf = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1_i2p1);
        const double kappa_1_nue_i0_i1_i2mhalf = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1_i2m1);

        const double kappa_0_anue_i0phalf_i1_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0p1_i1_i2);
        const double kappa_0_anue_i0mhalf_i1_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0m1_i1_i2);
        const double kappa_0_anue_i0_i1phalf_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1p1_i2);
        const double kappa_0_anue_i0_i1mhalf_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1m1_i2);
        const double kappa_0_anue_i0_i1_i2phalf = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1_i2p1);
        const double kappa_0_anue_i0_i1_i2mhalf = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1_i2m1);

        const double kappa_1_anue_i0phalf_i1_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0p1_i1_i2);
        const double kappa_1_anue_i0mhalf_i1_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0m1_i1_i2);
        const double kappa_1_anue_i0_i1phalf_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1p1_i2);
        const double kappa_1_anue_i0_i1mhalf_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1m1_i2);
        const double kappa_1_anue_i0_i1_i2phalf = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1_i2p1);
        const double kappa_1_anue_i0_i1_i2mhalf = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1_i2m1);

        const double kappa_0_nux_i0phalf_i1_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0p1_i1_i2);
        const double kappa_0_nux_i0mhalf_i1_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0m1_i1_i2);
        const double kappa_0_nux_i0_i1phalf_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1p1_i2);
        const double kappa_0_nux_i0_i1mhalf_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1m1_i2);
        const double kappa_0_nux_i0_i1_i2phalf = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1_i2p1);
        const double kappa_0_nux_i0_i1_i2mhalf = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1_i2m1);

        const double kappa_1_nux_i0phalf_i1_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0p1_i1_i2);
        const double kappa_1_nux_i0mhalf_i1_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0m1_i1_i2);
        const double kappa_1_nux_i0_i1phalf_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1p1_i2);
        const double kappa_1_nux_i0_i1mhalf_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1m1_i2);
        const double kappa_1_nux_i0_i1_i2phalf = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1_i2p1);
        const double kappa_1_nux_i0_i1_i2mhalf = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1_i2m1);


        // Step 7: Compute optical depth at neighboring points
        const double tau_0_nue_i0p1_i1_i2 = tau_0_nue_p[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_0_nue_i0phalf_i1_i2;
        const double tau_0_nue_i0m1_i1_i2 = tau_0_nue_p[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_0_nue_i0mhalf_i1_i2;
        const double tau_0_nue_i0_i1p1_i2 = tau_0_nue_p[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_0_nue_i0_i1phalf_i2;
        const double tau_0_nue_i0_i1m1_i2 = tau_0_nue_p[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_0_nue_i0_i1mhalf_i2;
        const double tau_0_nue_i0_i1_i2p1 = tau_0_nue_p[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_0_nue_i0_i1_i2phalf;
        const double tau_0_nue_i0_i1_i2m1 = tau_0_nue_p[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_0_nue_i0_i1_i2mhalf;

        const double tau_1_nue_i0p1_i1_i2 = tau_1_nue_p[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_1_nue_i0phalf_i1_i2;
        const double tau_1_nue_i0m1_i1_i2 = tau_1_nue_p[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_1_nue_i0mhalf_i1_i2;
        const double tau_1_nue_i0_i1p1_i2 = tau_1_nue_p[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_1_nue_i0_i1phalf_i2;
        const double tau_1_nue_i0_i1m1_i2 = tau_1_nue_p[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_1_nue_i0_i1mhalf_i2;
        const double tau_1_nue_i0_i1_i2p1 = tau_1_nue_p[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_1_nue_i0_i1_i2phalf;
        const double tau_1_nue_i0_i1_i2m1 = tau_1_nue_p[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_1_nue_i0_i1_i2mhalf;

        const double tau_0_anue_i0p1_i1_i2 = tau_0_anue_p[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_0_anue_i0phalf_i1_i2;
        const double tau_0_anue_i0m1_i1_i2 = tau_0_anue_p[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_0_anue_i0mhalf_i1_i2;
        const double tau_0_anue_i0_i1p1_i2 = tau_0_anue_p[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_0_anue_i0_i1phalf_i2;
        const double tau_0_anue_i0_i1m1_i2 = tau_0_anue_p[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_0_anue_i0_i1mhalf_i2;
        const double tau_0_anue_i0_i1_i2p1 = tau_0_anue_p[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_0_anue_i0_i1_i2phalf;
        const double tau_0_anue_i0_i1_i2m1 = tau_0_anue_p[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_0_anue_i0_i1_i2mhalf;

        const double tau_1_anue_i0p1_i1_i2 = tau_1_anue_p[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_1_anue_i0phalf_i1_i2;
        const double tau_1_anue_i0m1_i1_i2 = tau_1_anue_p[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_1_anue_i0mhalf_i1_i2;
        const double tau_1_anue_i0_i1p1_i2 = tau_1_anue_p[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_1_anue_i0_i1phalf_i2;
        const double tau_1_anue_i0_i1m1_i2 = tau_1_anue_p[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_1_anue_i0_i1mhalf_i2;
        const double tau_1_anue_i0_i1_i2p1 = tau_1_anue_p[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_1_anue_i0_i1_i2phalf;
        const double tau_1_anue_i0_i1_i2m1 = tau_1_anue_p[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_1_anue_i0_i1_i2mhalf;

        const double tau_0_nux_i0p1_i1_i2 = tau_0_nux_p[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_0_nux_i0phalf_i1_i2;
        const double tau_0_nux_i0m1_i1_i2 = tau_0_nux_p[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_0_nux_i0mhalf_i1_i2;
        const double tau_0_nux_i0_i1p1_i2 = tau_0_nux_p[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_0_nux_i0_i1phalf_i2;
        const double tau_0_nux_i0_i1m1_i2 = tau_0_nux_p[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_0_nux_i0_i1mhalf_i2;
        const double tau_0_nux_i0_i1_i2p1 = tau_0_nux_p[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_0_nux_i0_i1_i2phalf;
        const double tau_0_nux_i0_i1_i2m1 = tau_0_nux_p[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_0_nux_i0_i1_i2mhalf;

        const double tau_1_nux_i0p1_i1_i2 = tau_1_nux_p[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_1_nux_i0phalf_i1_i2;
        const double tau_1_nux_i0m1_i1_i2 = tau_1_nux_p[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_1_nux_i0mhalf_i1_i2;
        const double tau_1_nux_i0_i1p1_i2 = tau_1_nux_p[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_1_nux_i0_i1phalf_i2;
        const double tau_1_nux_i0_i1m1_i2 = tau_1_nux_p[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_1_nux_i0_i1mhalf_i2;
        const double tau_1_nux_i0_i1_i2p1 = tau_1_nux_p[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_1_nux_i0_i1_i2phalf;
        const double tau_1_nux_i0_i1_i2m1 = tau_1_nux_p[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_1_nux_i0_i1_i2mhalf;


        // Step 8: Select path of least resistance
        const double new_tau_0_nue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_0_nue_i0p1_i1_i2,tau_0_nue_i0m1_i1_i2),tau_0_nue_i0_i1p1_i2),tau_0_nue_i0_i1m1_i2),tau_0_nue_i0_i1_i2p1),tau_0_nue_i0_i1_i2m1);
        const double new_tau_1_nue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_1_nue_i0p1_i1_i2,tau_1_nue_i0m1_i1_i2),tau_1_nue_i0_i1p1_i2),tau_1_nue_i0_i1m1_i2),tau_1_nue_i0_i1_i2p1),tau_1_nue_i0_i1_i2m1);
        const double new_tau_0_anue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_0_anue_i0p1_i1_i2,tau_0_anue_i0m1_i1_i2),tau_0_anue_i0_i1p1_i2),tau_0_anue_i0_i1m1_i2),tau_0_anue_i0_i1_i2p1),tau_0_anue_i0_i1_i2m1);
        const double new_tau_1_anue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_1_anue_i0p1_i1_i2,tau_1_anue_i0m1_i1_i2),tau_1_anue_i0_i1p1_i2),tau_1_anue_i0_i1m1_i2),tau_1_anue_i0_i1_i2p1),tau_1_anue_i0_i1_i2m1);
        const double new_tau_0_nux_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_0_nux_i0p1_i1_i2,tau_0_nux_i0m1_i1_i2),tau_0_nux_i0_i1p1_i2),tau_0_nux_i0_i1m1_i2),tau_0_nux_i0_i1_i2p1),tau_0_nux_i0_i1_i2m1);
        const double new_tau_1_nux_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_1_nux_i0p1_i1_i2,tau_1_nux_i0m1_i1_i2),tau_1_nux_i0_i1p1_i2),tau_1_nux_i0_i1m1_i2),tau_1_nux_i0_i1_i2p1),tau_1_nux_i0_i1_i2m1);

         // Step 9: Write results to main memory
        tau_0_nue[i0_i1_i2] = new_tau_0_nue_i0_i1_i2;
        tau_1_nue[i0_i1_i2] = new_tau_1_nue_i0_i1_i2;
        tau_0_anue[i0_i1_i2] = new_tau_0_anue_i0_i1_i2;
        tau_1_anue[i0_i1_i2] = new_tau_1_anue_i0_i1_i2;
        tau_0_nux[i0_i1_i2] = new_tau_0_nux_i0_i1_i2;
        tau_1_nux[i0_i1_i2] = new_tau_1_nux_i0_i1_i2;

      } // for(int i0=Ng0;i0<N0-Ng0;i0++)
    } // for(int i1=Ng1;i1<N1-Ng1;i1++)
  } // for(int i2=Ng2;i2<N2-Ng2;i2++)
}

void constantdensitysphere_test(
      const eos_parameters *restrict eos,
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
  double **xx         = (double **)malloc(sizeof(double *)*3);
  xx[0]               = (double *)malloc(sizeof(double)*Nt0);
  xx[1]               = (double *)malloc(sizeof(double)*Nt1);
  xx[2]               = (double *)malloc(sizeof(double)*Nt2);
  double *gammaDD00   = (double *)malloc(sizeof(double)*Ntotal);
  double *gammaDD11   = (double *)malloc(sizeof(double)*Ntotal);
  double *gammaDD22   = (double *)malloc(sizeof(double)*Ntotal);
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
  double T_exterior   = 0.01;

  if( test_key == 1 ) {
    rho_interior *= (1+randf(-1,1)*1e-14);
    Y_e_interior *= (1+randf(-1,1)*1e-14);
    T_interior   *= (1+randf(-1,1)*1e-14);
    rho_exterior *= (1+randf(-1,1)*1e-14);
    Y_e_exterior *= (1+randf(-1,1)*1e-14);
    T_exterior   *= (1+randf(-1,1)*1e-14);
  }

  // Step 4: Compute opacities in the interior and exterior
  neutrino_optical_depths tau_in = {{0,0},{0,0},{0,0}};

  // Step 4.a: Interior
  neutrino_opacities kappa_interior;
  NRPyLeakage_compute_neutrino_opacities(eos,rho_interior,Y_e_interior,T_interior,
                                         &tau_in, &kappa_interior);

  // Step 4.b: Exterior
  neutrino_opacities kappa_exterior;
  NRPyLeakage_compute_neutrino_opacities(eos,rho_exterior,Y_e_exterior,T_exterior,
                                         &tau_in, &kappa_exterior);

  // Step 5: Print basic information
  grhayl_info("Test information:\n");
  grhayl_info("    Domain properties:\n");
  grhayl_info("        - Sphere radius = %22.15e\n",rSph);
  grhayl_info("        - xmin          = %22.15e\n",xmin);
  grhayl_info("        - xmax          = %22.15e\n",xmax);
  grhayl_info("        - ymin          = %22.15e\n",ymin);
  grhayl_info("        - ymax          = %22.15e\n",ymax);
  grhayl_info("        - zmin          = %22.15e\n",zmin);
  grhayl_info("        - zmax          = %22.15e\n",zmax);
  grhayl_info("        - Nx            = (%d) + 2x(%d)\n",N0,Ng0);
  grhayl_info("        - Ny            = (%d) + 2x(%d)\n",N1,Ng1);
  grhayl_info("        - Nz            = (%d) + 2x(%d)\n",N2,Ng2);
  grhayl_info("        - dx            = %22.15e\n",dx);
  grhayl_info("        - dy            = %22.15e\n",dy);
  grhayl_info("        - dz            = %22.15e\n",dz);
  grhayl_info("    Hydro quantities at sphere interior:\n");
  grhayl_info("        - rho           = %22.15e\n", rho_interior);
  grhayl_info("        - Y_e           = %22.15e\n", Y_e_interior);
  grhayl_info("        -  T            = %22.15e\n", T_interior);
  grhayl_info("        - kappa_0_nue   = %22.15e\n", kappa_interior.nue [0]);
  grhayl_info("        - kappa_1_nue   = %22.15e\n", kappa_interior.nue [1]);
  grhayl_info("        - kappa_0_anue  = %22.15e\n", kappa_interior.anue[0]);
  grhayl_info("        - kappa_1_anue  = %22.15e\n", kappa_interior.anue[1]);
  grhayl_info("        - kappa_0_nux   = %22.15e\n", kappa_interior.nux [0]);
  grhayl_info("        - kappa_1_nux   = %22.15e\n", kappa_interior.nux [1]);
  grhayl_info("    Hydro quantities at sphere exterior:\n");
  grhayl_info("        - rho           = %22.15e\n", rho_exterior);
  grhayl_info("        - Y_e           = %22.15e\n", Y_e_exterior);
  grhayl_info("        -  T            = %22.15e\n", T_exterior);
  grhayl_info("        - kappa_0_nue   = %22.15e\n", kappa_exterior.nue [0]);
  grhayl_info("        - kappa_1_nue   = %22.15e\n", kappa_exterior.nue [1]);
  grhayl_info("        - kappa_0_anue  = %22.15e\n", kappa_exterior.anue[0]);
  grhayl_info("        - kappa_1_anue  = %22.15e\n", kappa_exterior.anue[1]);
  grhayl_info("        - kappa_0_nux   = %22.15e\n", kappa_exterior.nux [0]);
  grhayl_info("        - kappa_1_nux   = %22.15e\n", kappa_exterior.nux [1]);

  // Step 6: Loop over the grid, set opacities
  //         and initialize optical depth to zero
#pragma omp parallel for
  for(int i2=0;i2<Nt2;i2++) {
    const double z = zmin + (i2-Ng2+0.5)*dz;
    xx[2][i2] = z;
    for(int i1=0;i1<Nt1;i1++) {
      const double y = ymin + (i1-Ng1+0.5)*dy;
      xx[1][i1] = y;
      for(int i0=0;i0<Nt0;i0++) {
        const double x = xmin + (i0-Ng0+0.5)*dx;
        xx[0][i0] = x;

        // Step 3.a: Set local index
        const int index = IDX3D(i0,i1,i2);

        // Step 3.b: Initialize metric to flat space
        gammaDD00[index] = 1.0;
        gammaDD11[index] = 1.0;
        gammaDD22[index] = 1.0;

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

    // Step 4.b: Update optical depth
    update_optical_depths_path_of_least_resistance(
                                                   Nt0,Nt1,Nt2,Ng0,Ng1,Ng2,
                                                   dx,dy,dz,
                                                   gammaDD00,gammaDD11,gammaDD22,
                                                   kappa_nue [0],kappa_nue [1],
                                                   kappa_anue[0],kappa_anue[1],
                                                   kappa_nux [0],kappa_nux [1],
                                                   tau_nue_p [0],tau_nue_p [1],
                                                   tau_anue_p[0],tau_anue_p[1],
                                                   tau_nux_p [0],tau_nux_p [1],
                                                   tau_nue   [0],tau_nue   [1],
                                                   tau_anue  [0],tau_anue  [1],
                                                   tau_nux   [0],tau_nux   [1]);
  }

  // Step 5: Dump the data
  if( test_key != 2 ) {
    FILE *fp;
    if( test_key )
      fp = fopen("nrpyleakage_constant_density_sphere_perturbed.bin", "wb");
    else
      fp = fopen("nrpyleakage_constant_density_sphere_unperturbed.bin", "wb");

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
    FILE *fp_unpert = fopen("nrpyleakage_constant_density_sphere_unperturbed.bin", "rb");
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
      grhayl_error("Error reading file nrpyleakage_constant_density_sphere_unperturbed.bin\n");

    // Read in the perturbed data
    FILE *fp_pert   = fopen("nrpyleakage_constant_density_sphere_perturbed.bin", "rb");
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
      grhayl_error("Error reading file nrpyleakage_constant_density_sphere_perturbed.bin\n");
    fclose(fp_pert);

    // Perform the validation
#pragma omp parallel for
    for(int i2=0;i2<Nt2;i2++) {
      for(int i1=0;i1<Nt1;i1++) {
        for(int i0=0;i0<Nt0;i0++) {
          const int index = IDX3D(i0,i1,i2);
          validate(kappa_nue_pert [0][index], kappa_nue [0][index], kappa_nue_unpert [0][index]);
          validate(kappa_nue_pert [1][index], kappa_nue [1][index], kappa_nue_unpert [1][index]);
          validate(kappa_anue_pert[0][index], kappa_anue[0][index], kappa_anue_unpert[0][index]);
          validate(kappa_anue_pert[1][index], kappa_anue[1][index], kappa_anue_unpert[1][index]);
          validate(kappa_nux_pert [0][index], kappa_nux [0][index], kappa_nux_unpert [0][index]);
          validate(kappa_nux_pert [1][index], kappa_nux [1][index], kappa_nux_unpert [1][index]);
          validate(tau_nue_pert   [0][index], tau_nue   [0][index], tau_nue_unpert   [0][index]);
          validate(tau_nue_pert   [1][index], tau_nue   [1][index], tau_nue_unpert   [1][index]);
          validate(tau_anue_pert  [0][index], tau_anue  [0][index], tau_anue_unpert  [0][index]);
          validate(tau_anue_pert  [1][index], tau_anue  [1][index], tau_anue_unpert  [1][index]);
          validate(tau_nux_pert   [0][index], tau_nux   [0][index], tau_nux_unpert   [0][index]);
          validate(tau_nux_pert   [1][index], tau_nux   [1][index], tau_nux_unpert   [1][index]);
        }
      }
    }
  }

  // Step 6: Free memory
  for(int i=0;i<3;i++) free(xx[i]);
  free(xx);
  free(gammaDD00);
  free(gammaDD11);
  free(gammaDD22);
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

int main(int argc, char **argv) {

  if( argc != 3 ) {
    grhayl_info("Correct usage is: %s <eos table> <test key>\n", argv[0]);
    grhayl_info("Available test keys:\n");
    grhayl_info("  0 : generate unperturbed\n");
    grhayl_info("  1 : generate perturbed\n");
    grhayl_info("  2 : perform validation test\n");
    exit(1);
  }

  const char *tablepath  = argv[1];
  const int test_key     = atoi(argv[2]);
  const double W_max     = 10.0;
  const double rho_b_atm = 1e-12;
  const double rho_b_min = -1;
  const double rho_b_max = -1;
  const double Y_e_atm   = 0.5;
  const double Y_e_min   = -1;
  const double Y_e_max   = -1;
  const double T_atm     = 1e-2;
  const double T_min     = -1;
  const double T_max     = -1;

  eos_parameters eos;
  initialize_tabulated_eos_functions_and_params(tablepath, W_max,
                                                rho_b_atm, rho_b_min, rho_b_max,
                                                Y_e_atm, Y_e_min, Y_e_max,
                                                T_atm, T_min, T_max, &eos);
  eos.root_finding_precision=1e-10;

  constantdensitysphere_test(&eos, test_key);

  eos.tabulated_free_memory(&eos);

  return 0;
}
