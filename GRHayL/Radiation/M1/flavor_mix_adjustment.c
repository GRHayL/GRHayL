#include "ghl_radiation.h"

void flavor_mix_adjustment(double rN_e,
                           double rN_ae,
                           double rN_x,
                           double rN_ax,
                           const double Nmin,
                           double rE_e,
                           double rE_ae,
                           double rE_x,
                           double rE_ax,
                           ghl_radiation_flux_vector rF_e,
                           ghl_radiation_flux_vector rF_ae,
                           ghl_radiation_flux_vector rF_x,
                           ghl_radiation_flux_vector rF_ax,
                           const bool mix_type){

  //neutrino number and number density
  const double ne = max(Nmin, rN_e);
  const double nae = max(Nmin, rN_ae);
  const double nx = max(Nmin, rN_x);
  const double nax = max(Nmin, rN_ax);

  const double N = ne + nae + nx + nax;
  const double Ne = ne - nae;
  const double Nx = nx - nax;

  //calculate maximally mixed state density
  double ne_mm, nae_mm, nx_mm, nax_mm;
  if (mix_type == 0){
    ne_mm = -N/6 + Ne/2 + sqrt(4*N*N + 12*Ne*Ne - 3*Nx*Nx)/6;
    nae_mm = ne_mm - Ne;
    nx_mm = 0.5*(N + Ne + Nx) - ne_mm;
    nax_mm = nx_mm - Nx;
  }
  else {
    ne_mm = N/6 + Ne/2;
    nae_mm = N/6 - Ne/2;
    nx_mm = N/3 + Nx/2;
    nax_mm = N/3 - Nx/2;
  }

  //declare a parameter to determine what neutrinos participate in mixing
  double alpha_nu = 0; //This is NOT the lapse, just to be clear...

  if(ne_mm < 0){
    alpha = max(alpha, -ne_mm/(ne - ne_mm));
  }
  if(nae_mm < 0){
    alpha = max(alpha, -nae_mm/(nae - nae_mm));
  }
  if(nx_mm < 0){
    alpha = max(alpha, -nx_mm/(nx - nx_mm));
  }
  if(nax_mm < 0){
    alpha = max(alpha, -nax_mm/(nax - nax_mm));
  }
  alpha = min(1.0, alpha);

  //Now find the mixing densities
  const double ne_mix = max(Nmin, alpha*ne + (1-alpha)*ne_mm);
  const double nae_mix = max(Nmin, alpha*nae + (1-alpha)*nae_mm);
  const double nx_mix = max(Nmin, alpha*nx + (1-alpha)*nx_mm);
  const double nax_mix = max(Nmin, alpha*nax + (1-alpha)*nax_mm);

  //Use these to find the transition probabilities
  const double Prob_ee = min(1., ne_mix/ne);
  const double Prob_ex = 1. - Prob_ee;
  const double Prob_xx = min(1., nx_mix/nx);
  const double Prob_xe = 1. - Prob_xx;
  const double Prob_a_ee = min(1., nae_mix/nae);
  const double Prob_a_ex = 1. - Prob_a_ee;
  const double Prob_a_xx = min(1., nax_mix/nax);
  const double Prob_a_xe = 1. - Prob_a_xx;

  //Now adjust Number, energy, and fluxes based on what neutrinos changed flavor.
  rN_e = ne_mix;
  rN_ae = nae_mix;
  rN_x = nx_mix;
  rN_ax = nax_mix;

  rE_e = rE_e*(Prob_ee - Prob_ex) + rE_x*Prob_xe;
  rE_ae = rE_ae*(Prob_a_ee - Prob_a_ex) + rE_ax*Prob_a_xe;
  rE_x = rE_e*Prob_ex + rE_x*(Prob_xx - Prob_xe);
  rE_ax = rE_ae*Prob_a_ex + rE_ax*(Prob_a_xx - Prob_a_xe);

  for(int i; i < 4, ++i){
    rF_e[i] = rF_e[i]*(Prob_ee - Prob_ex) + rF_x[i]*Prob_xe;
    rF_ae[i] = rF_ae[i]*(Prob_a_ee - Prob_a_ex) + rF_ax[i]*Prob_a_xe;
    rF_x[i] = rF_e[i]*Prob_ex + rF_x[i]*(Prob_xx - Prob_xe);
    rF_ax[i] = rF_ae[i]*Prob_a_ex + rF_ax[i]*(Prob_a_xx - Prob_a_xe);
  }
}
