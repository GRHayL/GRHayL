#include "cctk.h"
#include "con2prim_header.h"

// This subroutine calculates the eigenvalues of a real, symmetric 3x3
// matrix M={{M11,M12,M13},{M12,M22,M23},{M13,M23,M33}} based on the
// algorithm described in
// http://en.wikipedia.org/wiki/Eigenvalue_algorithm#Eigenvalues_of_3.C3.973_matrices
// which simply solve the cubic equation Det( M - lamnda I)=0 analytically.
// The eigenvalues are stored in lam1, lam2 and lam3.
inline void eigenvalues_3by3_real_sym_matrix(double *restrict  lam1, double *restrict  lam2, double *restrict  lam3,
                                      const double M11, const double M12, const double M13,
                                      const double M22, const double M23, const double M33) {
  double m = (M11 + M22 + M33)/3.0;
  double K11 = M11 - m, K12 = M12, K13 = M13, K22 = M22-m, K23 = M23, K33=M33-m;
  double q = 0.5* (K11*K22*K33 + K12*K23*K13 + K13*K12*K23 - K13*K22*K13
                      - K12*K12*K33 - K11*K23*K23);
  double p = ( SQR(K11) + SQR(K22) + SQR(K33) + 2.0*(SQR(K12) + SQR(K13) + SQR(K23) ) )/6.0;

  double phi;
  double p32 = sqrt(p*p*p);
  if (fabs(q) >= fabs(p32) ) {
    phi = 0.0;
  } else {
    phi = acos(q/p32)/3.0;
  }
  if (phi<0.0) phi += M_PI/3.0;

  double sqrtp = sqrt(p);
  double sqrtp_cosphi = sqrtp*cos(phi);
  double sqrtp_sqrt3_sinphi = sqrtp*sqrt(3.0)*sin(phi);
  *lam1 = m + 2.0*sqrtp_cosphi;
  *lam2 = m - sqrtp_cosphi - sqrtp_sqrt3_sinphi;
  *lam3 = m - sqrtp_cosphi + sqrtp_sqrt3_sinphi;
}

int apply_tau_floor(const GRMHD_parameters *restrict params, const eos_parameters *restrict eos,
                           metric_quantities *restrict metric, const primitive_quantities *restrict prims,
                           conservative_quantities *restrict cons, con2prim_diagnostics *restrict diagnostics) {

  double lam1,lam2,lam3;
  eigenvalues_3by3_real_sym_matrix(&lam1, &lam2, &lam3, metric->bssn_gxx, metric->bssn_gxy, metric->bssn_gxz, metric->bssn_gyy, metric->bssn_gyz, metric->bssn_gzz);
  if (lam1 < 0.0 || lam2 < 0.0 || lam3 < 0.0) {
    // Metric is not positive-defitive, reset the metric to be conformally-flat.
    metric->bssn_gxx = 1.0;
    metric->bssn_gxy = 0.0;
    metric->bssn_gxz = 0.0;
    metric->bssn_gyy = 1.0;
    metric->bssn_gyz = 0.0;
    metric->bssn_gzz = 1.0;
    metric->bssn_gupxx = metric->psi4inv;
    metric->bssn_gupxy = 0.0;
    metric->bssn_gupxz = 0.0;
    metric->bssn_gupyy = metric->psi4inv;
    metric->bssn_gupyz = 0.0;
    metric->bssn_gupzz = metric->psi4inv;
  }

  //Next, prepare for the tau and stilde fixes:

  double Bxbar = prims->Bx*ONE_OVER_SQRT_4PI;
  double Bybar = prims->By*ONE_OVER_SQRT_4PI;
  double Bzbar = prims->Bz*ONE_OVER_SQRT_4PI;

  double Bbar_x = metric->adm_gxx*Bxbar + metric->adm_gxy*Bybar + metric->adm_gxz*Bzbar;
  double Bbar_y = metric->adm_gxy*Bxbar + metric->adm_gyy*Bybar + metric->adm_gyz*Bzbar;
  double Bbar_z = metric->adm_gxz*Bxbar + metric->adm_gyz*Bybar + metric->adm_gzz*Bzbar;
  double Bbar2 = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  double Bbar = sqrt(Bbar2);

  double check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
    // need to compute Bbar specially to prevent floating-point underflow
    double Bmax = fabs(Bxbar);
    if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
    if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
    double Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
    double B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
    Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }

  double BbardotS = Bxbar*cons->S_x + Bybar*cons->S_y + Bzbar*cons->S_z;
  double hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;

  // Limit hatBbardotS
  //double max_gammav = 100.0;
  //double rhob_max = CONSERVS[RHOSTAR]/METRIC_LAP_PSI4[PSI6];
  //double hmax = 1.0 + 2.0*rhob_max;
  //double abs_hatBbardotS_max = sqrt(SQR(max_gammav)-1.0)*CONSERVS[RHOSTAR]*hmax;
  //if (fabs(hatBbardotS) > abs_hatBbardotS_max) {
  //   double fac_reduce = abs_hatBbardotS_max/fabs(hatBbardotS);
  //   double hatBbardotS_max = hatBbardotS*fac_reduce;
  //   double Bbar_inv = 1.0/Bbar;
  //   double hat_Bbar_x = Bbar_x*Bbar_inv;
  //   double hat_Bbar_y = Bbar_y*Bbar_inv;
  //   double hat_Bbar_z = Bbar_z*Bbar_inv;
  //   double sub_fact = hatBbardotS_max - hatBbardotS;
  //   CONSERVS[STILDEX] += sub_fact*hat_Bbar_x;
  //   CONSERVS[STILDEY] += sub_fact*hat_Bbar_y;
  //   CONSERVS[STILDEZ] += sub_fact*hat_Bbar_z;
  //   hatBbardotS = hatBbardotS_max;
  //   BbardotS *= fac_reduce;
  //   CONSERVS[STILDEX] = CONSERVS[STILDEX]; CONSERVS[STILDEY] = CONSERVS[STILDEY]; CONSERVS[STILDEZ] = CONSERVS[STILDEZ];
  //}


  double sdots= metric->adm_gupxx*SQR(cons->S_x)+metric->adm_gupyy*SQR(cons->S_y)+metric->adm_gupzz*SQR(cons->S_z)+2.0*
    (metric->adm_gupxy*cons->S_x*cons->S_y+metric->adm_gupxz*cons->S_x*cons->S_z+metric->adm_gupyz*cons->S_y*cons->S_z);

  double Wm = sqrt(SQR(hatBbardotS)+ SQR(cons->rho))/metric->psi6;
  double Sm2 = (SQR(Wm)*sdots + SQR(BbardotS)*(Bbar2+2.0*Wm))/SQR(Wm+Bbar2);
  double Wmin = sqrt(Sm2 + SQR(cons->rho))/metric->psi6;
  double sdots_fluid_max = sdots;

  //tau fix, applicable when B==0 and B!=0:
  if(cons->tau < 0.5*metric->psi6*Bbar2) {
    cons->tau = eos->tau_atm+0.5*metric->psi6*Bbar2;
    diagnostics->failure_checker+=1000000;
  }

  double tau_fluid_min = cons->tau - 0.5*metric->psi6*Bbar2 - (Bbar2*sdots - SQR(BbardotS))*0.5/(metric->psi6*SQR(Wmin+Bbar2));

  //Apply Stilde fix when B==0.

  if(check_B_small*check_B_small < eos->press_atm*1e-32) {
    double rhot=cons->tau*(cons->tau+2.0*cons->rho);
    double safetyfactor = 0.999999;
    //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) safetyfactor=0.99;

    if(sdots > safetyfactor*rhot) {
      double rfactm1 = sqrt((safetyfactor*rhot)/sdots);
      cons->S_x*=rfactm1;
      cons->S_y*=rfactm1;
      cons->S_z*=rfactm1;
      diagnostics->failure_checker+=10000000;
    }
  } else if(metric->psi6>params->psi6threshold) {
    //Apply new Stilde fix.
    if (tau_fluid_min < eos->tau_atm*1.001) {
      tau_fluid_min = eos->tau_atm*1.001;
      cons->tau = tau_fluid_min + 0.5*metric->psi6*Bbar2 + (Bbar2*sdots - SQR(BbardotS))*0.5/(metric->psi6*SQR(Wmin+Bbar2));
    }

    double LHS = tau_fluid_min*(tau_fluid_min+2.0*cons->rho);
    double RHS = sdots_fluid_max;

    double safetyfactor = 0.999999;
    if(safetyfactor*LHS < RHS) {
      double rfactm1 = sqrt((safetyfactor*LHS)/RHS);
      cons->S_x*=rfactm1;
      cons->S_y*=rfactm1;
      cons->S_z*=rfactm1;
      diagnostics->failure_checker+=100000000;
    }
  }

  return 0;
}
