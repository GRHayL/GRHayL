#include "unit_tests.h"

void compute_ccc_BSSN(
  const int dirlength,
  const double *restrict lapse,
  const double *restrict betax,
  const double *restrict betay,
  const double *restrict betaz,
  const double *restrict psi,
  const double *restrict gtupxx,
  const double *restrict gtupxy,
  const double *restrict gtupxz,
  const double *restrict gtupyy,
  const double *restrict gtupyz,
  const double *restrict gtupzz,
  const double *restrict phitilde,
  const double *restrict Ax,
  const double *restrict Ay,
  const double *restrict Az,
  double *restrict alpha_interp,
  double *restrict betax_interp,
  double *restrict betay_interp,
  double *restrict betaz_interp,
  double *restrict alpha_Phi_minus_betaj_A_j_interp,
  double *restrict sqrtg_Ax_interp,
  double *restrict sqrtg_Ay_interp,
  double *restrict sqrtg_Az_interp) {

  for(int k=1; k<dirlength-1; k++) {
    for(int j=1; j<dirlength-1; j++) {
      for(int i=1; i<dirlength-1; i++) {
        const int index = indexf(dirlength,i,j,k);

        metric_quantities metric_stencil[2][2][2];
        double psi_stencil[2][2][2];
        double Ax_stencil[3][3][3];
        double Ay_stencil[3][3][3];
        double Az_stencil[3][3][3];
        induction_interp_vars interp_vars;

        // Read in variable at interpolation stencil points from main memory.
        for(int iterz=0; iterz<2; iterz++) {
          for(int itery=0; itery<2; itery++) {
            for(int iterx=0; iterx<2; iterx++) {
              const int ind = indexf(dirlength,i+iterx,j+itery,k+iterz);
              metric_stencil[iterz][itery][iterx].lapse          = lapse[ind];
              metric_stencil[iterz][itery][iterx].betaU[0]       = betax[ind];
              metric_stencil[iterz][itery][iterx].betaU[1]       = betay[ind];
              metric_stencil[iterz][itery][iterx].betaU[2]       = betaz[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[0][0]  = gtupxx[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[0][1]  = gtupxy[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[0][2]  = gtupxz[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[1][1]  = gtupyy[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[1][2]  = gtupyz[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[2][2]  = gtupzz[ind];
              psi_stencil[iterz][itery][iterx] = psi[ind];
            }
          }
        }
        // A_x needs a stencil s.t. interp_limits={ 0,1,-1,1,-1,1}.
        // A_y needs a stencil s.t. interp_limits={-1,1, 0,1,-1,1}.
        // A_z needs a stencil s.t. interp_limits={-1,1,-1,1, 0,1}.
        // We could fill only the needed elements, but it is cleaner
        // to fill in the whole 3x3x3 array.
        // TODO: the old code explicitly only filled in the necessary
        // elements. If we want to remove ~15 memcopies, do that here.
        for(int iterz=-1; iterz<2; iterz++) {
          for(int itery=-1; itery<2; itery++) {
            for(int iterx=-1; iterx<2; iterx++) {
              const int ind = indexf(dirlength,i+iterx,j+itery,k+iterz);
              Ax_stencil[iterz+1][itery+1][iterx+1] = Ax[ind];
              Ay_stencil[iterz+1][itery+1][iterx+1] = Ay[ind];
              Az_stencil[iterz+1][itery+1][iterx+1] = Az[ind];
            }
          }
        }
        ghl_interpolate_with_cell_centered_BSSN(metric_stencil, psi_stencil, Ax_stencil, Ay_stencil, Az_stencil, phitilde[index], &interp_vars);

        alpha_interp[index] = interp_vars.alpha;
        betax_interp[index] = interp_vars.betai[0];
        betay_interp[index] = interp_vars.betai[1];
        betaz_interp[index] = interp_vars.betai[2];
        sqrtg_Ax_interp[index] = interp_vars.sqrtg_Ai[0];
        sqrtg_Ay_interp[index] = interp_vars.sqrtg_Ai[1];
        sqrtg_Az_interp[index] = interp_vars.sqrtg_Ai[2];
        alpha_Phi_minus_betaj_A_j_interp[index] = interp_vars.alpha_Phi_minus_betaj_A_j;
      }
    }
  }
}
