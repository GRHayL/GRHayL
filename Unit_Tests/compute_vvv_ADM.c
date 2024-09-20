#include "ghl_unit_tests.h"

void ghl_test_compute_vvv_ADM(
  const int dirlength,
  const double *restrict lapse,
  const double *restrict betax,
  const double *restrict betay,
  const double *restrict betaz,
  const double *restrict gxx,
  const double *restrict gxy,
  const double *restrict gxz,
  const double *restrict gyy,
  const double *restrict gyz,
  const double *restrict gzz,
  const double *restrict phitilde,
  const double *restrict Ax,
  const double *restrict Ay,
  const double *restrict Az,
  double *restrict alpha_Phi_minus_betaj_A_j_interp,
  double *restrict sqrtg_Ax_interp,
  double *restrict sqrtg_Ay_interp,
  double *restrict sqrtg_Az_interp) {

  for(int k=1; k<dirlength-1; k++) {
    for(int j=1; j<dirlength-1; j++) {
      for(int i=1; i<dirlength-1; i++) {
        const int index = indexf(dirlength,i,j,k);

        ghl_metric_quantities metric_stencil[2][2][2];
        double Ax_stencil[3][3][3];
        double Ay_stencil[3][3][3];
        double Az_stencil[3][3][3];
        induction_interp_vars interp_vars;

        // Read in variable at interpolation stencil points from main memory.
        for(int iterz=0; iterz<2; iterz++) {
          for(int itery=0; itery<2; itery++) {
            for(int iterx=0; iterx<2; iterx++) {
              const int ind = indexf(dirlength,i+iterx,j+itery,k+iterz);
              ghl_initialize_metric(
                    lapse[ind], betax[ind], betay[ind], betaz[ind],
                    gxx[ind], gxy[ind], gxz[ind],
                    gyy[ind], gyz[ind], gzz[ind],
                    &metric_stencil[iterz][itery][iterx]);
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
        ghl_interpolate_with_vertex_centered_ADM(metric_stencil, Ax_stencil, Ay_stencil, Az_stencil, phitilde[index], &interp_vars);

        sqrtg_Ax_interp[index] = interp_vars.sqrtg_Ai[0];
        sqrtg_Ay_interp[index] = interp_vars.sqrtg_Ai[1];
        sqrtg_Az_interp[index] = interp_vars.sqrtg_Ai[2];
        alpha_Phi_minus_betaj_A_j_interp[index] = interp_vars.alpha_Phi_minus_betaj_A_j;
      }
    }
  }
}
