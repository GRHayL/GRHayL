#include "unit_tests.h"

inline int indexf(const int gridmax, const int i, const int j, const int k) {
  return i + j*gridmax + k*gridmax*gridmax;
}

int rel_tol(const double reltol, const double x1, const double x2) {
  const double rel_diff = relative_error(x1, x2);
  if(rel_diff > reltol) return 1;
  return 0;
}
// Tolerance limit for numerical values
const double reltol = 1.0e-15;

int main(int argc, char **argv) {
  const int gridmin     = 0;
  const int gridmax     = 21;
  const int arraylength = gridmax*gridmax*gridmax;

  const double dX[3] = {0.1, 0.1, 0.1};

  const double Lorenz_damping_factor = 0.1;

  double *gupxx = (double*) malloc(sizeof(double)*arraylength);
  double *gupxy = (double*) malloc(sizeof(double)*arraylength);
  double *gupxz = (double*) malloc(sizeof(double)*arraylength);
  double *gupyy = (double*) malloc(sizeof(double)*arraylength);
  double *gupyz = (double*) malloc(sizeof(double)*arraylength);
  double *gupzz = (double*) malloc(sizeof(double)*arraylength);

  double *psi = (double*) malloc(sizeof(double)*arraylength);
  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde = (double*) malloc(sizeof(double)*arraylength);
  double *Ax = (double*) malloc(sizeof(double)*arraylength);
  double *Ay = (double*) malloc(sizeof(double)*arraylength);
  double *Az = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_interp = (double*) malloc(sizeof(double)*arraylength);
  double *shiftx_interp = (double*) malloc(sizeof(double)*arraylength);
  double *shifty_interp = (double*) malloc(sizeof(double)*arraylength);
  double *shiftz_interp = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_Phi_minus_betaj_A_j_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Ax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Ay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Az_interp = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Az_rhs = (double*) malloc(sizeof(double)*arraylength);

//  // Read in data from file to ensure portability
  FILE* infile;
  infile = fopen("gauge_rhs_initial_data.bin", "rb");
  check_file_was_successfully_open(infile, "gauge_rhs_initial_data.bin");

  int key;
  key  = fread(gupxx, sizeof(double), arraylength, infile);
  key += fread(gupxy, sizeof(double), arraylength, infile);
  key += fread(gupxz, sizeof(double), arraylength, infile);
  key += fread(gupyy, sizeof(double), arraylength, infile);
  key += fread(gupyz, sizeof(double), arraylength, infile);
  key += fread(gupzz, sizeof(double), arraylength, infile);

  key += fread(psi,   sizeof(double), arraylength, infile);
  key += fread(lapse, sizeof(double), arraylength, infile);
  key += fread(betax, sizeof(double), arraylength, infile);
  key += fread(betay, sizeof(double), arraylength, infile);
  key += fread(betaz, sizeof(double), arraylength, infile);

  key += fread(phitilde, sizeof(double), arraylength, infile);
  key += fread(Ax,       sizeof(double), arraylength, infile);
  key += fread(Ay,       sizeof(double), arraylength, infile);
  key += fread(Az,       sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*15)
    grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Data which should be written before it is used is poisoned to validate behavior.
  // RHSs for A are set to 0 because they are assumed to already contain the zero-gauge
  // contribution to the RHS before entering this function.
  const double poison = 1e200;

#pragma omp parallel for
  for(int k=gridmin; k<gridmax; k++)
    for(int j=gridmin; j<gridmax; j++)
      for(int i=gridmin; i<gridmax; i++) {
        const int index = indexf(gridmax,i,j,k);

        alpha_interp[index]  = poison;
        shiftx_interp[index] = poison;
        shifty_interp[index] = poison;
        shiftz_interp[index] = poison;

        alpha_Phi_minus_betaj_A_j_interp[index] = poison;
        alpha_sqrtg_Ax_interp[index]            = poison;
        alpha_sqrtg_Ay_interp[index]            = poison;
        alpha_sqrtg_Az_interp[index]            = poison;

        Ax_rhs[index]       = 0.0;
        Ay_rhs[index]       = 0.0;
        Az_rhs[index]       = 0.0;
        phitilde_rhs[index] = poison;
  }

#pragma omp parallel for
  for(int k=gridmin+1; k<gridmax-1; k++)
    for(int j=gridmin+1; j<gridmax-1; j++)
      for(int i=gridmin+1; i<gridmax-1; i++) {
        const int index = indexf(gridmax,i,j,k);

        A_gauge_vars gauge_vars;
        A_gauge_rhs_vars gauge_rhs_vars;

        // Read in variable at interpolation stencil points from main memory.
        gauge_vars.phitilde = phitilde[index];
        for(int iterz=0; iterz<2; iterz++)
          for(int itery=0; itery<2; itery++)
            for(int iterx=0; iterx<2; iterx++) {
              const int ind = indexf(gridmax,i+iterx,j+itery,k+iterz);
              gauge_vars.gupxx[iterz][itery][iterx]  = gupxx[ind];
              gauge_vars.gupxy[iterz][itery][iterx]  = gupxy[ind];
              gauge_vars.gupxz[iterz][itery][iterx]  = gupxz[ind];
              gauge_vars.gupyy[iterz][itery][iterx]  = gupyy[ind];
              gauge_vars.gupyz[iterz][itery][iterx]  = gupyz[ind];
              gauge_vars.gupzz[iterz][itery][iterx]  = gupzz[ind];
              gauge_vars.lapse[iterz][itery][iterx]  = lapse[ind];
              gauge_vars.psi[iterz][itery][iterx]    = psi[ind];
              gauge_vars.shiftx[iterz][itery][iterx] = betax[ind];
              gauge_vars.shifty[iterz][itery][iterx] = betay[ind];
              gauge_vars.shiftz[iterz][itery][iterx] = betaz[ind];
        }
        // A_x needs a stencil s.t. interp_limits={ 0,1,-1,1,-1,1}.
        // A_y needs a stencil s.t. interp_limits={-1,1, 0,1,-1,1}.
        // A_z needs a stencil s.t. interp_limits={-1,1,-1,1, 0,1}.
        // We could fill only the needed elements, but it is cleaner
        // to fill in the whole 3x3x3 array.
        // TODO: the old code explicitly only filled in the necessary
        // elements. If we want to remove ~15 memcopies, do that here.
        for(int iterz=-1; iterz<2; iterz++)
          for(int itery=-1; itery<2; itery++)
            for(int iterx=-1; iterx<2; iterx++) {
              const int ind = indexf(gridmax,i+iterx,j+itery,k+iterz);
              gauge_vars.A_x[iterz+1][itery+1][iterx+1] = Ax[ind];
              gauge_vars.A_y[iterz+1][itery+1][iterx+1] = Ay[ind];
              gauge_vars.A_z[iterz+1][itery+1][iterx+1] = Az[ind];
        }
// This code should only copy the needed data that isn't copied in the loop for other variables, but it is untested.
//        for(int iter2=0; iter2<2; iter2++)
//        for(int iter1=0; iter1<2; iter1++) {
//          gauge_vars.A_x[iter2+1][0][iter1+1] = in_vars[A_XI][indexf(gridmax, i+iter1,     j-1, k+iter2)]; // { (0,1),    -1, (0,1)}
//          gauge_vars.A_x[0][iter2+1][iter1+1] = in_vars[A_XI][indexf(gridmax, i+iter1, j+iter2, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_y[iter2+1][iter1+1][0] = in_vars[A_YI][indexf(gridmax,     i-1, j+iter1, k+iter2)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_y[0][iter1+1][iter2+1] = in_vars[A_YI][indexf(gridmax, i+iter2, j+iter1, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_z[iter1+1][iter2+1][0] = in_vars[A_ZI][indexf(gridmax,     i-1, j+iter2, k+iter1)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_z[iter1+1][0][iter2+1] = in_vars[A_ZI][indexf(gridmax, i+iter2,     j-1, k+iter1)]; // { (0,1),    -1, (0,1)}
//        }

        interpolate_for_A_gauge_rhs(&gauge_vars, &gauge_rhs_vars);

        alpha_interp[index] = gauge_rhs_vars.alpha_interp;
        alpha_sqrtg_Ax_interp[index] = gauge_rhs_vars.alpha_sqrtg_Ax_interp[0];
        alpha_sqrtg_Ay_interp[index] = gauge_rhs_vars.alpha_sqrtg_Ay_interp[0];
        alpha_sqrtg_Az_interp[index] = gauge_rhs_vars.alpha_sqrtg_Az_interp[0];
        alpha_Phi_minus_betaj_A_j_interp[index] = gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[0];
        shiftx_interp[index] = gauge_rhs_vars.shiftx_interp[0];
        shifty_interp[index] = gauge_rhs_vars.shifty_interp[0];
        shiftz_interp[index] = gauge_rhs_vars.shiftz_interp[0];
      }

  const double dxinv[3] = {1.0/dX[0], 1.0/dX[1], 1.0/dX[2]};

  // This loop requires two additional ghostzones in every direction. Hence the following loop definition:
#pragma omp parallel for
  for(int k=gridmin+3; k<gridmax-3; k++)
    for(int j=gridmin+3; j<gridmax-3; j++)
      for(int i=gridmin+3; i<gridmax-3; i++) {
        const int index = indexf(gridmax,i,j,k);

        // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
        //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
        A_gauge_rhs_vars gauge_rhs_vars;
    
        gauge_rhs_vars.alpha_interp = alpha_interp[index];
    
        gauge_rhs_vars.dxi[0] = dxinv[0];
        gauge_rhs_vars.dxi[1] = dxinv[1];
        gauge_rhs_vars.dxi[2] = dxinv[2];
    
        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[0] = alpha_Phi_minus_betaj_A_j_interp[index];
        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[1] = alpha_Phi_minus_betaj_A_j_interp[indexf(gridmax,i-1,j,k)];
        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[2] = alpha_Phi_minus_betaj_A_j_interp[indexf(gridmax,i,j-1,k)];
        gauge_rhs_vars.alpha_Phi_minus_betaj_A_j_interp[3] = alpha_Phi_minus_betaj_A_j_interp[indexf(gridmax,i,j,k-1)];
    
        gauge_rhs_vars.alpha_sqrtg_Ax_interp[0] = alpha_sqrtg_Ax_interp[index];
        gauge_rhs_vars.alpha_sqrtg_Ay_interp[0] = alpha_sqrtg_Ay_interp[index];
        gauge_rhs_vars.alpha_sqrtg_Az_interp[0] = alpha_sqrtg_Az_interp[index];
        gauge_rhs_vars.alpha_sqrtg_Ax_interp[1] = alpha_sqrtg_Ax_interp[indexf(gridmax,i+1,j,k)];
        gauge_rhs_vars.alpha_sqrtg_Ay_interp[1] = alpha_sqrtg_Ay_interp[indexf(gridmax,i,j+1,k)];
        gauge_rhs_vars.alpha_sqrtg_Az_interp[1] = alpha_sqrtg_Az_interp[indexf(gridmax,i,j,k+1)];
    
        for(int iter=-2; iter<3; iter++) {
          const int indexx = indexf(gridmax,i+iter,j,     k     );
          const int indexy = indexf(gridmax,i,     j+iter,k     );
          const int indexz = indexf(gridmax,i,     j,     k+iter);
          gauge_rhs_vars.shiftx_interp[iter+2] = shiftx_interp[indexx];
          gauge_rhs_vars.shifty_interp[iter+2] = shifty_interp[indexy];
          gauge_rhs_vars.shiftz_interp[iter+2] = shiftz_interp[indexz];
          gauge_rhs_vars.phitildex[iter+2] = phitilde[indexx];
          gauge_rhs_vars.phitildey[iter+2] = phitilde[indexy];
          gauge_rhs_vars.phitildez[iter+2] = phitilde[indexz];
        }
        calculate_phitilde_and_A_gauge_rhs(Lorenz_damping_factor, &gauge_rhs_vars);
    
        phitilde_rhs[index] = gauge_rhs_vars.phitilde_rhs;
        Ax_rhs[index] += gauge_rhs_vars.A_x_gauge_rhs;
        Ay_rhs[index] += gauge_rhs_vars.A_y_gauge_rhs;
        Az_rhs[index] += gauge_rhs_vars.A_z_gauge_rhs;
  }


  infile = fopen("phitilde_and_A_gauge_rhs.bin", "rb");
  check_file_was_successfully_open(infile, "phitilde_and_A_gauge_rhs.bin");

  double *phitilde_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_trusted       = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_trusted       = (double*) malloc(sizeof(double)*arraylength);
  double *Az_trusted       = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(phitilde_trusted, sizeof(double), arraylength, infile);
  key += fread(Ax_trusted,       sizeof(double), arraylength, infile);
  key += fread(Ay_trusted,       sizeof(double), arraylength, infile);
  key += fread(Az_trusted,       sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*4)
    grhayl_error("An error has occured with reading in trusted data. Please check that comparison data\n"
                 "is up-to-date with current test version.\n");

#pragma omp parallel for
  for(int k=gridmin+3; k<gridmax-3; k++)
    for(int j=gridmin+3; j<gridmax-3; j++)
      for(int i=gridmin+3; i<gridmax-3; i++) {
        const int index = indexf(gridmax,i,j,k);
        if(rel_tol(reltol, phitilde_trusted[index], phitilde_rhs[index]))
          grhayl_error("Test unit_test_gauge_rhs has failed for variable phitilde_rhs at index (%d,%d,%d).\n", i, j, k);

        if(rel_tol(reltol, Ax_trusted[index], Ax_rhs[index]))
          grhayl_error("Test unit_test_gauge_rhs has failed for variable Ax_rhs at index (%d,%d,%d).\n", i, j, k);

        if(rel_tol(reltol, Ay_trusted[index], Ay_rhs[index]))
          grhayl_error("Test unit_test_gauge_rhs has failed for variable Ay_rhs at index (%d,%d,%d).\n", i, j, k);

        if(rel_tol(reltol, Az_trusted[index], Az_rhs[index]))
          grhayl_error("Test unit_test_gauge_rhs has failed for variable Az_rhs at index (%d,%d,%d).\n", i, j, k);
  }
  grhayl_info("Induction equation gauge RHS test has passed!\n");
  return 0;
}
