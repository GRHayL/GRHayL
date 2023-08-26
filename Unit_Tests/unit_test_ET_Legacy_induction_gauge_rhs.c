#include "unit_tests.h"

int main(int argc, char **argv) {
  FILE* infile = fopen_with_check("ET_Legacy_induction_gauge_rhs_input.bin","rb");

  int dirlength;
  int key = fread(&dirlength, sizeof(int), 1, infile);
  if( key != 1 || dirlength < 1 )
    ghl_error("An error has occured with reading the grid size. "
                 "Please check that Noble2D_initial_data.bin"
                 "is up-to-date with current test version.\n");
  const int arraylength = dirlength*dirlength*dirlength;

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
  double *betax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *betaz_interp = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_Phi_minus_betaj_A_j_interp = (double*) malloc(sizeof(double)*arraylength);
  double *sqrtg_Ax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *sqrtg_Ay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *sqrtg_Az_interp = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Az_rhs = (double*) malloc(sizeof(double)*arraylength);

  // Read in data from file to ensure portability
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
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Data which should be written before it is used is poisoned to ghl_pert_test_fail behavior.
  // RHSs for A are set to 0 because they are assumed to already contain the zero-gauge
  // contribution to the RHS before entering this function.
  const double poison = 1e300;

  for(int k=0; k<dirlength; k++) {
    for(int j=0; j<dirlength; j++) {
      for(int i=0; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        alpha_interp[index]  = poison;
        betax_interp[index] = poison;
        betay_interp[index] = poison;
        betaz_interp[index] = poison;

        alpha_Phi_minus_betaj_A_j_interp[index] = poison;
        sqrtg_Ax_interp[index]                  = poison;
        sqrtg_Ay_interp[index]                  = poison;
        sqrtg_Az_interp[index]                  = poison;

        Ax_rhs[index]       = 0.0;
        Ay_rhs[index]       = 0.0;
        Az_rhs[index]       = 0.0;
        phitilde_rhs[index] = poison;
      }
    }
  }

  for(int k=1; k<dirlength-1; k++) {
    for(int j=1; j<dirlength-1; j++) {
      for(int i=1; i<dirlength-1; i++) {
        const int index = indexf(dirlength,i,j,k);

        ghl_metric_quantities metric_stencil[2][2][2];
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
              metric_stencil[iterz][itery][iterx].gammaUU[0][0]  = gupxx[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[0][1]  = gupxy[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[0][2]  = gupxz[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[1][1]  = gupyy[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[1][2]  = gupyz[ind];
              metric_stencil[iterz][itery][iterx].gammaUU[2][2]  = gupzz[ind];
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
// This code should only copy the needed data that isn't copied in the loop for other variables, but it is untested.
//        for(int iter2=0; iter2<2; iter2++)
//        for(int iter1=0; iter1<2; iter1++) {
//          gauge_vars.A_x[iter2+1][0][iter1+1] = in_vars[A_XI][indexf(dirlength, i+iter1,     j-1, k+iter2)]; // { (0,1),    -1, (0,1)}
//          gauge_vars.A_x[0][iter2+1][iter1+1] = in_vars[A_XI][indexf(dirlength, i+iter1, j+iter2, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_y[iter2+1][iter1+1][0] = in_vars[A_YI][indexf(dirlength,     i-1, j+iter1, k+iter2)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_y[0][iter1+1][iter2+1] = in_vars[A_YI][indexf(dirlength, i+iter2, j+iter1, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_z[iter1+1][iter2+1][0] = in_vars[A_ZI][indexf(dirlength,     i-1, j+iter2, k+iter1)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_z[iter1+1][0][iter2+1] = in_vars[A_ZI][indexf(dirlength, i+iter2,     j-1, k+iter1)]; // { (0,1),    -1, (0,1)}
//        }

        ghl_interpolate_with_cell_centered_BSSN(metric_stencil, psi_stencil, Ax_stencil, Ay_stencil, Az_stencil, phitilde[index], &interp_vars);

        alpha_interp[index] = interp_vars.alpha;
        sqrtg_Ax_interp[index] = interp_vars.sqrtg_Ai[0];
        sqrtg_Ay_interp[index] = interp_vars.sqrtg_Ai[1];
        sqrtg_Az_interp[index] = interp_vars.sqrtg_Ai[2];
        alpha_Phi_minus_betaj_A_j_interp[index] = interp_vars.alpha_Phi_minus_betaj_A_j;
        betax_interp[index] = interp_vars.betai[0];
        betay_interp[index] = interp_vars.betai[1];
        betaz_interp[index] = interp_vars.betai[2];
      }
    }
  }

  const double dxinv[3] = {1.0/dX[0], 1.0/dX[1], 1.0/dX[2]};

  // This loop requires two additional ghostzones in every direction. Hence the following loop definition:
  for(int k=3; k<dirlength-3; k++) {
    for(int j=3; j<dirlength-3; j++) {
      for(int i=3; i<dirlength-3; i++) {
        const int index = indexf(dirlength,i,j,k);

        // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
        //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
        Ax_rhs[index] += dxinv[0]*(alpha_Phi_minus_betaj_A_j_interp[indexf(dirlength,i-1,j,k)] - alpha_Phi_minus_betaj_A_j_interp[index]);
        Ay_rhs[index] += dxinv[1]*(alpha_Phi_minus_betaj_A_j_interp[indexf(dirlength,i,j-1,k)] - alpha_Phi_minus_betaj_A_j_interp[index]);
        Az_rhs[index] += dxinv[2]*(alpha_Phi_minus_betaj_A_j_interp[indexf(dirlength,i,j,k-1)] - alpha_Phi_minus_betaj_A_j_interp[index]);

        double betax[5], betay[5], betaz[5];
        double phitilde_stencil[3][5], sqrtg_Ai_stencil[3][2];

        sqrtg_Ai_stencil[0][0] = sqrtg_Ax_interp[index];
        sqrtg_Ai_stencil[1][0] = sqrtg_Ay_interp[index];
        sqrtg_Ai_stencil[2][0] = sqrtg_Az_interp[index];

        sqrtg_Ai_stencil[0][1] = sqrtg_Ax_interp[indexf(dirlength,i+1,j,k)];
        sqrtg_Ai_stencil[1][1] = sqrtg_Ay_interp[indexf(dirlength,i,j+1,k)];
        sqrtg_Ai_stencil[2][1] = sqrtg_Az_interp[indexf(dirlength,i,j,k+1)];

        for(int iter=-2; iter<3; iter++) {
          const int indexx = indexf(dirlength,i+iter,j,     k     );
          const int indexy = indexf(dirlength,i,     j+iter,k     );
          const int indexz = indexf(dirlength,i,     j,     k+iter);
          betax[iter+2] = betax_interp[indexx];
          betay[iter+2] = betay_interp[indexy];
          betaz[iter+2] = betaz_interp[indexz];
          phitilde_stencil[0][iter+2] = phitilde[indexx];
          phitilde_stencil[1][iter+2] = phitilde[indexy];
          phitilde_stencil[2][iter+2] = phitilde[indexz];
        }
        phitilde_rhs[index] = ghl_calculate_phitilde_rhs(dxinv, Lorenz_damping_factor, alpha_interp[index], betax, betay, betaz, sqrtg_Ai_stencil, phitilde_stencil);
      }
    }
  }


  infile = fopen_with_check("ET_Legacy_induction_gauge_rhs_output.bin", "rb");

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
    ghl_error("An error has occured with reading in trusted data. Please check that comparison data\n"
                 "is up-to-date with current test version.\n");

  infile = fopen_with_check("ET_Legacy_induction_gauge_rhs_output_pert.bin", "rb");

  double *phitilde_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_pert       = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_pert       = (double*) malloc(sizeof(double)*arraylength);
  double *Az_pert       = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(phitilde_pert, sizeof(double), arraylength, infile);
  key += fread(Ax_pert,       sizeof(double), arraylength, infile);
  key += fread(Ay_pert,       sizeof(double), arraylength, infile);
  key += fread(Az_pert,       sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*4)
    ghl_error("An error has occured with reading in perturbed data. Please check that comparison data\n"
                 "is up-to-date with current test version.\n");

  for(int k=3; k<dirlength-3; k++)
    for(int j=3; j<dirlength-3; j++)
      for(int i=3; i<dirlength-3; i++) {
        const int index = indexf(dirlength,i,j,k);

        if( ghl_pert_test_fail(phitilde_trusted[index], phitilde_rhs[index], phitilde_pert[index]) )
          ghl_error("Test unit_test_ET_Legacy_induction_gauge_rhs has failed for variable phitilde_rhs.\n"
                       "  phitilde trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", phitilde_trusted[index], phitilde_rhs[index], phitilde_pert[index],
                                                   relative_error(phitilde_trusted[index], phitilde_rhs[index]),
                                                   relative_error(phitilde_trusted[index], phitilde_pert[index]));

        if( ghl_pert_test_fail(Ax_trusted[index], Ax_rhs[index], Ax_pert[index]) )
          ghl_error("Test unit_test_ET_Legacy_induction_gauge_rhs has failed for variable Ax_rhs.\n"
                       "  Ax trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", Ax_trusted[index], Ax_rhs[index], Ax_pert[index],
                                                   relative_error(Ax_trusted[index], Ax_rhs[index]),
                                                   relative_error(Ax_trusted[index], Ax_pert[index]));

        if( ghl_pert_test_fail(Ay_trusted[index], Ay_rhs[index], Ay_pert[index]) )
          ghl_error("Test unit_test_ET_Legacy_induction_gauge_rhs has failed for variable Ay_rhs.\n"
                       "  Ay trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", Ay_trusted[index], Ay_rhs[index], Ay_pert[index],
                                                   relative_error(Ay_trusted[index], Ay_rhs[index]),
                                                   relative_error(Ay_trusted[index], Ay_pert[index]));

        if( ghl_pert_test_fail(Az_trusted[index], Az_rhs[index], Az_pert[index]) )
          ghl_error("Test unit_test_ET_Legacy_induction_gauge_rhs has failed for variable Az_rhs.\n"
                       "  Az trusted %.14e computed %.14e perturbed %.14e\n"
                       "  rel.err. %.14e %.14e\n", Az_trusted[index], Az_rhs[index], Az_pert[index],
                                                   relative_error(Az_trusted[index], Az_rhs[index]),
                                                   relative_error(Az_trusted[index], Az_pert[index]));
  }
  ghl_info("ET_Legacy induction gauge RHS test has passed!\n");
  free(gupxx); free(gupxy); free(gupxz);
  free(gupyy); free(gupyz); free(gupzz);
  free(psi); free(lapse);
  free(betax); free(betay); free(betaz);
  free(phitilde);
  free(Ax); free(Ay); free(Az);
  free(alpha_interp);
  free(betax_interp); free(betay_interp); free(betaz_interp);
  free(alpha_Phi_minus_betaj_A_j_interp);
  free(sqrtg_Ax_interp); free(sqrtg_Ay_interp); free(sqrtg_Az_interp);
  free(phitilde_rhs);
  free(Ax_rhs); free(Ay_rhs); free(Az_rhs);
  free(phitilde_trusted);
  free(Ax_trusted); free(Ay_trusted); free(Az_trusted);
  free(phitilde_pert);
  free(Ax_pert); free(Ay_pert); free(Az_pert);
}
