#include "ghl.h"

/* Function    : ghl_compute_ADM_auxiliaries()
 * Description : Initialize the metric struct from user input
 *
 * Inputs      : lapse          - value of the lapse
 *             : gij            - value of the (i,j) component of the
 *                                cartesian ADM metric g_ij
 *             : betax          - value of the x component of the shift
 *             : betay          - value of the y component of the shift
 *             : betaz          - value of the z component of the shift
 *
 * Outputs     : metric         - returns metric_quantities struct containing
 *                                the inputs and additional auxiliary data computed
 *                                from input
 */

void ghl_compute_ADM_auxiliaries(
      const metric_quantities *restrict ADM_metric,
      ADM_aux_quantities *restrict metric_aux) {

  metric_aux->phi = (1.0/12.0) * log(ADM_metric->gijdet);

  metric_aux->psi2 = exp(2.0*metric_aux->phi);
  metric_aux->psi4 = SQR(metric_aux->psi2);
  metric_aux->psi6 = metric_aux->psi4*metric_aux->psi2;
  metric_aux->psi4inv = 1.0/metric_aux->psi4;

  double betaD[3];
  ghl_lower_vector_3D(ADM_metric->gammaDD, ADM_metric->betaU, betaD);

  const double shift2 = betaD[0]*ADM_metric->betaU[0]
                      + betaD[1]*ADM_metric->betaU[1]
                      + betaD[2]*ADM_metric->betaU[2];

  const double shift_over_lapse2[3] = {ADM_metric->betaU[0]*ADM_metric->lapseinv2,
                                       ADM_metric->betaU[1]*ADM_metric->lapseinv2,
                                       ADM_metric->betaU[2]*ADM_metric->lapseinv2};

  metric_aux->g4DD[0][0]                          = -SQR(ADM_metric->lapse) + shift2;
  metric_aux->g4DD[0][1] = metric_aux->g4DD[1][0] = betaD[0];
  metric_aux->g4DD[0][2] = metric_aux->g4DD[2][0] = betaD[1];
  metric_aux->g4DD[0][3] = metric_aux->g4DD[3][0] = betaD[2];
  metric_aux->g4DD[1][1]                          = ADM_metric->gammaDD[0][0];
  metric_aux->g4DD[1][2] = metric_aux->g4DD[2][1] = ADM_metric->gammaDD[0][1];
  metric_aux->g4DD[1][3] = metric_aux->g4DD[3][1] = ADM_metric->gammaDD[0][2];
  metric_aux->g4DD[2][2]                          = ADM_metric->gammaDD[1][1];
  metric_aux->g4DD[2][3] = metric_aux->g4DD[3][2] = ADM_metric->gammaDD[1][2];
  metric_aux->g4DD[3][3]                          = ADM_metric->gammaDD[2][2];

  metric_aux->g4UU[0][0]                          = -ADM_metric->lapseinv2;
  metric_aux->g4UU[0][1] = metric_aux->g4UU[1][0] = shift_over_lapse2[0];
  metric_aux->g4UU[0][2] = metric_aux->g4UU[2][0] = shift_over_lapse2[1];
  metric_aux->g4UU[0][3] = metric_aux->g4UU[3][0] = shift_over_lapse2[2];
  metric_aux->g4UU[1][1]                          = ADM_metric->gammaUU[0][0] - ADM_metric->betaU[0]*shift_over_lapse2[0];
  metric_aux->g4UU[1][2] = metric_aux->g4UU[2][1] = ADM_metric->gammaUU[0][1] - ADM_metric->betaU[0]*shift_over_lapse2[1];
  metric_aux->g4UU[1][3] = metric_aux->g4UU[3][1] = ADM_metric->gammaUU[0][2] - ADM_metric->betaU[0]*shift_over_lapse2[2];
  metric_aux->g4UU[2][2]                          = ADM_metric->gammaUU[1][1] - ADM_metric->betaU[1]*shift_over_lapse2[1];
  metric_aux->g4UU[2][3] = metric_aux->g4UU[3][2] = ADM_metric->gammaUU[1][2] - ADM_metric->betaU[1]*shift_over_lapse2[2];
  metric_aux->g4UU[3][3]                          = ADM_metric->gammaUU[2][2] - ADM_metric->betaU[2]*shift_over_lapse2[2];
}
