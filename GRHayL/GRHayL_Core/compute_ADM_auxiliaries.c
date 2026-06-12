#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Compute ADM auxiliary metric quantities from an ADM metric struct
 *
 * @details
 * This function takes a ghl_metric_quantities struct with ADM metric
 * data and computes all the elements of the ghl_ADM_aux_quantities
 * struct.
 *
 * @param[in] metric_adm pointer to a ghl_metric_quantities struct.
 *                         This contains the input ADM metric data.
 *
 * @param[out] metric_aux pointer to a ghl_ADM_aux_quantities struct.
 *                         This will be filled with computed auxiliary quantities.
 *
 */
void ghl_compute_ADM_auxiliaries(
      const ghl_metric_quantities *restrict metric_adm,
      ghl_ADM_aux_quantities *restrict metric_aux) {

  double betaD[3];
  ghl_raise_lower_vector_3D(metric_adm->gammaDD, metric_adm->betaU, betaD);

  const double shift2 = betaD[0]*metric_adm->betaU[0]
                      + betaD[1]*metric_adm->betaU[1]
                      + betaD[2]*metric_adm->betaU[2];

  const double shift_over_lapse2[3] = {metric_adm->betaU[0]*metric_adm->lapseinv2,
                                       metric_adm->betaU[1]*metric_adm->lapseinv2,
                                       metric_adm->betaU[2]*metric_adm->lapseinv2};

  metric_aux->g4DD[0][0]                          = -SQR(metric_adm->lapse) + shift2;
  metric_aux->g4DD[0][1] = metric_aux->g4DD[1][0] = betaD[0];
  metric_aux->g4DD[0][2] = metric_aux->g4DD[2][0] = betaD[1];
  metric_aux->g4DD[0][3] = metric_aux->g4DD[3][0] = betaD[2];
  metric_aux->g4DD[1][1]                          = metric_adm->gammaDD[0][0];
  metric_aux->g4DD[1][2] = metric_aux->g4DD[2][1] = metric_adm->gammaDD[0][1];
  metric_aux->g4DD[1][3] = metric_aux->g4DD[3][1] = metric_adm->gammaDD[0][2];
  metric_aux->g4DD[2][2]                          = metric_adm->gammaDD[1][1];
  metric_aux->g4DD[2][3] = metric_aux->g4DD[3][2] = metric_adm->gammaDD[1][2];
  metric_aux->g4DD[3][3]                          = metric_adm->gammaDD[2][2];

  metric_aux->g4UU[0][0]                          = -metric_adm->lapseinv2;
  metric_aux->g4UU[0][1] = metric_aux->g4UU[1][0] = shift_over_lapse2[0];
  metric_aux->g4UU[0][2] = metric_aux->g4UU[2][0] = shift_over_lapse2[1];
  metric_aux->g4UU[0][3] = metric_aux->g4UU[3][0] = shift_over_lapse2[2];
  metric_aux->g4UU[1][1]                          = metric_adm->gammaUU[0][0] - metric_adm->betaU[0]*shift_over_lapse2[0];
  metric_aux->g4UU[1][2] = metric_aux->g4UU[2][1] = metric_adm->gammaUU[0][1] - metric_adm->betaU[0]*shift_over_lapse2[1];
  metric_aux->g4UU[1][3] = metric_aux->g4UU[3][1] = metric_adm->gammaUU[0][2] - metric_adm->betaU[0]*shift_over_lapse2[2];
  metric_aux->g4UU[2][2]                          = metric_adm->gammaUU[1][1] - metric_adm->betaU[1]*shift_over_lapse2[1];
  metric_aux->g4UU[2][3] = metric_aux->g4UU[3][2] = metric_adm->gammaUU[1][2] - metric_adm->betaU[1]*shift_over_lapse2[2];
  metric_aux->g4UU[3][3]                          = metric_adm->gammaUU[2][2] - metric_adm->betaU[2]*shift_over_lapse2[2];
}
