#include "../../utils_Noble.h"

GHL_DEVICE
ghl_error_codes_t ghl_initialize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      harm_aux_vars_struct *restrict harm_aux,
      double *restrict Z_ptr) {

  harm_aux->n_iter = 0;
  harm_aux->max_iterations = params->con2prim_max_iterations;
  harm_aux->solver_tolerance = params->con2prim_solver_tolerance;

  harm_aux->Bsq = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims->BU);

  const double uu = - cons_undens->tau*ADM_metric->lapse
                    - ADM_metric->lapse*cons_undens->rho
                    + ADM_metric->betaU[0]*cons_undens->SD[0]
                    + ADM_metric->betaU[1]*cons_undens->SD[1]
                    + ADM_metric->betaU[2]*cons_undens->SD[2];

  const double QD[4] = {uu,
                        cons_undens->SD[0],
                        cons_undens->SD[1],
                        cons_undens->SD[2]};

  ghl_raise_lower_vector_4D(metric_aux->g4UU, QD, harm_aux->QU);
  const double Qsq = QD[0]*harm_aux->QU[0]
                   + QD[1]*harm_aux->QU[1]
                   + QD[2]*harm_aux->QU[2]
                   + QD[3]*harm_aux->QU[3];

  harm_aux->QdotB = QD[1]*prims->BU[0]
                  + QD[2]*prims->BU[1]
                  + QD[3]*prims->BU[2];

  harm_aux->QdotBsq = harm_aux->QdotB*harm_aux->QdotB;

  // n_{\mu}Q^{\mu} = -alpha Q^{0}, since n_{\mu} = (-alpha,0,0,0)
  harm_aux->Qdotn = -ADM_metric->lapse*harm_aux->QU[0];

  harm_aux->Qtsq = Qsq + harm_aux->Qdotn*harm_aux->Qdotn;

  harm_aux->D    = cons_undens->rho;

  /* calculate Z from last timestep and use for guess */
  const double utU_guess[3] = {prims->u0*(prims->vU[0] + ADM_metric->betaU[0]),
                               prims->u0*(prims->vU[1] + ADM_metric->betaU[1]),
                               prims->u0*(prims->vU[2] + ADM_metric->betaU[2])};
  const double tmp_u = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, utU_guess);

  double utsq = tmp_u/(1.0-tmp_u);

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) {
    utsq = fabs(utsq);
  } else if(utsq < 0.0 || utsq > 10.0) {
    return ghl_error_invalid_utsq;
  }

  const double Wsq = 1.0 + utsq;   // Lorentz factor squared
  const double W = sqrt(Wsq);

  // Always calculate rho from D and W so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(utsq).
  const double rho0 = harm_aux->D / W;

  /*  Assumes hybrid EOS */
  double w = rho0*(1.0 + prims->eps) + prims->press;
  /*  Assumes hybrid EOS */

  double Z_last = w*Wsq;

  // Make sure that Z is large enough so that v^2 < 1 :
  int i_increase = 0;
  while( (( Z_last*Z_last*Z_last * ( Z_last + 2.*harm_aux->Bsq )
            - harm_aux->QdotBsq*(2.*Z_last + harm_aux->Bsq) ) <= Z_last*Z_last*(harm_aux->Qtsq-harm_aux->Bsq*harm_aux->Bsq))
         && (i_increase < 10) ) {
    Z_last *= 10.;
    i_increase++;
  }
  *Z_ptr = Z_last;
  return ghl_success;
}

GHL_DEVICE
ghl_error_codes_t ghl_initialize_Noble_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      harm_aux_vars_struct *restrict harm_aux,
      double *restrict rho_ptr,
      double *restrict Z_ptr) {

  harm_aux->n_iter = 0;
  harm_aux->max_iterations = params->con2prim_max_iterations;
  harm_aux->solver_tolerance = params->con2prim_solver_tolerance;

  harm_aux->Bsq = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims->BU);

  harm_aux->W_times_S = cons_undens->entropy;

  const double uu = - cons_undens->tau*ADM_metric->lapse
                    - ADM_metric->lapse*cons_undens->rho
                    + ADM_metric->betaU[0]*cons_undens->SD[0]
                    + ADM_metric->betaU[1]*cons_undens->SD[1]
                    + ADM_metric->betaU[2]*cons_undens->SD[2];

  const double QD[4] = {uu,
                        cons_undens->SD[0],
                        cons_undens->SD[1],
                        cons_undens->SD[2]};

  ghl_raise_lower_vector_4D(metric_aux->g4UU, QD, harm_aux->QU);
  const double Qsq = QD[0]*harm_aux->QU[0]
                   + QD[1]*harm_aux->QU[1]
                   + QD[2]*harm_aux->QU[2]
                   + QD[3]*harm_aux->QU[3];

  harm_aux->QdotB = QD[1]*prims->BU[0]
                  + QD[2]*prims->BU[1]
                  + QD[3]*prims->BU[2];

  harm_aux->QdotBsq = harm_aux->QdotB*harm_aux->QdotB;

  // n_{\mu}Q^{\mu} = -alpha Q^{0}, since n_{\mu} = (-alpha,0,0,0)
  harm_aux->Qdotn = -ADM_metric->lapse*harm_aux->QU[0];

  harm_aux->Qtsq = Qsq + harm_aux->Qdotn*harm_aux->Qdotn;

  harm_aux->D    = cons_undens->rho;

  /* calculate Z from last timestep and use for guess */
  const double utU_guess[3] = {prims->u0*(prims->vU[0] + ADM_metric->betaU[0]),
                               prims->u0*(prims->vU[1] + ADM_metric->betaU[1]),
                               prims->u0*(prims->vU[2] + ADM_metric->betaU[2])};
  const double tmp_u = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, utU_guess);

  double utsq = tmp_u/(1.0-tmp_u);

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) {
    utsq = fabs(utsq);
  } else if(utsq < 0.0 || utsq > 10.0) {
    return ghl_error_invalid_utsq;
  }

  const double Wsq = 1.0 + utsq;   // Lorentz factor squared
  const double W = sqrt(Wsq);

  // Always calculate rho from D and W so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq).
  const double rho0 = harm_aux->D / W;
  const double Gamma_ppoly = eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(eos, rho0)];
  const double Gm1        = Gamma_ppoly - 1.0;

  const double rho_Gm1     = pow(rho0,Gm1);     // HARM auxiliary variable

  /* The definition of the entropy density, S, is
   *
   * S = P / rho^(Gamma - 1)
   *
   * Thus we have
   * .-------------------------.
   * | P = rho^(Gamma - 1) * S |
   * .-------------------------.
   */
  const double p = cons_undens->entropy * rho_Gm1/W;
  const double w = rho0 + p/Gm1 + p;

  *rho_ptr = rho0;
  *Z_ptr = w*Wsq;
  return ghl_success;
}
