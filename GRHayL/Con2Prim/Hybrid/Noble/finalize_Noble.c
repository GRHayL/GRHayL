#include "../../utils_Noble.h"

int ghl_finalize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const harm_aux_vars_struct *restrict harm_aux,
      const double Z,
      const double vsq,
      ghl_primitive_quantities *restrict prims) {

  const double gtmp = sqrt(1. - vsq);
  const double W = 1.0/gtmp;
  const double w = Z * (1.0 - vsq);

  const double nU[4] = {ADM_metric->lapseinv,
                       -ADM_metric->lapseinv*ADM_metric->betaU[0],
                       -ADM_metric->lapseinv*ADM_metric->betaU[1],
                       -ADM_metric->lapseinv*ADM_metric->betaU[2]};

  const double g_o_ZBsq = W/(Z+harm_aux->Bsq);
  const double QdB_o_Z  = harm_aux->QdotB / Z;

  const double Qtcon[4] = {harm_aux->QU[0] + nU[0] * harm_aux->Qdotn,
                           harm_aux->QU[1] + nU[1] * harm_aux->Qdotn,
                           harm_aux->QU[2] + nU[2] * harm_aux->Qdotn,
                           harm_aux->QU[3] + nU[3] * harm_aux->Qdotn};

  double utU[3] = {g_o_ZBsq * ( Qtcon[1] + QdB_o_Z*harm_aux->BbarU[0] ),
                   g_o_ZBsq * ( Qtcon[2] + QdB_o_Z*harm_aux->BbarU[1] ),
                   g_o_ZBsq * ( Qtcon[3] + QdB_o_Z*harm_aux->BbarU[2] )};

  const int speed_limited = ghl_limit_utilde_and_compute_v(params, ADM_metric, utU, prims);

  prims->rho = harm_aux->D * gtmp;
  if(speed_limited==1)
    prims->rho = cons_undens->rho/(ADM_metric->lapse*prims->u0);

  /*  Assumes hybrid EOS */
  prims->press = ghl_pressure_rho0_w(eos, prims->rho, w);
  prims->eps = ghl_hybrid_compute_epsilon(eos, prims->rho, prims->press);
  if( params->evolve_entropy )
    prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);
  /*  Assumes hybrid EOS */

  return speed_limited;
}
