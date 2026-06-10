#include "../../utils_Noble.h"
#include "ghl_nrpyeos_tabulated.h"

#include <float.h>
#include <limits.h>

static int imin(int a, int b) { return a < b ? a : b; }
static int imax(int a, int b) { return a > b ? a : b; }

static int get_min_max(
      const int nx,
      const int ny,
      const int nz,
      const double *restrict x_arr,
      const double *restrict y_arr,
      const double *restrict f_arr,
      const double x,
      const double y,
      double *restrict f_min,
      double *restrict f_max) {

  if(x < x_arr[0] || x > x_arr[nx - 1] || y < y_arr[0] || y > y_arr[ny - 1]) {
    return 1;
  }

  const double dx = x_arr[1] - x_arr[0];
  const double dy = y_arr[1] - y_arr[0];

  const int ix = imin(imax((x - x_arr[0]) / dx, 0), nx - 2);
  const int iy = imin(imax((y - y_arr[0]) / dy, 0), ny - 2);

  const double tx = (x - x_arr[ix]) / dx;
  const double ty = (y - y_arr[iy]) / dy;

  const double w00 = (1.0 - tx) * (1.0 - ty);
  const double w10 = tx * (1.0 - ty);
  const double w01 = (1.0 - tx) * ty;
  const double w11 = tx * ty;

  double min = +DBL_MAX;
  double max = -DBL_MAX;

#define IDX(ix, iy, iz) ((ix) + nx*((iy) + ny*(iz)))

  for(int iz = 0; iz < nz; iz++) {
    const int i00 = IDX(ix + 0, iy + 0, iz);
    const int i01 = IDX(ix + 0, iy + 1, iz);
    const int i10 = IDX(ix + 1, iy + 0, iz);
    const int i11 = IDX(ix + 1, iy + 1, iz);

    const double f = w00 * f_arr[i00] +
                     w10 * f_arr[i10] +
                     w01 * f_arr[i01] +
                     w11 * f_arr[i11];

    if(f < min) {
      min = f;
    }
    if(f > max) {
      max = f;
    }
  }

  *f_min = min;
  *f_max = max;

  return 0;
}

static void enforce_bounds_rho_Ye_h(
  const ghl_eos_parameters *restrict eos,
  double *rho,
  double *Y_e,
  double *h
) {

  *rho = fmin(fmax(*rho, eos->rho_min), eos->rho_max);
  *Y_e = fmin(fmax(*Y_e, eos->Y_e_min), eos->Y_e_max);
  
  double logh_min = 0.0;
  double logh_max = 0.0;
  int err = get_min_max(eos->N_rho, eos->N_Ye, eos->N_T,
                        eos->table_logrho, eos->table_Y_e, eos->table_logh,
                        log(*rho), *Y_e, &logh_min, &logh_max);
  if(err) {
    ghl_error("Could not enforce table bounds\n");
  }

  *h = fmin(fmax(*h, exp(logh_min)), exp(logh_max));
}

static ghl_error_codes_t tabulated_initialize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const ghl_primitive_quantities *restrict prims,
      harm_aux_vars_struct *restrict harm_aux,
      double *restrict Z_ptr) {

  ghl_error_codes_t error = ghl_initialize_Noble(
        params, eos, ADM_metric, metric_aux, cons_undens, prims, harm_aux, Z_ptr);
  if(error) {
    return error;
  }

  const double v_valenciaU[3] = {
    ADM_metric->lapseinv * (prims->vU[0] + ADM_metric->betaU[0]),
    ADM_metric->lapseinv * (prims->vU[1] + ADM_metric->betaU[1]),
    ADM_metric->lapseinv * (prims->vU[2] + ADM_metric->betaU[2]),
  };
  double vsq = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, v_valenciaU);
  vsq = fmin(fmax(vsq, 0.0), 1.0 - 1e-15);
  double Wsq = 1.0 / (1.0 - vsq);
  double W = sqrt(Wsq);

  double rho = cons_undens->rho / W;
  double Y_e = cons_undens->Y_e / cons_undens->rho;
  double T = prims->temperature;
  ghl_tabulated_enforce_bounds_rho_Ye_T(eos, &rho, &Y_e, &T);

  double press = 0.0;
  double eps = 0.0;
  error = ghl_tabulated_compute_P_eps_from_T(eos, rho, Y_e, T, &press, &eps);
  if(error != ghl_success || !isfinite(press) || !isfinite(eps)) {
    return ghl_error_table_bisection;
  }

  harm_aux->Y_e = Y_e;
  harm_aux->temp_guess = T;

  const double w = rho * (1.0 + eps) + press;
  double Z_last = w * Wsq;

  *Z_ptr = Z_last;

  return ghl_success;
}

static void tabulated_func_2D(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const double dummy,
      const double x[],
      double dx[],
      double *restrict f,
      double *restrict df) {
  (void)dummy;

  const double Z = x[0];
  const double vsq = x[1];
  const double Wsq = 1.0 / (1.0 - vsq);
  const double W = sqrt(Wsq);
  const double rho = harm_aux->D / W;
  const double Y_e = harm_aux->Y_e;
  const double h = fabs(Z / (rho * Wsq));

  double press = 0.0;
  double eps = 0.0;
  double dPdrho = 0.0;
  double dPdeps = 0.0;
  double T = harm_aux->temp_guess;
  const ghl_error_codes_t error = ghl_tabulated_compute_P_eps_dPdrho_dPdeps_T_from_h(
                                     eos, rho, Y_e, h, &press, &eps, &dPdrho, &dPdeps, &T);

  double dPdvsq = NAN;
  double dPdZ = NAN;
  if(error == ghl_success) {
    harm_aux->temp_guess = T;
    const double dPdeps_over_rho = dPdeps / rho;

    dPdZ = (dPdeps_over_rho / (1.0 + dPdeps_over_rho)) / Wsq;

    const double dPdvsq_1 = -0.5 * harm_aux->D * W * dPdrho;
    const double dPdvsq_2 = -0.5 * (Z + press * Wsq) / rho;
    dPdvsq = (dPdvsq_1 + dPdeps * dPdvsq_2) / (1.0 + dPdeps_over_rho);
  }

  const double t2  = -0.5 * harm_aux->Bsq + dPdvsq;
  const double t3  = harm_aux->Bsq + Z;
  const double t4  = t3 * t3;
  const double t9  = 1.0 / (Z * Z);
  const double t11 = harm_aux->Qtsq - vsq * t4
                   + harm_aux->QdotBsq * (harm_aux->Bsq + 2.0 * Z) * t9;
  const double t16 = harm_aux->QdotBsq * t9;
  const double t18 = -harm_aux->Qdotn - 0.5 * harm_aux->Bsq * (1.0 + vsq)
                   + 0.5 * t16 - Z + press;
  const double t21 = 1.0 / t3;
  const double t23 = 1.0 / Z;
  const double t24 = t16 * t23;
  const double t25 = -1.0 + dPdZ - t24;
  const double t35 = t25 * t3
                   + (harm_aux->Bsq - 2.0 * dPdvsq)
                   * (harm_aux->QdotBsq + vsq * Z * Z * Z) * t9 * t23;
  const double t36 = 1.0 / t35;
  const double t40 = (vsq + t24) * t3;

  dx[0] = -(t2 * t11 + t4 * t18) * t21 * t36;
  dx[1] = -(-t25 * t11 - 2.0 * t40 * t18) * t21 * t36;

  *df = -t11 * t11 - t18 * t18;
  *f = -0.5 * (*df);

}

static ghl_error_codes_t tabulated_finalize_Noble(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      const harm_aux_vars_struct *restrict harm_aux,
      const double Z,
      const double vsq,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {
  (void)metric_aux;

  const double gtmp = sqrt(1.0 - vsq);
  double W = 1.0 / gtmp;
  const double w = Z * (1.0 - vsq);

  const double nU[4] = {
    ADM_metric->lapseinv,
   -ADM_metric->lapseinv * ADM_metric->betaU[0],
   -ADM_metric->lapseinv * ADM_metric->betaU[1],
   -ADM_metric->lapseinv * ADM_metric->betaU[2],
  };

  const double g_o_ZBsq = W / (Z + harm_aux->Bsq);
  const double QdB_o_Z = harm_aux->QdotB / Z;
  const double Qtcon[4] = {
    harm_aux->QU[0] + nU[0] * harm_aux->Qdotn,
    harm_aux->QU[1] + nU[1] * harm_aux->Qdotn,
    harm_aux->QU[2] + nU[2] * harm_aux->Qdotn,
    harm_aux->QU[3] + nU[3] * harm_aux->Qdotn,
  };
  double utU[3] = {
    g_o_ZBsq * (Qtcon[1] + QdB_o_Z * prims->BU[0]),
    g_o_ZBsq * (Qtcon[2] + QdB_o_Z * prims->BU[1]),
    g_o_ZBsq * (Qtcon[3] + QdB_o_Z * prims->BU[2]),
  };

  diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(params, ADM_metric, utU, prims);

  if(diagnostics->speed_limited) {
    // Recompute W so its compatible with new velocity
    W = ADM_metric->lapse * prims->u0;
  }

  prims->rho = harm_aux->D / W;
  prims->Y_e = cons_undens->Y_e / cons_undens->rho;

  double h = w / prims->rho;
  double press = 0.0;
  double eps = 0.0;
  double S = 0.0;
  double T = harm_aux->temp_guess;
  enforce_bounds_rho_Ye_h(eos, &prims->rho, &prims->Y_e, &h);
  const ghl_error_codes_t error = ghl_tabulated_compute_P_eps_S_T_from_h(
                                     eos, prims->rho, prims->Y_e, h,
                                     &press, &eps, &S, &T);
  if(error != ghl_success) {
    return ghl_error_table_bisection;
  }

  prims->press = press;
  prims->eps = eps;
  prims->entropy = S;
  prims->temperature = T;

  return ghl_success;
}

ghl_error_codes_t ghl_tabulated_Noble2D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double gnr_out[2];
  harm_aux_vars_struct harm_aux;

  ghl_error_codes_t error = tabulated_initialize_Noble(
        params, eos, ADM_metric, metric_aux, cons_undens, prims, &harm_aux, &gnr_out[0]);
  if(error) {
    return error;
  }

  const double v_valenciaU[3] = {
    ADM_metric->lapseinv * (prims->vU[0] + ADM_metric->betaU[0]),
    ADM_metric->lapseinv * (prims->vU[1] + ADM_metric->betaU[1]),
    ADM_metric->lapseinv * (prims->vU[2] + ADM_metric->betaU[2]),
  };
  double vsq = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, v_valenciaU);
  gnr_out[1] = fmin(fmax(vsq, 0.0), 1.0);

  const ghl_error_codes_t retval = ghl_general_newton_raphson(
        eos, &harm_aux, 2, 0.0, gnr_out, ghl_validate_2D, tabulated_func_2D);
  const double Z = gnr_out[0];

  if(retval != ghl_success) {
    return retval;
  }
  if(Z <= 0.0 || Z > 1e20) {
    return ghl_error_invalid_Z;
  }

  vsq = gnr_out[1];
  if(vsq >= 1.0) {
    vsq = 1.0 - 2.e-16;
  }
  else if(vsq < 0.0) {
    return ghl_error_neg_vsq;
  }

  error = tabulated_finalize_Noble(
        params, eos, ADM_metric, metric_aux, cons_undens, &harm_aux, Z, vsq, prims, diagnostics);
  if(error) {
    return error;
  }
  if(prims->press <= 0.0) {
    return ghl_error_neg_pressure;
  }

  diagnostics->n_iter = harm_aux.n_iter;
  diagnostics->which_routine = Noble2D;

  return ghl_success;
}
