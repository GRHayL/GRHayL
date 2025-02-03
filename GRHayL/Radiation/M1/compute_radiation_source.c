#include "ghl.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

// gsl_multiroots: https://www.gnu.org/software/gsl/doc/html/multiroots.html
void init_params(ghl_m1_powell_params *p,
                 const ghl_m1_closure_t closure,
                 const gsl_root_fsolver *gsl_solver_1d,
                 const double cdt,
                 const ghl_metric_quantities *metric,
                 const ghl_ADM_aux_quantities *adm_aux,
                 const ghl_radiation_metric_tensor *proj4,
                 const ghl_primitive_quantities *prims,
                 const double * u4D,
                 const double * u4U,
                 const double * n4D,
                 const double * n4U,
                 const double J,
                 const ghl_radiation_flux_vector *F4,
                 const ghl_radiation_flux_vector *H4,
                 const ghl_radiation_pressure_tensor *P4,
                 const ghl_stress_energy *rT4DD,
                 const ghl_radiation_con_source_vector *S4,
                 const double W,
                 const double E_star,
                 const ghl_radiation_flux_vector *F4_star,
                 const double Edot,
                 const double chi,
                 const double eta,
                 const double kabs,
                 const double kscat,
                 const double E,
                 const double rE_source,
                 const ghl_radiation_con_flux_vector *rF_source,
                 const double GE_source,
                 const ghl_radiation_con_source_vector * F_src) {
    p->closure = closure;
    p->gsl_solver_1d = gsl_solver_1d;
    p->cdt = cdt;
    p->metric = metric;
    p->adm_aux = adm_aux;
    p->proj4 = proj4;
    p->prims = prims;
    p->u4D = u4D;
    p->u4U = u4U;
    p->n4D = n4D;
    p->n4U = n4U;
    p->J = J;
    p->F4 = F4;
    p->H4 = H4;
    p->P4 = P4;
    p->rT4DD = rT4DD;
    p->W = W;
    p->E_star = E_star;
    p->F4_star = F4_star;
    p->Edot = Edot;
    p->chi = chi;
    p->eta = eta;
    p->kabs = kabs;
    p->kscat = kscat;
    p->E = E;
    p->rE_source = rE_source;
    p->rF_source = rF_source;
    p->GE_source = GE_source;
    p->F_src = F_src;
}

// rad_E_floor rad_eps are parameters in param.ccl
// this function modifies E and F4 given the floor values
void apply_floor(
      const ghl_ADM_aux_quantities *adm_aux,
      double *E,
      ghl_radiation_flux_vector *F4,
      const double rad_E_floor,
      const double rad_eps) {
  *E = fmax(rad_E_floor, *E);

  double F2 = 0; // needs a better name
  for(int a = 0; a < 4; ++a) {
    for(int b = 0; b < 4; ++b) {
      F2 += adm_aux->g4UU[a][b] * F4->D[a] * F4->D[b];
    }
  }
  const double lim = (*E) * (*E) * (1 - rad_eps);
  if(F2 > lim) {
    double fac = lim / F2;
    for(int a = 0; a < 4; ++a) {
      F4->D[a] *= fac;
    }
  }
}
// // rad_eps "Impose F_a F^a < (1 - rad_E_eps) E2"
// // rad_E_floor "Radiation energy density floor"
// // TODO: CCTK_param in THC: what to do here?
// double rad_eps = 1.0e-15; // default value in THC
// double rad_E_floor = 1.0e-15;
// apply_floor(adm_aux, &E, &F4, rad_E_floor, rad_eps); // TODO: what does this step do?

/* return dthin */
inline double radM1_set_dthin(double chi) { return 1.5 * chi - 0.5; }

/* return dthick */
inline double radM1_set_dthick(double chi) { return 1.5 * (1 - chi); }

// Low level kernel computing the Jacobian matrix
void __source_jacobian_low_level(
      double *qpre,
      double Fup[4],
      double F2,
      double chi,
      double kapa,
      double kaps,
      double vup[4],
      double vdown[4],
      double v2,
      double W,
      double alpha,
      double cdt,
      double *qstar,
      ghl_radiation_con_source_tensor *J) {
  const double kapas = kapa + kaps;
  const double alpW = alpha * W;

  const double dthin = radM1_set_dthin(chi);
  const double dthick = radM1_set_dthick(chi);

  const double vx = vdown[1];
  const double vy = vdown[2];
  const double vz = vdown[3];
  const double W2 = W * W;
  const double W3 = W2 * W;

  const double vdotF = Fup[1] * vdown[1] + Fup[2] * vdown[2] + Fup[3] * vdown[3];
  const double normF = sqrt(F2);
  const double inormF = (normF > 0 ? 1 / normF : 0);
  const double vdothatf = vdotF * inormF;
  const double vdothatf2 = vdothatf * vdothatf;
  const double hatfx = qpre[1] * inormF; // hatf_i
  const double hatfy = qpre[2] * inormF;
  const double hatfz = qpre[3] * inormF;
  const double hatfupx = Fup[1] * inormF; // hatf^i
  const double hatfupy = Fup[2] * inormF;
  const double hatfupz = Fup[3] * inormF;
  const double e = qpre[0];
  const double eonormF = min(e * inormF, 1.0); // with factor dthin ...

  // drvts of J
  // Eq (65)
  double JdE = W2 + dthin * vdothatf2 * W2
               + (dthick * (3 - 2 * W2) * (-1 + W2)) / (1 + 2 * W2);
  // Eq (69)
  double JdFv = 2 * W2
                * (-1 + (dthin * eonormF * vdothatf)
                   + (2 * dthick * (-1 + W2)) / (1 + 2 * W2));
  // Eq (70)
  double JdFf = (-2 * dthin * eonormF * vdothatf2 * W2);
  // Eq (66)
  double JdFx = JdFv * vup[1] + JdFf * hatfupx;
  double JdFy = JdFv * vup[2] + JdFf * hatfupy;
  double JdFz = JdFv * vup[3] + JdFf * hatfupz;

  // drvts of Hi
  // Eq (71)
  double HdEv = W3 * (-1 - dthin * vdothatf2 + (dthick * (-3 + 2 * W2)) / (1 + 2 * W2));
  // Eq (72)
  double HdEf = -(dthin * vdothatf * W);
  // Eq (67)
  double HxdE = HdEv * vx + HdEf * hatfx;
  double HydE = HdEv * vy + HdEf * hatfy;
  double HzdE = HdEv * vz + HdEf * hatfz;
  // Eq (73)
  double HdFdelta = (1 - dthick * v2 - (dthin * eonormF * vdothatf)) * W;
  // Eq (74)
  double HdFvv = (2 * (1 - dthin * eonormF * vdothatf) * W3)
                 + dthick * W * (2 - 2 * W2 + 1 / (-1 - 2 * W2));
  // Eq (75)
  double HdFff = (2 * dthin * eonormF * vdothatf * W);
  // Eq (76)
  double HdFvf = (2 * dthin * eonormF * vdothatf2 * W3);
  // Eq (77)
  double HdFfv = -(dthin * eonormF * W);
  // Eq (68)
  double HxdFx = HdFdelta + HdFvv * vx * vup[1] + HdFff * hatfx * hatfupx
                 + HdFvf * vx * hatfupx + HdFfv * hatfx * vup[1];
  double HydFx = HdFvv * vy * vup[1] + HdFff * hatfy * hatfupx + HdFvf * vy * hatfupx
                 + HdFfv * hatfy * vup[1];
  double HzdFx = HdFvv * vz * vup[1] + HdFff * hatfz * hatfupx + HdFvf * vz * hatfupx
                 + HdFfv * hatfz * vup[1];

  double HxdFy = HdFvv * vx * vup[2] + HdFff * hatfx * hatfupy + HdFvf * vx * hatfupy
                 + HdFfv * hatfx * vup[2];
  double HydFy = HdFdelta + HdFvv * vy * vup[2] + HdFff * hatfy * hatfupy
                 + HdFvf * vy * hatfupy + HdFfv * hatfy * vup[2];
  double HzdFy = HdFvv * vz * vup[2] + HdFff * hatfz * hatfupy + HdFvf * vz * hatfupy
                 + HdFfv * hatfz * vup[2];

  double HxdFz = HdFvv * vx * vup[3] + HdFff * hatfx * hatfupz + HdFvf * vx * hatfupz
                 + HdFfv * hatfx * vup[3];
  double HydFz = HdFvv * vy * vup[3] + HdFff * hatfy * hatfupz + HdFvf * vy * hatfupz
                 + HdFfv * hatfy * vup[3];
  double HzdFz = HdFdelta + HdFvv * vz * vup[3] + HdFff * hatfz * hatfupz
                 + HdFvf * vz * hatfupz + HdFfv * hatfz * vup[3];

  // Build the Jacobian
  // Eq (61)
  double J00 = -alpW * (kapas - kaps * JdE);
  // Eq (62)
  double J0x = +alpW * kaps * JdFx + alpW * kapas * vup[1];
  double J0y = +alpW * kaps * JdFy + alpW * kapas * vup[2];
  double J0z = +alpW * kaps * JdFz + alpW * kapas * vup[3];
  // Eq (63)
  double Jx0 = -alpha * (kapas * HxdE + W * kapa * vx * JdE);
  double Jy0 = -alpha * (kapas * HydE + W * kapa * vy * JdE);
  double Jz0 = -alpha * (kapas * HzdE + W * kapa * vz * JdE);
  // Eq (64)
  double Jxx = -alpha * (kapas * HxdFx + W * kapa * vx * JdFx);
  double Jxy = -alpha * (kapas * HxdFy + W * kapa * vx * JdFy);
  double Jxz = -alpha * (kapas * HxdFz + W * kapa * vx * JdFz);

  double Jyy = -alpha * (kapas * HydFx + W * kapa * vy * JdFx);
  double Jyx = -alpha * (kapas * HydFy + W * kapa * vy * JdFy);
  double Jyz = -alpha * (kapas * HydFz + W * kapa * vy * JdFz);

  double Jzx = -alpha * (kapas * HzdFx + W * kapa * vz * JdFx);
  double Jzy = -alpha * (kapas * HzdFy + W * kapa * vz * JdFy);
  double Jzz = -alpha * (kapas * HzdFz + W * kapa * vz * JdFz);

  // Store Jacobian into J
  double A_data[4][4] = {
    1 - cdt * J00, -cdt * J0x, -cdt * J0y, -cdt * J0z,    -cdt * Jx0,    1 - cdt * Jxx,
    -cdt * Jxy,    -cdt * Jxz, -cdt * Jy0, -cdt * Jyx,    1 - cdt * Jyy, -cdt * Jyz,
    -cdt * Jz0,    -cdt * Jzx, -cdt * Jzy, 1 - cdt * Jzz,
  };
  for(int a = 0; a < 4; ++a) {
    for(int b = 0; b < 4; ++b) {
      // gsl_matrix_set(J, a, b, A_data[a][b]);
      J->DD[a][b] = A_data[a][b];
    }
  }
}

double dot(double *a, double *b, int length) {
  double result = 0.0;
  for(int i = 0; i < length; i++) {
    result += a[i] * b[i];
  }
  return result;
}

int prepare_sources(gsl_vector const * q, ghl_m1_powell_params * p) {
    assemble_rT_lab_frame(p->n4D,p->E,p->F4,p->P4,p->rT4DD);

  p->J = calc_J_from_rT(p->u4U, p->proj4, p->rT4DD);
  calc_H4D_from_rT(p->u4U, p->proj4, p->rT4DD, &p->H4);

  calc_rad_sources(p->eta, p->kabs, p->kscat, p->u4U, p->J, p->H4, &p->S4);

  p->Edot = calc_rE_source(p->metric, p->S4);
  calc_rF_source(p->metric, p->adm_aux, p->S4, &p->rF_source);
  p->rE_source = calc_rE_source(p->metric, p->S4);

  return GSL_SUCCESS;
}

int prepare(const gsl_vector *q, ghl_m1_powell_params *p) {
  // Johnny: removed chi passing around to just p->chi
  //         this can lead to error here if prepare_closure fails
  int ierr = prepare_closure(q, p);
  if(ierr != GSL_SUCCESS) {
    return ierr;
  }
  return GSL_SUCCESS;
}

// Jacobian of the implicit function
int impl_func_jac(gsl_vector const * q, void * params, gsl_matrix * J) {
    ghl_m1_powell_params * p = (ghl_m1_powell_params *)params;
    // Johnny: maybe it is best to separate param struct used in powell and in the large source func
    int ierr = prepare(q, p);
    if (ierr != GSL_SUCCESS) {
        return ierr;
    }

#define EVALUATE_ZJAC                                                               \
  double m_q[] = { p->E, p->F4->D[1], p->F4->D[2], p->F4->D[3] };                      \
  double m_Fup[] = { p->F4->U[0], p->F4->U[1], p->F4->U[2], p->F4->U[3] };              \
  double m_F2 = dot(p->F4->D, p->F4->U);                                              \
  double m_chi = p->chi;                                                            \
  double m_kscat = p->kscat;                                                        \
  double m_kabs = p->kabs;                                                          \
  double m_vup[] = { p->v->U[0], p->v->U[1], p->v->U[2], p->v->U[3] };              \
  double m_vdw[] = { p->v->D[0], p->v->D[1], p->v->D[2], p->v->D[3] };              \
  double m_v2 = dot(p->v->U, p->v->D);                                              \
  double m_W = p->W;                                                                \
  double m_alpha = p->alp;                                                          \
  double m_cdt = p->cdt;                                                            \
  double m_qstar[] = { p->E_star, p->F4_star->D[1], p->F4_star->D[2], p->F4_star->D[3] };  \
                                                                                    \
  __source_jacobian_low_level(                                                      \
        m_q, m_Fup, m_F2, m_chi, m_kscat, m_kabs, m_vup, m_vdw, m_v2, m_W, m_alpha, \
        m_cdt, m_qstar, J);

  EVALUATE_ZJAC

  return GSL_SUCCESS;
}

// Function to rootfind for
//    f(q) = q - q^* - dt S[q]
int impl_func_val(const gsl_vector *q, void *params, gsl_vector *f) {
  ghl_m1_powell_params *p = (ghl_m1_powell_params *)params;
  int ierr = prepare(q, p);
  if(ierr != GSL_SUCCESS) {
    return ierr;
  }
#define EVALUATE_ZFUNC \
    gsl_vector_set(f, 0, gsl_vector_get(q, 0) - p->E_star        - p->cdt * p->rE_source); \
    gsl_vector_set(f, 1, gsl_vector_get(q, 1) - p->F4_star->D[1] - p->cdt * p->rF_source->D[1]); \
    gsl_vector_set(f, 2, gsl_vector_get(q, 2) - p->F4_star->D[2] - p->cdt * p->rF_source->D[2]); \
    gsl_vector_set(f, 3, gsl_vector_get(q, 3) - p->F4_star->D[3] - p->cdt * p->rF_source->D[3]);

  EVALUATE_ZFUNC

  return GSL_SUCCESS;
}

// Function and Jacobian evaluation
int impl_func_val_jac(gsl_vector const * q, void * params, gsl_vector * f, gsl_matrix * J) {
    ghl_m1_powell_params * p = (ghl_m1_powell_params *)params;
    int ierr = prepare(q, p);
    if (ierr != 0) {
        return ierr;
    }
  EVALUATE_ZFUNC
  EVALUATE_ZJAC

  return GSL_SUCCESS;
}
#undef EVALUATE_ZFUNC
#undef EVALUATE_ZJAC

// Explicit update
//  F^k+1 = F^* + dt   [ A[F^m] + S[F^k+1] ]
void explicit_update(
        ghl_m1_powell_params * p,
        double * E_new,
        ghl_radiation_flux_vector * F4_new) {
    double alp = p->metric->lapse;
    double n4U[4] = {p->metric->lapseinv, 
                    -p->metric->betaU[0]*p->metric->lapseinv,
                    -p->metric->betaU[1]*p->metric->lapseinv,
                    -p->metric->betaU[2]*p->metric->lapseinv};
    *E_new      = p->E_star       + p->cdt * p->rE_source;
    
    // rF_source only has U, need to lower indicies
    double rF_source_D[4] = {0.0,0.0,0.0,0.0};
    for(int a = 0; a < 4; ++a) {
      for (int b = 0; b < 4; ++b) {
        rF_source_D[a] = p->rF_source->U[b] * p->adm_aux->g4DD[a][b];
      }
    }
    F4_new->D[1] = p->F4_star->D[1] + p->cdt * rF_source_D[1];
    F4_new->D[2] = p->F4_star->D[2] + p->cdt * rF_source_D[2];
    F4_new->D[3] = p->F4_star->D[3] + p->cdt * rF_source_D[3];
    // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
    F4_new->D[0] = - alp * p->n4U[1] * F4_new->D[1]
                   - alp * p->n4U[2] * F4_new->D[2]
                   - alp * p->n4U[3] * F4_new->D[3];
} //explicit_update

int source_update(
      const double cdt,
      ghl_m1_closure_t closure,
      gsl_multiroot_fdfsolver *gsl_solver_nd,
      const double E_old,
      double *E_new,
      const ghl_radiation_flux_vector *F4_old,
      ghl_radiation_flux_vector *F4_new,
      const ghl_m1_thc_params *thc_params,
      ghl_m1_powell_params *p) {

  gsl_multiroot_function_fdf zfunc
        = { impl_func_val, impl_func_jac, impl_func_val_jac, 4, &p };

  double alp = p->metric->lapse;
  double n4U[4] = { p->metric->lapseinv, 
                    -p->metric->betaU[0] * p->metric->lapseinv,
                    -p->metric->betaU[1] * p->metric->lapseinv,
                    -p->metric->betaU[2] * p->metric->lapseinv };

  // Old solution
  double qold[4] = { E_old, F4_old->D[1], F4_old->D[2], F4_old->D[3] };
  gsl_vector_view xold = gsl_vector_view_array(qold, 4);

  // Non stiff limit, use explicit update
  if(cdt * p->kabs < 1 && cdt * p->kscat < 1) {
    prepare(&xold.vector, &p);
    explicit_update(&p, E_new, F4_new);

    double q[4] = { *E_new, F4_new->D[1], F4_new->D[2], F4_new->D[3] };
    gsl_vector_view x = gsl_vector_view_array(q, 4);
    prepare_closure(&x.vector, &p);
    return -1;
    // return THC_M1_SOURCE_THIN; // TODO
  }

  // Our scheme cannot capture this dynamics (tau << dt), so we go
  // directly to the equilibrium
  if (thc_params->source_thick_limit > 0 &&
          SQ(cdt)*(p->kabs*(p->kabs + p->kscat)) > SQ(thc_params->source_thick_limit)) {
      return -2; //TODO
      //return THC_M1_SOURCE_EQUIL; //TODO
  }
  // This handles the scattering dominated limit
  if (thc_params->source_scat_limit > 0 && cdt*p->kscat > thc_params->source_scat_limit) {
      return -3; //TODO
      //return THC_M1_SOURCE_SCAT; //TODO
  }

  // Initial guess for the solution
  double q[4] = { *E_new, F4_new->D[1], F4_new->D[2], F4_new->D[3] };
  gsl_vector_view x = gsl_vector_view_array(q, 4);

  int ierr = gsl_multiroot_fdfsolver_set(gsl_solver_nd, &zfunc, &x.vector);
  int iter = 0;
  do {
    if(iter < thc_params->source_maxiter) {
      ierr = gsl_multiroot_fdfsolver_iterate(gsl_solver_nd);
      iter++;
    }
    // TODO: removed error handling for now
    ierr = gsl_multiroot_test_delta(
          gsl_solver_nd->dx, gsl_solver_nd->x, thc_params->source_epsabs, thc_params->source_epsrel);
  } while(ierr == GSL_CONTINUE);

  *E_new = gsl_vector_get(gsl_solver_nd->x, 0);
  F4_new->D[1] = gsl_vector_get(gsl_solver_nd->x, 1);
  F4_new->D[2] = gsl_vector_get(gsl_solver_nd->x, 2);
  F4_new->D[3] = gsl_vector_get(gsl_solver_nd->x, 3);
  // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  F4_new->D[0] = - alp * n4U[1] * F4_new->D[1]
                 - alp * n4U[2] * F4_new->D[2]
                 - alp * n4U[3] * F4_new->D[3];

  prepare_closure(gsl_solver_nd->x, &p);
  return -4; //TODO
  //return THC_M1_SOURCE_OK; //TODO



} // source_update
