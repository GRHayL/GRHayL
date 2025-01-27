#include "ghl.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

// gsl_multiroots: https://www.gnu.org/software/gsl/doc/html/multiroots.html

// TODO: move this to ghl_radiation.h
typedef struct {
  ghl_m1_closure_t closure;
  gsl_root_fsolver *gsl_solver_1d;
  double cdt;
  ghl_metric_quantities *metric;
  ghl_ADM_aux_quantities *adm_aux;
  ghl_radiation_metric_tensor *proj4;
  ghl_primitive_quantities *prims;

  // remove these
  double *u4D;
  double *u4U;
  double *n4D;
  double *n4U;

  double J;
  ghl_radiation_flux_vector *F4;
  ghl_radiation_flux_vector *H4;

  ghl_radiation_pressure_tensor *P4;
  ghl_stress_energy *rT4DD;
  ghl_radiation_con_source_vector *S4;

  ghl_radiation_con_source_vector *F_src;
  double W;
  double Estar;
  ghl_radiation_flux_vector *Fstar;
  double Edot;
  double chi;
  double eta;
  double kabs;
  double kscat;

  double E;
} Params;

void init_params(
      Params *p,
      const ghl_m1_closure_t closure,
      const gsl_root_fsolver *gsl_solver_1d,
      const double cdt,
      const ghl_metric_quantities *metric,
      const ghl_ADM_aux_quantities *adm_aux,
      const ghl_radiation_metric_tensor *proj4,
      const ghl_primitive_quantities *prims,
      const double *u4D,
      const double *u4U,
      const double *n4D,
      const double *n4U,
      const double J,
      const ghl_radiation_flux_vector *F4,
      const ghl_radiation_flux_vector *H4,
      const ghl_radiation_pressure_tensor *P4,
      const ghl_radiation_con_source_vector *F_src,
      const double W,
      const double Estar,
      const ghl_radiation_flux_vector *Fstar,
      const double Edot,
      const ghl_stress_energy *rT4DD,
      const ghl_radiation_con_source_vector *S4,
      const double chi,
      const double eta,
      const double kabs,
      const double kscat) {
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
  p->P4 = P4;
  p->F4 = F4;
  p->H4 = H4;
  p->F_src = F_src;
  p->W = W;
  p->Estar = Estar;
  p->Fstar = Fstar;
  p->Edot = Edot;
  p->chi = chi;
  p->eta = eta;
  p->kabs = kabs;
  p->kscat = kscat;
  p->E = 0.0; // Initialize E to a default value, if needed
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

int prepare_sources(const gsl_vector *q, Params *p) {
  assemble_rT_lab_frame(p->n4D, p->E, p->F4, p->P4, p->rT4DD);

  p->J = calc_J_from_rT(p->u4U, p->proj4, p->rT4DD);
  calc_H4D_from_rT(p->u4U, p->proj4, p->rT4DD, &p->H4);

  calc_rad_sources(p->eta, p->kabs, p->kscat, p->u4U, p->J, p->H4, &p->S4);

  p->Edot = calc_rE_source(p->metric, p->S4);
  calc_rF_source(p->metric, p->adm_aux, p->S4, &p->F_src);

  return GSL_SUCCESS;
}

int prepare(const gsl_vector *q, Params *p) {
  int ierr = prepare_closure(q, p);
  if(ierr != GSL_SUCCESS) {
    return ierr;
  }

  ierr = prepare_sources(q, p);
  if(ierr != GSL_SUCCESS) {
    return ierr;
  }

  return GSL_SUCCESS;
}

// Jacobian of the implicit function
int impl_func_jac(const gsl_vector *q, void *params, gsl_matrix *J) {
  Params *p = (Params *)params;
  int ierr = prepare(q, p);
  if(ierr != GSL_SUCCESS) {
    return ierr;
  }

#define EVALUATE_ZJAC                                                               \
  double m_q[] = { p->E, p->F->D[1], p->F->D[2], p->F->D[3] };                      \
  double m_Fup[] = { p->F->U[0], p->F->U[1], p->F->U[2], p->F->U[3] };              \
  double m_F2 = dot(p->F->D, p->F->U);                                              \
  double m_chi = p->chi;                                                            \
  double m_kscat = p->kscat;                                                        \
  double m_kabs = p->kabs;                                                          \
  double m_vup[] = { p->v->U[0], p->v->U[1], p->v->U[2], p->v->U[3] };              \
  double m_vdw[] = { p->v->D[0], p->v->D[1], p->v->D[2], p->v->D[3] };              \
  double m_v2 = dot(p->v->U, p->v->D);                                              \
  double m_W = p->W;                                                                \
  double m_alpha = p->alp;                                                          \
  double m_cdt = p->cdt;                                                            \
  double m_qstar[] = { p->Estar, p->Fstar->D[1], p->Fstar->D[2], p->Fstar->D[3] };  \
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
  Params *p = (Params *)params;
  int ierr = prepare(q, p);
  if(ierr != GSL_SUCCESS) {
    return ierr;
  }

#define EVALUATE_ZFUNC                                                                \
  gsl_vector_set(f, 0, gsl_vector_get(q, 0) - p->Estar - p->cdt * p->Edot);           \
  gsl_vector_set(f, 1, gsl_vector_get(q, 1) - p->Fstar->D[1] - p->cdt * p->tS->D[1]); \
  gsl_vector_set(f, 2, gsl_vector_get(q, 2) - p->Fstar->D[2] - p->cdt * p->tS->D[2]); \
  gsl_vector_set(f, 3, gsl_vector_get(q, 3) - p->Fstar->D[3] - p->cdt * p->tS->D[3]);

  EVALUATE_ZFUNC

  return GSL_SUCCESS;
}

// Function and Jacobian evaluation
int impl_func_val_jac(const gsl_vector *q, void *params, gsl_vector *f, gsl_matrix *J) {
  Params *p = (Params *)params;
  int ierr = prepare(q, p);
  if(ierr != 0) {
    return ierr;
  }

  EVALUATE_ZFUNC
  EVALUATE_ZJAC

  return GSL_SUCCESS;
}

#undef EVALUATE_ZFUNC
#undef EVALUATE_ZJAC

void explicit_update(Params *p, double *Enew, ghl_radiation_flux_vector *Fnew) {
  double alp = p->metric->lapse;
  double n4U[4] = { p->metric->lapseinv, -p->metric->betaU[0] * p->metric->lapseinv,
                    -p->metric->betaU[1] * p->metric->lapseinv,
                    -p->metric->betaU[2] * p->metric->lapseinv };
  *Enew = p->Estar + p->cdt * p->Edot;
  Fnew->D[1] = p->Fstar->D[1] + p->cdt * p->S4->D[1];
  Fnew->D[2] = p->Fstar->D[2] + p->cdt * p->S4->D[2];
  Fnew->D[3] = p->Fstar->D[3] + p->cdt * p->S4->D[3];
  // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  Fnew->D[0] = -alp * p->n4U[1] * Fnew->D[1] - alp * p->n4U[2] * Fnew->D[2]
               - alp * p->n4U[3] * Fnew->D[3];
} // explicit_update

// Function to evaluate f(q)
int impl_func_val(const gsl_vector *q, void *params, gsl_vector *f) {
  // Cast the generic pointer to Params structure
  Params *p = (Params *)params;

  // Call the preparation function
  int ierr = prepare(q, p);
  if(ierr != GSL_SUCCESS) {
    return ierr; // Propagate the error if preparation fails
  }

  // Define the root-finding equations
  gsl_vector_set(f, 0, gsl_vector_get(q, 0) - p->Estar - p->cdt * p->F_src->D[0]);
  gsl_vector_set(f, 1, gsl_vector_get(q, 1) - p->Fstar->D[1] - p->cdt * p->F_src->D[1]);
  gsl_vector_set(f, 2, gsl_vector_get(q, 2) - p->Fstar->D[2] - p->cdt * p->F_src->D[2]);
  gsl_vector_set(f, 3, gsl_vector_get(q, 3) - p->Fstar->D[3] - p->cdt * p->F_src->D[3]);

  return GSL_SUCCESS; // Indicate success
} // impt_func_val

int source_update(
      ghl_m1_closure_t closure,
      gsl_multiroot_fdfsolver *gsl_solver_nd,
      const double eta,
      const double kabs,
      const double kscat,
      const double Eold,
      const double cdt,
      double *Enew,
      const ghl_radiation_flux_vector *Fold,
      ghl_radiation_flux_vector *Fnew,
      double *chi,
      double source_thick_limit,
      double source_scat_limit,
      int source_maxiter,
      double source_epsabs,
      double source_epsrel,
      Params *p) {

  gsl_multiroot_function_fdf zfunc
        = { impl_func_val, impl_func_jac, impl_func_val_jac, 4, &p };

  double alp = p->metric->lapse;
  double n4U[4] = { p->metric->lapseinv, -p->metric->betaU[0] * p->metric->lapseinv,
                    -p->metric->betaU[1] * p->metric->lapseinv,
                    -p->metric->betaU[2] * p->metric->lapseinv };

  // Old solution
  double qold[4] = { Eold, Fold->D[1], Fold->D[2], Fold->D[3] };
  gsl_vector_view xold = gsl_vector_view_array(qold, 4);

  // Non stiff limit, use explicit update
  if(cdt * kabs < 1 && cdt * kscat < 1) {
    prepare(&xold.vector, &p);
    explicit_update(&p, Enew, Fnew);

    double q[4] = { *Enew, Fnew->D[1], Fnew->D[2], Fnew->D[3] };
    gsl_vector_view x = gsl_vector_view_array(q, 4);
    prepare_closure(&x.vector, &p);
    *chi = p.chi;

    return THC_M1_SOURCE_THIN; // TODO
  }

  // Our scheme cannot capture this dynamics (tau << dt), so we go
  // directly to the equilibrium
  if(source_thick_limit > 0
     && SQ(cdt) * (kabs * (kabs + kscat)) > SQ(source_thick_limit)) {
    return THC_M1_SOURCE_EQUIL; // TODO
  }

  // This handles the scattering dominated limit
  if(source_scat_limit > 0 && cdt * kscat > source_scat_limit) {
    return THC_M1_SOURCE_SCAT; // TODO
  }

  // Initial guess for the solution
  double q[4] = { *Enew, Fnew->D[1], Fnew->D[2], Fnew->D[3] };
  gsl_vector_view x = gsl_vector_view_array(q, 4);

  int ierr = gsl_multiroot_fdfsolver_set(gsl_solver_nd, &zfunc, &x.vector);
  int iter = 0;
  do {
    if(iter < source_maxiter) {
      ierr = gsl_multiroot_fdfsolver_iterate(gsl_solver_nd);
      iter++;
    }
    // TODO: removed error handling for now
    ierr = gsl_multiroot_test_delta(
          gsl_solver_nd->dx, gsl_solver_nd->x, source_epsabs, source_epsrel);
  } while(ierr == GSL_CONTINUE);

  *Enew = gsl_vector_get(gsl_solver_nd->x, 0);
  Fnew->D[1] = gsl_vector_get(gsl_solver_nd->x, 1);
  Fnew->D[2] = gsl_vector_get(gsl_solver_nd->x, 2);
  Fnew->D[3] = gsl_vector_get(gsl_solver_nd->x, 3);
  // F_0 = g_0i F^i = beta_i F^i = beta^i F_i
  Fnew->D[0] = -alp * n4U[1] * Fnew->D[1] - alp * n4U[2] * Fnew->D[2]
               - alp * n4U[3] * Fnew->D[3];

  prepare_closure(gsl_solver_nd->x, &p);
  *chi = p.chi;

  return THC_M1_SOURCE_OK; // TODO

} // source_update
