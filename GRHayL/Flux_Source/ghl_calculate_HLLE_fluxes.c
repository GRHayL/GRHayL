#include "ghl_flux_source.h"
/*
 *  * Compute the HLLE-derived fluxes on the left face for all components.
 *   */

/* NRPYSTART
import sympy as sp
import nrpy.indexedexp as ixp
import nrpy.c_codegen as ccg
import ghl_generate_fluxes as ghl

S_flux_rhs = ixp.zerorank1()
rho_flux_rhs, tau_flux_rhs, S_flux_rhs, Ye_flux_rhs, ent_flux_rhs = ghl.calculate_fluxes_rhs()


with open("fluxes_hybrid.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2]], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]"]))

with open("fluxes_hybrid_entropy.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2], ent_flux_rhs], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]", "cons_rhs->entropy"]))

with open("fluxes_tabulated.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2], Ye_flux_rhs], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]", "cons_rhs->Y_e"]))

with open("fluxes_tabulated_entropy.h", "w") as file:
    file.write(ccg.c_codegen([rho_flux_rhs, tau_flux_rhs, S_flux_rhs[0], S_flux_rhs[1], S_flux_rhs[2], Ye_flux_rhs, ent_flux_rhs], ["cons_rhs->rho", "cons_rhs->tau", "cons_rhs->SD[0]", "cons_rhs->SD[1]", "cons_rhs->SD[2]", "cons_rhs->Y_e", "cons_rhs->entropy"]))
NRPYEND */

void ghl_calculate_HLLE_fluxes_hybrid2(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

  double smallb_zeror;
  double smallb_zerol;
  {
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  const double sqrt_detg_over_csum = sqrt_detg/(cmax + cmin);
  const double cmul = cmax*cmin;

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);
  const double betaD[3] = {metric_aux.g4DD[0][1],
                           metric_aux.g4DD[0][2],
                           metric_aux.g4DD[0][3]};

  double vlD[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaDD, prims_l->vU, vlD);
  const double Bvtilde_l = prims_l->u0*(prims_l->BU[0]*(vlD[0] + betaD[0])
                                      + prims_l->BU[1]*(vlD[1] + betaD[1])
                                      + prims_l->BU[2]*(vlD[2] + betaD[2]));

  const double tmp_4vec_lU[4] = {ADM_metric->lapseinv*Bvtilde_l,
                                ADM_metric->lapseinv*(prims_l->BU[0]/prims_l->u0 + Bvtilde_l*prims_l->vU[0]),
                                ADM_metric->lapseinv*(prims_l->BU[1]/prims_l->u0 + Bvtilde_l*prims_l->vU[1]),
                                ADM_metric->lapseinv*(prims_l->BU[2]/prims_l->u0 + Bvtilde_l*prims_l->vU[2])};
  const double tmp_vec2_l = ghl_compute_vec2_from_vec4D( metric_aux.g4DD, tmp_4vec_lU);

  const double tmp_1 = (tmp_vec2_l + h_l*prims_l->rho)*prims_l->u0*prims_l->u0;
  const double tmp_2 = prims_l->press + tmp_vec2_l/2.0;

  // These variables represent various elements of the solution, some of which
  // change w.r.t. the direction of the flux.
  const double soln_scal_l    =  -tmp_4vec_lU[0]*tmp_4vec_lU[0] + tmp_1   /* u^0 */      /* u^0 */    + tmp_2*metric_aux.g4UU[0][0];
  const double soln_vec1_l[3] = {-tmp_4vec_lU[0]*tmp_4vec_lU[1] + tmp_1*prims_l->vU[0]   /* u^0 */    + tmp_2*metric_aux.g4UU[0][1],
                                 -tmp_4vec_lU[0]*tmp_4vec_lU[2] + tmp_1*prims_l->vU[1]   /* u^0 */    + tmp_2*metric_aux.g4UU[0][2],
                                 -tmp_4vec_lU[0]*tmp_4vec_lU[3] + tmp_1*prims_l->vU[2]   /* u^0 */    + tmp_2*metric_aux.g4UU[0][3]};
  const double soln_vec2_l[3] = {-tmp_4vec_lU[direction+1]*tmp_4vec_lU[1] + tmp_1*prims_l->vU[direction]*prims_l->vU[0] + tmp_2*metric_aux.g4UU[direction+1][1],
                                 -tmp_4vec_lU[direction+1]*tmp_4vec_lU[2] + tmp_1*prims_l->vU[direction]*prims_l->vU[1] + tmp_2*metric_aux.g4UU[direction+1][2],
                                 -tmp_4vec_lU[direction+1]*tmp_4vec_lU[3] + tmp_1*prims_l->vU[direction]*prims_l->vU[2] + tmp_2*metric_aux.g4UU[direction+1][3]};

  double vrD[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaDD, prims_r->vU, vrD);
  const double Bvtilde_r = prims_r->u0*(prims_r->BU[0]*(vrD[0] + betaD[0])
                                      + prims_r->BU[1]*(vrD[1] + betaD[1])
                                      + prims_r->BU[2]*(vrD[2] + betaD[2]));

  const double tmp_4vec_rU[4] = {ADM_metric->lapseinv*Bvtilde_r,
                                ADM_metric->lapseinv*(prims_r->BU[0]/prims_r->u0 + Bvtilde_r*prims_r->vU[0]),
                                ADM_metric->lapseinv*(prims_r->BU[1]/prims_r->u0 + Bvtilde_r*prims_r->vU[1]),
                                ADM_metric->lapseinv*(prims_r->BU[2]/prims_r->u0 + Bvtilde_r*prims_r->vU[2])};
  const double tmp_vec2_r = ghl_compute_vec2_from_vec4D( metric_aux.g4DD, tmp_4vec_rU);

  const double tmp_3 = (tmp_vec2_r + h_r*prims_r->rho)*prims_r->u0*prims_r->u0;
  const double tmp_4 = prims_r->press + tmp_vec2_r/2.0;

  const double soln_scal_r    =  -tmp_4vec_rU[0]*tmp_4vec_rU[0] + tmp_3   /* u^0 */      /* u^0 */    + tmp_4*metric_aux.g4UU[0][0];
  const double soln_vec1_r[3] = {-tmp_4vec_rU[0]*tmp_4vec_rU[1] + tmp_3*prims_r->vU[0]   /* u^0 */    + tmp_4*metric_aux.g4UU[0][1],
                                 -tmp_4vec_rU[0]*tmp_4vec_rU[2] + tmp_3*prims_r->vU[1]   /* u^0 */    + tmp_4*metric_aux.g4UU[0][2],
                                 -tmp_4vec_rU[0]*tmp_4vec_rU[3] + tmp_3*prims_r->vU[2]   /* u^0 */    + tmp_4*metric_aux.g4UU[0][3]};
  const double soln_vec2_r[3] = {-tmp_4vec_rU[direction+1]*tmp_4vec_rU[1] + tmp_3*prims_r->vU[direction]*prims_r->vU[0] + tmp_4*metric_aux.g4UU[direction+1][1],
                                 -tmp_4vec_rU[direction+1]*tmp_4vec_rU[2] + tmp_3*prims_r->vU[direction]*prims_r->vU[1] + tmp_4*metric_aux.g4UU[direction+1][2],
                                 -tmp_4vec_rU[direction+1]*tmp_4vec_rU[3] + tmp_3*prims_r->vU[direction]*prims_r->vU[2] + tmp_4*metric_aux.g4UU[direction+1][3]};

  const double soln_scal_r2    =  0;//metric_aux.g4UU[0][0]*tmp_vec2_r/2.0;
  const double soln_scal_l2    =  0;//metric_aux.g4UU[0][0]*tmp_vec2_l/2.0;
  cons_rhs->tau = ADM_metric->lapse*sqrt_detg_over_csum*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] - cmul*(soln_scal_r2 - soln_scal_l2));
  //cons_rhs->tau = ADM_metric->lapse*sqrt_detg_over_csum*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] - cmul*(soln_scal_r - soln_scal_l)) - cons_rhs->rho;

  cons_rhs->SD[0] = sqrt_detg_over_csum*(metric_aux.g4DD[1][1]*(cmax*soln_vec2_l[0] + cmin*soln_vec2_r[0] + cmul*(soln_vec1_l[0] - soln_vec1_r[0]))
                                   + metric_aux.g4DD[1][2]*(cmax*soln_vec2_l[1] + cmin*soln_vec2_r[1] + cmul*(soln_vec1_l[1] - soln_vec1_r[1]))
                                   + metric_aux.g4DD[1][3]*(cmax*soln_vec2_l[2] + cmin*soln_vec2_r[2] + cmul*(soln_vec1_l[2] - soln_vec1_r[2]))
                                   + betaD[0]*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] + cmul*(soln_scal_l - soln_scal_r)));

  cons_rhs->SD[1] = sqrt_detg_over_csum*(metric_aux.g4DD[1][2]*(cmax*soln_vec2_l[0] + cmin*soln_vec2_r[0] + cmul*(soln_vec1_l[0] - soln_vec1_r[0]))
                                   + metric_aux.g4DD[2][2]*(cmax*soln_vec2_l[1] + cmin*soln_vec2_r[1] + cmul*(soln_vec1_l[1] - soln_vec1_r[1]))
                                   + metric_aux.g4DD[2][3]*(cmax*soln_vec2_l[2] + cmin*soln_vec2_r[2] + cmul*(soln_vec1_l[2] - soln_vec1_r[2]))
                                   + betaD[1]*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] + cmul*(soln_scal_l - soln_scal_r)));

  cons_rhs->SD[2] = sqrt_detg_over_csum*(metric_aux.g4DD[1][3]*(cmax*soln_vec2_l[0] + cmin*soln_vec2_r[0] + cmul*(soln_vec1_l[0] - soln_vec1_r[0]))
                                   + metric_aux.g4DD[2][3]*(cmax*soln_vec2_l[1] + cmin*soln_vec2_r[1] + cmul*(soln_vec1_l[1] - soln_vec1_r[1]))
                                   + metric_aux.g4DD[3][3]*(cmax*soln_vec2_l[2] + cmin*soln_vec2_r[2] + cmul*(soln_vec1_l[2] - soln_vec1_r[2]))
                                   + betaD[2]*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] + cmul*(soln_scal_l - soln_scal_r)));
smallb_zeror = Bvtilde_r;
smallb_zerol = Bvtilde_l;
  }

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

const double g4DD00 = metric_aux.g4DD[0][0];
const double g4DD01 = metric_aux.g4DD[0][1];
const double g4DD02 = metric_aux.g4DD[0][2];
const double g4DD03 = metric_aux.g4DD[0][3];
const double g4DD11 = metric_aux.g4DD[1][1];
const double g4DD12 = metric_aux.g4DD[1][2];
const double g4DD13 = metric_aux.g4DD[1][3];
const double g4DD22 = metric_aux.g4DD[2][2];
const double g4DD23 = metric_aux.g4DD[3][3];
const double g4DD33 = metric_aux.g4DD[3][3];

const double g4UU00 = metric_aux.g4UU[0][0];
const double g4UU01 = metric_aux.g4UU[0][1];
const double g4UU02 = metric_aux.g4UU[0][2];
const double g4UU03 = metric_aux.g4UU[0][3];
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;

const double tmp0 = (1.0/(cmax + cmin));
const double tmp1 = prims_l->rho*prims_l->u0*sqrt_detg;
const double tmp3 = prims_r->rho*prims_r->u0*sqrt_detg;
const double tmp5 = cmax*cmin;
const double tmp7 = prims_l->u0*prims_l->vU[0];
const double tmp8 = prims_l->u0*prims_l->vU[1];
const double tmp9 = prims_l->u0*prims_l->vU[2];
const double tmp18 = (1.0/(prims_l->u0));
const double tmp33 = ADM_metric->lapse*sqrt_detg;
const double tmp34 = prims_r->u0*prims_r->vU[0];
const double tmp35 = prims_r->u0*prims_r->vU[1];
const double tmp36 = prims_r->u0*prims_r->vU[2];
const double tmp41 = (1.0/(prims_r->u0));
const double tmp58 = g4DD01*metric_aux.g4UU[0][direction+1] + g4DD11*metric_aux.g4UU[1][direction+1] + g4DD12*metric_aux.g4UU[2][direction+1] + g4DD13*metric_aux.g4UU[3][direction+1];
const double tmp61 = cmax*sqrt_detg;
const double tmp63 = cmin*sqrt_detg;
const double tmp70 = g4DD02*metric_aux.g4UU[0][direction+1] + g4DD12*metric_aux.g4UU[1][direction+1] + g4DD22*metric_aux.g4UU[2][direction+1] + g4DD23*metric_aux.g4UU[3][direction+1];
const double tmp71 = g4DD03*metric_aux.g4UU[0][direction+1] + g4DD13*metric_aux.g4UU[1][direction+1] + g4DD23*metric_aux.g4UU[2][direction+1] + g4DD33*metric_aux.g4UU[3][direction+1];
const double tmp10 = g4DD01*prims_l->u0 + g4DD11*tmp7 + g4DD12*tmp8 + g4DD13*tmp9;
const double tmp13 = g4DD02*prims_l->u0 + g4DD12*tmp7 + g4DD22*tmp8 + g4DD23*tmp9;
  //const double Bvtilde_r = prims_r->u0*(prims_r->BU[0]*(vrD[0] + betaD[0])
  //                                    + prims_r->BU[1]*(vrD[1] + betaD[1])
  //                                    + prims_r->BU[2]*(vrD[2] + betaD[2]));
const double tmp40 = smallb_zeror*ADM_metric->lapseinv;
// TODO: check if there's a missing u0 in these
//ADM_metric->lapseinv*(prims_r->BU[0]*(g4DD01*prims_r->u0 + g4DD11*prims_r->vU[0] + g4DD12*prims_r->vU[1] + g4DD13*prims_r->vU[2])
//                    + prims_r->BU[1]*(g4DD02*prims_r->u0 + g4DD12*prims_r->vU[0] + g4DD22*prims_r->vU[1] + g4DD23*prims_r->vU[2])
//                    + prims_r->BU[2]*(g4DD03*prims_r->u0 + g4DD13*prims_r->vU[0] + g4DD23*prims_r->vU[1] + g4DD33*prims_r->vU[2]));
const double tmp59 = g4DD00*prims_l->u0 + g4DD01*prims_l->u0*prims_l->vU[0] + g4DD02*prims_l->u0*prims_l->vU[1] + g4DD03*prims_l->u0*prims_l->vU[2];
const double tmp17 = smallb_zerol*ADM_metric->lapseinv;
//const double tmp17 = ADM_metric->lapseinv*(prims_l->BU[0]*tmp10
//                                         + prims_l->BU[1]*tmp13
//                                         + prims_l->BU[2]*(g4DD03*prims_l->u0 + g4DD13*tmp7 + g4DD23*tmp8 + g4DD33*tmp9));
const double tmp42 = ADM_metric->lapseinv*prims_r->BU[direction]*tmp41 + prims_r->vU[direction]*tmp40;
const double tmp43 = tmp41*(ADM_metric->lapseinv*prims_r->BU[0] + tmp34*tmp40);
const double tmp44 = tmp41*(ADM_metric->lapseinv*prims_r->BU[1] + tmp35*tmp40);
const double tmp45 = tmp41*(ADM_metric->lapseinv*prims_r->BU[2] + tmp36*tmp40);
const double tmp19 = ADM_metric->lapseinv*prims_l->BU[direction]*tmp18 + prims_l->vU[direction]*tmp17;
const double tmp20 = tmp18*(ADM_metric->lapseinv*prims_l->BU[0] + tmp17*tmp7);
const double tmp21 = tmp18*(ADM_metric->lapseinv*prims_l->BU[1] + tmp17*tmp8);
const double tmp22 = tmp18*(ADM_metric->lapseinv*prims_l->BU[2] + tmp17*tmp9);
const double tmp46 = g4DD00*tmp40 + g4DD01*tmp43 + g4DD02*tmp44 + g4DD03*tmp45;
const double tmp48 = g4DD01*tmp40 + g4DD11*tmp43 + g4DD12*tmp44 + g4DD13*tmp45;
const double tmp50 = g4DD02*tmp40 + g4DD12*tmp43 + g4DD22*tmp44 + g4DD23*tmp45;
const double tmp52 = tmp45*(g4DD03*tmp40 + g4DD13*tmp43 + g4DD23*tmp44 + g4DD33*tmp45);
const double tmp23 = g4DD00*tmp17 + g4DD01*tmp20 + g4DD02*tmp21 + g4DD03*tmp22;
const double tmp25 = g4DD01*tmp17 + g4DD11*tmp20 + g4DD12*tmp21 + g4DD13*tmp22;
const double tmp27 = g4DD02*tmp17 + g4DD12*tmp20 + g4DD22*tmp21 + g4DD23*tmp22;
const double tmp29 = tmp22*(g4DD03*tmp17 + g4DD13*tmp20 + g4DD23*tmp21 + g4DD33*tmp22);

const double tmp_4vec_rU[4] = {ADM_metric->lapseinv*smallb_zeror,
                              ADM_metric->lapseinv*(prims_r->BU[0]/prims_r->u0 + smallb_zeror*prims_r->vU[0]),
                              ADM_metric->lapseinv*(prims_r->BU[1]/prims_r->u0 + smallb_zeror*prims_r->vU[1]),
                              ADM_metric->lapseinv*(prims_r->BU[2]/prims_r->u0 + smallb_zeror*prims_r->vU[2])};
const double tmp_vec2_r = ghl_compute_vec2_from_vec4D( metric_aux.g4DD, tmp_4vec_rU);
// TODO: see why commented terms do not pass while tmp_vec2_r/2.0 does; similarly for tmp30 below
const double tmp53 = prims_r->press + tmp_vec2_r/2.0;//+ 0.5*tmp45*tmp46 + 0.5*tmp45*tmp48 + 0.5*tmp45*tmp50 + 0.5*tmp52;
const double tmp54 = h_r*prims_r->rho + tmp_vec2_r;

const double tmp_4vec_lU[4] = {ADM_metric->lapseinv*smallb_zerol,
                              ADM_metric->lapseinv*(prims_l->BU[0]/prims_l->u0 + smallb_zerol*prims_l->vU[0]),
                              ADM_metric->lapseinv*(prims_l->BU[1]/prims_l->u0 + smallb_zerol*prims_l->vU[1]),
                              ADM_metric->lapseinv*(prims_l->BU[2]/prims_l->u0 + smallb_zerol*prims_l->vU[2])};
const double tmp_vec2_l = ghl_compute_vec2_from_vec4D( metric_aux.g4DD, tmp_4vec_lU);
const double tmp30 = prims_l->press + tmp_vec2_l/2.0;//+ 0.5*tmp22*tmp23 + 0.5*tmp22*tmp25 + 0.5*tmp22*tmp27 + 0.5*tmp29;
const double tmp31 = h_l*prims_l->rho + tmp_vec2_l;
const double soln_scal_r2    =  metric_aux.g4UU[0][0]*tmp_vec2_r/2.0;

const double tmp55 = ((prims_r->u0)*(prims_r->u0))*tmp54;
const double tmp62 = prims_r->u0*prims_r->vU[direction]*tmp54;
const double tmp32 = ((prims_l->u0)*(prims_l->u0))*tmp31;
const double tmp56 = g4UU00*tmp53 - ((tmp40)*(tmp40)) + tmp55;
const double tmp60 = prims_l->u0*prims_l->vU[direction]*tmp31;
const double tmp64 = g4UU01*tmp53 + prims_r->vU[0]*tmp55 - tmp40*tmp43;
const double tmp65 = g4UU02*tmp53 + prims_r->vU[1]*tmp55 - tmp40*tmp44;
const double tmp66 = g4UU03*tmp53 + prims_r->vU[2]*tmp55 - tmp40*tmp45;
const double tmp57 = g4UU00*tmp30 - ((tmp17)*(tmp17)) + tmp32;
const double tmp67 = g4UU01*tmp30 + prims_l->vU[0]*tmp32 - tmp17*tmp20;
const double tmp68 = g4UU02*tmp30 + prims_l->vU[1]*tmp32 - tmp17*tmp21;
const double tmp69 = g4UU03*tmp30 + prims_l->vU[2]*tmp32 - tmp17*tmp22;


cons_rhs->rho = tmp0*(cmax*prims_l->vU[direction]*tmp1 + cmin*prims_r->vU[direction]*tmp3 - tmp5*(-tmp1 + tmp3));
cons_rhs->tau += - tmp0*(tmp5*(tmp33*tmp56 - tmp33*tmp57)) - cons_rhs->rho;

//cons_rhs->tau = tmp0*(cmax*(tmp33*(metric_aux.g4UU[0][direction]*tmp30 + prims_l->vU[direction]*tmp32 - tmp17*tmp19))
//                    + cmin*(tmp33*(metric_aux.g4UU[0][direction]*tmp53 + prims_r->vU[direction]*tmp55 - tmp40*tmp42))
//                    - tmp5*(tmp33*tmp56 - tmp33*tmp57)) - cons_rhs->rho;
//cons_rhs->SD[0] = tmp0*(-tmp5*(sqrt_detg*(g4DD01*tmp56 + g4DD11*tmp64 + g4DD12*tmp65 + g4DD13*tmp66) - sqrt_detg*(g4DD01*tmp57 + g4DD11*tmp67 + g4DD12*tmp68 + g4DD13*tmp69)) + tmp61*(-tmp19*tmp23 + tmp30*tmp58 + tmp59*tmp60) + tmp63*(-tmp42*tmp46 + tmp53*tmp58 + tmp59*tmp62));
//cons_rhs->SD[1] = tmp0*(-tmp5*(sqrt_detg*(g4DD02*tmp56 + g4DD12*tmp64 + g4DD22*tmp65 + g4DD23*tmp66) - sqrt_detg*(g4DD02*tmp57 + g4DD12*tmp67 + g4DD22*tmp68 + g4DD23*tmp69)) + tmp61*(tmp10*tmp60 - tmp19*tmp25 + tmp30*tmp70) + tmp63*(tmp10*tmp62 - tmp42*tmp48 + tmp53*tmp70));
//cons_rhs->SD[2] = tmp0*(-tmp5*(sqrt_detg*(g4DD03*tmp56 + g4DD13*tmp64 + g4DD23*tmp65 + g4DD33*tmp66) - sqrt_detg*(g4DD03*tmp57 + g4DD13*tmp67 + g4DD23*tmp68 + g4DD33*tmp69)) + tmp61*(tmp13*tmp60 - tmp19*tmp27 + tmp30*tmp71) + tmp63*(tmp13*tmp62 - tmp42*tmp50 + tmp53*tmp71));
}

void ghl_calculate_HLLE_fluxes_hybrid_entropy2(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

const double g4DD00 = metric_aux.g4DD[0][0];
const double g4DD01 = metric_aux.g4DD[0][1];
const double g4DD02 = metric_aux.g4DD[0][2];
const double g4DD03 = metric_aux.g4DD[0][3];
const double g4DD11 = metric_aux.g4DD[1][1];
const double g4DD12 = metric_aux.g4DD[1][2];
const double g4DD13 = metric_aux.g4DD[1][3];
const double g4DD22 = metric_aux.g4DD[2][2];
const double g4DD23 = metric_aux.g4DD[3][3];
const double g4DD33 = metric_aux.g4DD[3][3];

const double g4UU00 = metric_aux.g4UU[0][0];
const double g4UU01 = metric_aux.g4UU[0][1];
const double g4UU02 = metric_aux.g4UU[0][2];
const double g4UU03 = metric_aux.g4UU[0][3];
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;

const double tmp0 = (1.0/(cmax + cmin));
const double tmp7 = cmax*cmin;
const double tmp9 = prims_l->u0*prims_l->vU[0];
const double tmp10 = prims_l->u0*prims_l->vU[1];
const double tmp11 = prims_l->u0*prims_l->vU[2];
const double tmp20 = (1.0/(prims_l->u0));
const double tmp35 = ADM_metric->lapse*sqrt_detg;
const double tmp36 = prims_r->u0*prims_r->vU[0];
const double tmp37 = prims_r->u0*prims_r->vU[1];
const double tmp38 = prims_r->u0*prims_r->vU[2];
const double tmp43 = (1.0/(prims_r->u0));
const double tmp60 = g4DD01*metric_aux.g4UU[0][direction] + g4DD11*metric_aux.g4UU[1][direction] + g4DD12*metric_aux.g4UU[2][direction] + g4DD13*metric_aux.g4UU[3][direction];
const double tmp63 = cmax*sqrt_detg;
const double tmp65 = cmin*sqrt_detg;
const double tmp72 = g4DD02*metric_aux.g4UU[0][direction] + g4DD12*metric_aux.g4UU[1][direction] + g4DD22*metric_aux.g4UU[2][direction] + g4DD23*metric_aux.g4UU[3][direction];
const double tmp73 = g4DD03*metric_aux.g4UU[0][direction] + g4DD13*metric_aux.g4UU[1][direction] + g4DD23*metric_aux.g4UU[2][direction] + g4DD33*metric_aux.g4UU[3][direction];
const double tmp2 = prims_l->rho*prims_l->u0*sqrt_detg;
const double tmp5 = prims_r->rho*prims_r->u0*sqrt_detg;
const double tmp12 = g4DD01*prims_l->u0 + g4DD11*tmp9 + g4DD12*tmp10 + g4DD13*tmp11;
const double tmp15 = g4DD02*prims_l->u0 + g4DD12*tmp9 + g4DD22*tmp10 + g4DD23*tmp11;
const double tmp42 = ADM_metric->lapseinv*prims_r->BU[0]*(g4DD01*prims_r->u0 + g4DD11*tmp36 + g4DD12*tmp37 + g4DD13*tmp38) + ADM_metric->lapseinv*prims_r->BU[1]*(g4DD02*prims_r->u0 + g4DD12*tmp36 + g4DD22*tmp37 + g4DD23*tmp38) + ADM_metric->lapseinv*prims_r->BU[2]*(g4DD03*prims_r->u0 + g4DD13*tmp36 + g4DD23*tmp37 + g4DD33*tmp38);
const double tmp61 = g4DD00*prims_l->u0 + g4DD01*prims_l->u0*prims_l->vU[0] + g4DD02*prims_l->u0*prims_l->vU[1] + g4DD03*prims_l->u0*prims_l->vU[2];
const double tmp74 = prims_l->entropy*prims_l->u0*sqrt_detg;
const double tmp75 = prims_r->entropy*prims_r->u0*sqrt_detg;
const double tmp19 = ADM_metric->lapseinv*prims_l->BU[0]*tmp12 + ADM_metric->lapseinv*prims_l->BU[1]*tmp15 + ADM_metric->lapseinv*prims_l->BU[2]*(g4DD03*prims_l->u0 + g4DD13*tmp9 + g4DD23*tmp10 + g4DD33*tmp11);
const double tmp44 = ADM_metric->lapseinv*prims_r->BU[direction]*tmp43 + prims_r->vU[direction]*tmp42;
const double tmp45 = tmp43*(ADM_metric->lapseinv*prims_r->BU[0] + tmp36*tmp42);
const double tmp46 = tmp43*(ADM_metric->lapseinv*prims_r->BU[1] + tmp37*tmp42);
const double tmp47 = tmp43*(ADM_metric->lapseinv*prims_r->BU[2] + tmp38*tmp42);
const double tmp21 = ADM_metric->lapseinv*prims_l->BU[direction]*tmp20 + prims_l->vU[direction]*tmp19;
const double tmp22 = tmp20*(ADM_metric->lapseinv*prims_l->BU[0] + tmp19*tmp9);
const double tmp23 = tmp20*(ADM_metric->lapseinv*prims_l->BU[1] + tmp10*tmp19);
const double tmp24 = tmp20*(ADM_metric->lapseinv*prims_l->BU[2] + tmp11*tmp19);
const double tmp48 = g4DD00*tmp42 + g4DD01*tmp45 + g4DD02*tmp46 + g4DD03*tmp47;
const double tmp50 = g4DD01*tmp42 + g4DD11*tmp45 + g4DD12*tmp46 + g4DD13*tmp47;
const double tmp52 = g4DD02*tmp42 + g4DD12*tmp45 + g4DD22*tmp46 + g4DD23*tmp47;
const double tmp54 = tmp47*(g4DD03*tmp42 + g4DD13*tmp45 + g4DD23*tmp46 + g4DD33*tmp47);
const double tmp25 = g4DD00*tmp19 + g4DD01*tmp22 + g4DD02*tmp23 + g4DD03*tmp24;
const double tmp27 = g4DD01*tmp19 + g4DD11*tmp22 + g4DD12*tmp23 + g4DD13*tmp24;
const double tmp29 = g4DD02*tmp19 + g4DD12*tmp22 + g4DD22*tmp23 + g4DD23*tmp24;
const double tmp31 = tmp24*(g4DD03*tmp19 + g4DD13*tmp22 + g4DD23*tmp23 + g4DD33*tmp24);
const double tmp55 = prims_r->press + 0.5*tmp47*tmp48 + 0.5*tmp47*tmp50 + 0.5*tmp47*tmp52 + 0.5*tmp54;
const double tmp56 = h_r*prims_r->rho + tmp47*tmp48 + tmp47*tmp50 + tmp47*tmp52 + tmp54;
const double tmp32 = prims_l->press + 0.5*tmp24*tmp25 + 0.5*tmp24*tmp27 + 0.5*tmp24*tmp29 + 0.5*tmp31;
const double tmp33 = h_l*prims_l->rho + tmp24*tmp25 + tmp24*tmp27 + tmp24*tmp29 + tmp31;
const double tmp57 = ((prims_r->u0)*(prims_r->u0))*tmp56;
const double tmp64 = prims_r->u0*prims_r->vU[direction]*tmp56;
const double tmp34 = ((prims_l->u0)*(prims_l->u0))*tmp33;
const double tmp58 = g4UU00*tmp55 - ((tmp42)*(tmp42)) + tmp57;
const double tmp62 = prims_l->u0*prims_l->vU[direction]*tmp33;
const double tmp66 = g4UU01*tmp55 + prims_r->vU[0]*tmp57 - tmp42*tmp45;
const double tmp67 = g4UU02*tmp55 + prims_r->vU[1]*tmp57 - tmp42*tmp46;
const double tmp68 = g4UU03*tmp55 + prims_r->vU[2]*tmp57 - tmp42*tmp47;
const double tmp59 = g4UU00*tmp32 - ((tmp19)*(tmp19)) + tmp34;
const double tmp69 = g4UU01*tmp32 + prims_l->vU[0]*tmp34 - tmp19*tmp22;
const double tmp70 = g4UU02*tmp32 + prims_l->vU[1]*tmp34 - tmp19*tmp23;
const double tmp71 = g4UU03*tmp32 + prims_l->vU[2]*tmp34 - tmp19*tmp24;
cons_rhs->rho = tmp0*(cmax*prims_l->vU[direction]*tmp2 + cmin*prims_r->vU[direction]*tmp5 - tmp7*(-tmp2 + tmp5));
cons_rhs->tau = tmp0*(cmax*(-prims_l->vU[direction]*tmp2 + tmp35*(metric_aux.g4UU[0][direction]*tmp32 + prims_l->vU[direction]*tmp34 - tmp19*tmp21)) + cmin*(-prims_r->vU[direction]*tmp5 + tmp35*(metric_aux.g4UU[0][direction]*tmp55 + prims_r->vU[direction]*tmp57 - tmp42*tmp44)) - tmp7*(tmp2 + tmp35*tmp58 - tmp35*tmp59 - tmp5));
cons_rhs->SD[0] = tmp0*(tmp63*(-tmp21*tmp25 + tmp32*tmp60 + tmp61*tmp62) + tmp65*(-tmp44*tmp48 + tmp55*tmp60 + tmp61*tmp64) - tmp7*(sqrt_detg*(g4DD01*tmp58 + g4DD11*tmp66 + g4DD12*tmp67 + g4DD13*tmp68) - sqrt_detg*(g4DD01*tmp59 + g4DD11*tmp69 + g4DD12*tmp70 + g4DD13*tmp71)));
cons_rhs->SD[1] = tmp0*(tmp63*(tmp12*tmp62 - tmp21*tmp27 + tmp32*tmp72) + tmp65*(tmp12*tmp64 - tmp44*tmp50 + tmp55*tmp72) - tmp7*(sqrt_detg*(g4DD02*tmp58 + g4DD12*tmp66 + g4DD22*tmp67 + g4DD23*tmp68) - sqrt_detg*(g4DD02*tmp59 + g4DD12*tmp69 + g4DD22*tmp70 + g4DD23*tmp71)));
cons_rhs->SD[2] = tmp0*(tmp63*(tmp15*tmp62 - tmp21*tmp29 + tmp32*tmp73) + tmp65*(tmp15*tmp64 - tmp44*tmp52 + tmp55*tmp73) - tmp7*(sqrt_detg*(g4DD03*tmp58 + g4DD13*tmp66 + g4DD23*tmp67 + g4DD33*tmp68) - sqrt_detg*(g4DD03*tmp59 + g4DD13*tmp69 + g4DD23*tmp70 + g4DD33*tmp71)));
cons_rhs->entropy = tmp0*(cmax*prims_l->vU[direction]*tmp74 + cmin*prims_r->vU[direction]*tmp75 - tmp7*(-tmp74 + tmp75));
}

void ghl_calculate_HLLE_fluxes_tabulated(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

const double g4DD00 = metric_aux.g4DD[0][0];
const double g4DD01 = metric_aux.g4DD[0][1];
const double g4DD02 = metric_aux.g4DD[0][2];
const double g4DD03 = metric_aux.g4DD[0][3];
const double g4DD11 = metric_aux.g4DD[1][1];
const double g4DD12 = metric_aux.g4DD[1][2];
const double g4DD13 = metric_aux.g4DD[1][3];
const double g4DD22 = metric_aux.g4DD[2][2];
const double g4DD23 = metric_aux.g4DD[3][3];
const double g4DD33 = metric_aux.g4DD[3][3];

const double g4UU00 = metric_aux.g4UU[0][0];
const double g4UU01 = metric_aux.g4UU[0][1];
const double g4UU02 = metric_aux.g4UU[0][2];
const double g4UU03 = metric_aux.g4UU[0][3];
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;

const double tmp0 = (1.0/(cmax + cmin));
const double tmp1 = prims_l->rho*prims_l->u0*sqrt_detg;
const double tmp4 = prims_r->rho*prims_r->u0*sqrt_detg;
const double tmp7 = cmax*cmin;
const double tmp9 = prims_l->u0*prims_l->vU[0];
const double tmp10 = prims_l->u0*prims_l->vU[1];
const double tmp11 = prims_l->u0*prims_l->vU[2];
const double tmp20 = (1.0/(prims_l->u0));
const double tmp35 = ADM_metric->lapse*sqrt_detg;
const double tmp36 = prims_r->u0*prims_r->vU[0];
const double tmp37 = prims_r->u0*prims_r->vU[1];
const double tmp38 = prims_r->u0*prims_r->vU[2];
const double tmp43 = (1.0/(prims_r->u0));
const double tmp60 = g4DD01*metric_aux.g4UU[0][direction] + g4DD11*metric_aux.g4UU[1][direction] + g4DD12*metric_aux.g4UU[2][direction] + g4DD13*metric_aux.g4UU[3][direction];
const double tmp63 = cmax*sqrt_detg;
const double tmp65 = cmin*sqrt_detg;
const double tmp72 = g4DD02*metric_aux.g4UU[0][direction] + g4DD12*metric_aux.g4UU[1][direction] + g4DD22*metric_aux.g4UU[2][direction] + g4DD23*metric_aux.g4UU[3][direction];
const double tmp73 = g4DD03*metric_aux.g4UU[0][direction] + g4DD13*metric_aux.g4UU[1][direction] + g4DD23*metric_aux.g4UU[2][direction] + g4DD33*metric_aux.g4UU[3][direction];
const double tmp12 = g4DD01*prims_l->u0 + g4DD11*tmp9 + g4DD12*tmp10 + g4DD13*tmp11;
const double tmp15 = g4DD02*prims_l->u0 + g4DD12*tmp9 + g4DD22*tmp10 + g4DD23*tmp11;
const double tmp42 = ADM_metric->lapseinv*prims_r->BU[0]*(g4DD01*prims_r->u0 + g4DD11*tmp36 + g4DD12*tmp37 + g4DD13*tmp38) + ADM_metric->lapseinv*prims_r->BU[1]*(g4DD02*prims_r->u0 + g4DD12*tmp36 + g4DD22*tmp37 + g4DD23*tmp38) + ADM_metric->lapseinv*prims_r->BU[2]*(g4DD03*prims_r->u0 + g4DD13*tmp36 + g4DD23*tmp37 + g4DD33*tmp38);
const double tmp61 = g4DD00*prims_l->u0 + g4DD01*prims_l->u0*prims_l->vU[0] + g4DD02*prims_l->u0*prims_l->vU[1] + g4DD03*prims_l->u0*prims_l->vU[2];
const double tmp3 = cmax*prims_l->vU[direction]*tmp1;
const double tmp6 = cmin*prims_r->vU[direction]*tmp4;
const double tmp19 = ADM_metric->lapseinv*prims_l->BU[0]*tmp12 + ADM_metric->lapseinv*prims_l->BU[1]*tmp15 + ADM_metric->lapseinv*prims_l->BU[2]*(g4DD03*prims_l->u0 + g4DD13*tmp9 + g4DD23*tmp10 + g4DD33*tmp11);
const double tmp44 = ADM_metric->lapseinv*prims_r->BU[direction]*tmp43 + prims_r->vU[direction]*tmp42;
const double tmp45 = tmp43*(ADM_metric->lapseinv*prims_r->BU[0] + tmp36*tmp42);
const double tmp46 = tmp43*(ADM_metric->lapseinv*prims_r->BU[1] + tmp37*tmp42);
const double tmp47 = tmp43*(ADM_metric->lapseinv*prims_r->BU[2] + tmp38*tmp42);
const double tmp21 = ADM_metric->lapseinv*prims_l->BU[direction]*tmp20 + prims_l->vU[direction]*tmp19;
const double tmp22 = tmp20*(ADM_metric->lapseinv*prims_l->BU[0] + tmp19*tmp9);
const double tmp23 = tmp20*(ADM_metric->lapseinv*prims_l->BU[1] + tmp10*tmp19);
const double tmp24 = tmp20*(ADM_metric->lapseinv*prims_l->BU[2] + tmp11*tmp19);
const double tmp48 = g4DD00*tmp42 + g4DD01*tmp45 + g4DD02*tmp46 + g4DD03*tmp47;
const double tmp50 = g4DD01*tmp42 + g4DD11*tmp45 + g4DD12*tmp46 + g4DD13*tmp47;
const double tmp52 = g4DD02*tmp42 + g4DD12*tmp45 + g4DD22*tmp46 + g4DD23*tmp47;
const double tmp54 = tmp47*(g4DD03*tmp42 + g4DD13*tmp45 + g4DD23*tmp46 + g4DD33*tmp47);
const double tmp25 = g4DD00*tmp19 + g4DD01*tmp22 + g4DD02*tmp23 + g4DD03*tmp24;
const double tmp27 = g4DD01*tmp19 + g4DD11*tmp22 + g4DD12*tmp23 + g4DD13*tmp24;
const double tmp29 = g4DD02*tmp19 + g4DD12*tmp22 + g4DD22*tmp23 + g4DD23*tmp24;
const double tmp31 = tmp24*(g4DD03*tmp19 + g4DD13*tmp22 + g4DD23*tmp23 + g4DD33*tmp24);
const double tmp55 = prims_r->press + 0.5*tmp47*tmp48 + 0.5*tmp47*tmp50 + 0.5*tmp47*tmp52 + 0.5*tmp54;
const double tmp56 = h_r*prims_r->rho + tmp47*tmp48 + tmp47*tmp50 + tmp47*tmp52 + tmp54;
const double tmp32 = prims_l->press + 0.5*tmp24*tmp25 + 0.5*tmp24*tmp27 + 0.5*tmp24*tmp29 + 0.5*tmp31;
const double tmp33 = h_l*prims_l->rho + tmp24*tmp25 + tmp24*tmp27 + tmp24*tmp29 + tmp31;
const double tmp57 = ((prims_r->u0)*(prims_r->u0))*tmp56;
const double tmp64 = prims_r->u0*prims_r->vU[direction]*tmp56;
const double tmp34 = ((prims_l->u0)*(prims_l->u0))*tmp33;
const double tmp58 = g4UU00*tmp55 - ((tmp42)*(tmp42)) + tmp57;
const double tmp62 = prims_l->u0*prims_l->vU[direction]*tmp33;
const double tmp66 = g4UU01*tmp55 + prims_r->vU[0]*tmp57 - tmp42*tmp45;
const double tmp67 = g4UU02*tmp55 + prims_r->vU[1]*tmp57 - tmp42*tmp46;
const double tmp68 = g4UU03*tmp55 + prims_r->vU[2]*tmp57 - tmp42*tmp47;
const double tmp59 = g4UU00*tmp32 - ((tmp19)*(tmp19)) + tmp34;
const double tmp69 = g4UU01*tmp32 + prims_l->vU[0]*tmp34 - tmp19*tmp22;
const double tmp70 = g4UU02*tmp32 + prims_l->vU[1]*tmp34 - tmp19*tmp23;
const double tmp71 = g4UU03*tmp32 + prims_l->vU[2]*tmp34 - tmp19*tmp24;
cons_rhs->rho = tmp0*(tmp3 + tmp6 - tmp7*(-tmp1 + tmp4));
cons_rhs->tau = tmp0*(cmax*(-prims_l->vU[direction]*tmp1 + tmp35*(metric_aux.g4UU[0][direction]*tmp32 + prims_l->vU[direction]*tmp34 - tmp19*tmp21)) + cmin*(-prims_r->vU[direction]*tmp4 + tmp35*(metric_aux.g4UU[0][direction]*tmp55 + prims_r->vU[direction]*tmp57 - tmp42*tmp44)) - tmp7*(tmp1 + tmp35*tmp58 - tmp35*tmp59 - tmp4));
cons_rhs->SD[0] = tmp0*(tmp63*(-tmp21*tmp25 + tmp32*tmp60 + tmp61*tmp62) + tmp65*(-tmp44*tmp48 + tmp55*tmp60 + tmp61*tmp64) - tmp7*(sqrt_detg*(g4DD01*tmp58 + g4DD11*tmp66 + g4DD12*tmp67 + g4DD13*tmp68) - sqrt_detg*(g4DD01*tmp59 + g4DD11*tmp69 + g4DD12*tmp70 + g4DD13*tmp71)));
cons_rhs->SD[1] = tmp0*(tmp63*(tmp12*tmp62 - tmp21*tmp27 + tmp32*tmp72) + tmp65*(tmp12*tmp64 - tmp44*tmp50 + tmp55*tmp72) - tmp7*(sqrt_detg*(g4DD02*tmp58 + g4DD12*tmp66 + g4DD22*tmp67 + g4DD23*tmp68) - sqrt_detg*(g4DD02*tmp59 + g4DD12*tmp69 + g4DD22*tmp70 + g4DD23*tmp71)));
cons_rhs->SD[2] = tmp0*(tmp63*(tmp15*tmp62 - tmp21*tmp29 + tmp32*tmp73) + tmp65*(tmp15*tmp64 - tmp44*tmp52 + tmp55*tmp73) - tmp7*(sqrt_detg*(g4DD03*tmp58 + g4DD13*tmp66 + g4DD23*tmp67 + g4DD33*tmp68) - sqrt_detg*(g4DD03*tmp59 + g4DD13*tmp69 + g4DD23*tmp70 + g4DD33*tmp71)));
cons_rhs->Y_e = tmp0*(prims_l->Y_e*tmp3 + prims_r->Y_e*tmp6 - tmp7*(-prims_l->Y_e*tmp1 + prims_r->Y_e*tmp4));
}

void ghl_calculate_HLLE_fluxes_tabulated_entropy(
          const int direction,
          const ghl_eos_parameters *restrict eos,
          const ghl_metric_quantities *restrict ADM_metric,
          ghl_primitive_quantities *restrict prims_r,
          ghl_primitive_quantities *restrict prims_l,
          const double cmin,
          const double cmax,
          ghl_conservative_quantities *restrict cons_rhs) {

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

const double g4DD00 = metric_aux.g4DD[0][0];
const double g4DD01 = metric_aux.g4DD[0][1];
const double g4DD02 = metric_aux.g4DD[0][2];
const double g4DD03 = metric_aux.g4DD[0][3];
const double g4DD11 = metric_aux.g4DD[1][1];
const double g4DD12 = metric_aux.g4DD[1][2];
const double g4DD13 = metric_aux.g4DD[1][3];
const double g4DD22 = metric_aux.g4DD[2][2];
const double g4DD23 = metric_aux.g4DD[3][3];
const double g4DD33 = metric_aux.g4DD[3][3];

const double g4UU00 = metric_aux.g4UU[0][0];
const double g4UU01 = metric_aux.g4UU[0][1];
const double g4UU02 = metric_aux.g4UU[0][2];
const double g4UU03 = metric_aux.g4UU[0][3];
  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;

const double tmp0 = (1.0/(cmax + cmin));
const double tmp9 = cmax*cmin;
const double tmp11 = prims_l->u0*prims_l->vU[0];
const double tmp12 = prims_l->u0*prims_l->vU[1];
const double tmp13 = prims_l->u0*prims_l->vU[2];
const double tmp22 = (1.0/(prims_l->u0));
const double tmp37 = ADM_metric->lapse*sqrt_detg;
const double tmp38 = prims_r->u0*prims_r->vU[0];
const double tmp39 = prims_r->u0*prims_r->vU[1];
const double tmp40 = prims_r->u0*prims_r->vU[2];
const double tmp45 = (1.0/(prims_r->u0));
const double tmp62 = g4DD01*metric_aux.g4UU[0][direction] + g4DD11*metric_aux.g4UU[1][direction] + g4DD12*metric_aux.g4UU[2][direction] + g4DD13*metric_aux.g4UU[3][direction];
const double tmp65 = cmax*sqrt_detg;
const double tmp67 = cmin*sqrt_detg;
const double tmp74 = g4DD02*metric_aux.g4UU[0][direction] + g4DD12*metric_aux.g4UU[1][direction] + g4DD22*metric_aux.g4UU[2][direction] + g4DD23*metric_aux.g4UU[3][direction];
const double tmp75 = g4DD03*metric_aux.g4UU[0][direction] + g4DD13*metric_aux.g4UU[1][direction] + g4DD23*metric_aux.g4UU[2][direction] + g4DD33*metric_aux.g4UU[3][direction];
const double tmp2 = prims_l->rho*prims_l->u0*sqrt_detg;
const double tmp6 = prims_r->rho*prims_r->u0*sqrt_detg;
const double tmp14 = g4DD01*prims_l->u0 + g4DD11*tmp11 + g4DD12*tmp12 + g4DD13*tmp13;
const double tmp17 = g4DD02*prims_l->u0 + g4DD12*tmp11 + g4DD22*tmp12 + g4DD23*tmp13;
const double tmp44 = ADM_metric->lapseinv*prims_r->BU[0]*(g4DD01*prims_r->u0 + g4DD11*tmp38 + g4DD12*tmp39 + g4DD13*tmp40) + ADM_metric->lapseinv*prims_r->BU[1]*(g4DD02*prims_r->u0 + g4DD12*tmp38 + g4DD22*tmp39 + g4DD23*tmp40) + ADM_metric->lapseinv*prims_r->BU[2]*(g4DD03*prims_r->u0 + g4DD13*tmp38 + g4DD23*tmp39 + g4DD33*tmp40);
const double tmp63 = g4DD00*prims_l->u0 + g4DD01*prims_l->u0*prims_l->vU[0] + g4DD02*prims_l->u0*prims_l->vU[1] + g4DD03*prims_l->u0*prims_l->vU[2];
const double tmp76 = prims_l->entropy*prims_l->u0*sqrt_detg;
const double tmp77 = prims_r->entropy*prims_r->u0*sqrt_detg;
const double tmp21 = ADM_metric->lapseinv*prims_l->BU[0]*tmp14 + ADM_metric->lapseinv*prims_l->BU[1]*tmp17 + ADM_metric->lapseinv*prims_l->BU[2]*(g4DD03*prims_l->u0 + g4DD13*tmp11 + g4DD23*tmp12 + g4DD33*tmp13);
const double tmp46 = ADM_metric->lapseinv*prims_r->BU[direction]*tmp45 + prims_r->vU[direction]*tmp44;
const double tmp47 = tmp45*(ADM_metric->lapseinv*prims_r->BU[0] + tmp38*tmp44);
const double tmp48 = tmp45*(ADM_metric->lapseinv*prims_r->BU[1] + tmp39*tmp44);
const double tmp49 = tmp45*(ADM_metric->lapseinv*prims_r->BU[2] + tmp40*tmp44);
const double tmp4 = cmax*prims_l->vU[direction]*tmp2;
const double tmp8 = cmin*prims_r->vU[direction]*tmp6;
const double tmp23 = ADM_metric->lapseinv*prims_l->BU[direction]*tmp22 + prims_l->vU[direction]*tmp21;
const double tmp24 = tmp22*(ADM_metric->lapseinv*prims_l->BU[0] + tmp11*tmp21);
const double tmp25 = tmp22*(ADM_metric->lapseinv*prims_l->BU[1] + tmp12*tmp21);
const double tmp26 = tmp22*(ADM_metric->lapseinv*prims_l->BU[2] + tmp13*tmp21);
const double tmp50 = g4DD00*tmp44 + g4DD01*tmp47 + g4DD02*tmp48 + g4DD03*tmp49;
const double tmp52 = g4DD01*tmp44 + g4DD11*tmp47 + g4DD12*tmp48 + g4DD13*tmp49;
const double tmp54 = g4DD02*tmp44 + g4DD12*tmp47 + g4DD22*tmp48 + g4DD23*tmp49;
const double tmp56 = tmp49*(g4DD03*tmp44 + g4DD13*tmp47 + g4DD23*tmp48 + g4DD33*tmp49);
const double tmp27 = g4DD00*tmp21 + g4DD01*tmp24 + g4DD02*tmp25 + g4DD03*tmp26;
const double tmp29 = g4DD01*tmp21 + g4DD11*tmp24 + g4DD12*tmp25 + g4DD13*tmp26;
const double tmp31 = g4DD02*tmp21 + g4DD12*tmp24 + g4DD22*tmp25 + g4DD23*tmp26;
const double tmp33 = tmp26*(g4DD03*tmp21 + g4DD13*tmp24 + g4DD23*tmp25 + g4DD33*tmp26);
const double tmp57 = prims_r->press + 0.5*tmp49*tmp50 + 0.5*tmp49*tmp52 + 0.5*tmp49*tmp54 + 0.5*tmp56;
const double tmp58 = h_r*prims_r->rho + tmp49*tmp50 + tmp49*tmp52 + tmp49*tmp54 + tmp56;
const double tmp34 = prims_l->press + 0.5*tmp26*tmp27 + 0.5*tmp26*tmp29 + 0.5*tmp26*tmp31 + 0.5*tmp33;
const double tmp35 = h_l*prims_l->rho + tmp26*tmp27 + tmp26*tmp29 + tmp26*tmp31 + tmp33;
const double tmp59 = ((prims_r->u0)*(prims_r->u0))*tmp58;
const double tmp66 = prims_r->u0*prims_r->vU[direction]*tmp58;
const double tmp36 = ((prims_l->u0)*(prims_l->u0))*tmp35;
const double tmp60 = g4UU00*tmp57 - ((tmp44)*(tmp44)) + tmp59;
const double tmp64 = prims_l->u0*prims_l->vU[direction]*tmp35;
const double tmp68 = g4UU01*tmp57 + prims_r->vU[0]*tmp59 - tmp44*tmp47;
const double tmp69 = g4UU02*tmp57 + prims_r->vU[1]*tmp59 - tmp44*tmp48;
const double tmp70 = g4UU03*tmp57 + prims_r->vU[2]*tmp59 - tmp44*tmp49;
const double tmp61 = g4UU00*tmp34 - ((tmp21)*(tmp21)) + tmp36;
const double tmp71 = g4UU01*tmp34 + prims_l->vU[0]*tmp36 - tmp21*tmp24;
const double tmp72 = g4UU02*tmp34 + prims_l->vU[1]*tmp36 - tmp21*tmp25;
const double tmp73 = g4UU03*tmp34 + prims_l->vU[2]*tmp36 - tmp21*tmp26;
cons_rhs->rho = tmp0*(tmp4 + tmp8 - tmp9*(-tmp2 + tmp6));
cons_rhs->tau = tmp0*(cmax*(-prims_l->vU[direction]*tmp2 + tmp37*(metric_aux.g4UU[0][direction]*tmp34 + prims_l->vU[direction]*tmp36 - tmp21*tmp23)) + cmin*(-prims_r->vU[direction]*tmp6 + tmp37*(metric_aux.g4UU[0][direction]*tmp57 + prims_r->vU[direction]*tmp59 - tmp44*tmp46)) - tmp9*(tmp2 + tmp37*tmp60 - tmp37*tmp61 - tmp6));
cons_rhs->SD[0] = tmp0*(tmp65*(-tmp23*tmp27 + tmp34*tmp62 + tmp63*tmp64) + tmp67*(-tmp46*tmp50 + tmp57*tmp62 + tmp63*tmp66) - tmp9*(sqrt_detg*(g4DD01*tmp60 + g4DD11*tmp68 + g4DD12*tmp69 + g4DD13*tmp70) - sqrt_detg*(g4DD01*tmp61 + g4DD11*tmp71 + g4DD12*tmp72 + g4DD13*tmp73)));
cons_rhs->SD[1] = tmp0*(tmp65*(tmp14*tmp64 - tmp23*tmp29 + tmp34*tmp74) + tmp67*(tmp14*tmp66 - tmp46*tmp52 + tmp57*tmp74) - tmp9*(sqrt_detg*(g4DD02*tmp60 + g4DD12*tmp68 + g4DD22*tmp69 + g4DD23*tmp70) - sqrt_detg*(g4DD02*tmp61 + g4DD12*tmp71 + g4DD22*tmp72 + g4DD23*tmp73)));
cons_rhs->SD[2] = tmp0*(tmp65*(tmp17*tmp64 - tmp23*tmp31 + tmp34*tmp75) + tmp67*(tmp17*tmp66 - tmp46*tmp54 + tmp57*tmp75) - tmp9*(sqrt_detg*(g4DD03*tmp60 + g4DD13*tmp68 + g4DD23*tmp69 + g4DD33*tmp70) - sqrt_detg*(g4DD03*tmp61 + g4DD13*tmp71 + g4DD23*tmp72 + g4DD33*tmp73)));
cons_rhs->Y_e = tmp0*(prims_l->Y_e*tmp4 + prims_r->Y_e*tmp8 - tmp9*(-prims_l->Y_e*tmp2 + prims_r->Y_e*tmp6));
cons_rhs->entropy = tmp0*(cmax*prims_l->vU[direction]*tmp76 + cmin*prims_r->vU[direction]*tmp77 - tmp9*(-tmp76 + tmp77));
}
