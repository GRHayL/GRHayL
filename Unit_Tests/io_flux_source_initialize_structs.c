#include "unit_tests.h"
#include "flux_source_unit_test.h"

/*
 * Initialize structs used for flux and source calculations
 */
void initialize_structs(const int idx, 
    const double *restrict auxevol_gfs, 
    primitive_quantities *restrict prims, 
    primitive_quantities *restrict prims_r, 
    primitive_quantities *restrict prims_l, 
    metric_quantities *restrict metric, 
    metric_quantities *restrict metric_face,
    extrinsic_curvature *restrict curv, 
    metric_derivatives *restrict metric_derivs) {

const double U4U0  =  auxevol_gfs[IDX4ptS(U4U0GF, idx)];
const double U4RU0 =  auxevol_gfs[IDX4ptS(U4RU0GF, idx)];
const double U4LU0 =  auxevol_gfs[IDX4ptS(U4LU0GF, idx)];
    
prims->u0 = U4U0;
prims_r->u0 = U4RU0;
prims_l->u0 = U4LU0;
    
prims->vx = auxevol_gfs[IDX4ptS(U4U1GF, idx)] / U4U0;
prims_r->vx = auxevol_gfs[IDX4ptS(U4RU1GF, idx)] / U4RU0;
prims_l->vx = auxevol_gfs[IDX4ptS(U4LU1GF, idx)] / U4LU0;
    
prims->vy = auxevol_gfs[IDX4ptS(U4U2GF, idx)] / U4U0;
prims_r->vy = auxevol_gfs[IDX4ptS(U4RU2GF, idx)] / U4RU0;
prims_l->vy = auxevol_gfs[IDX4ptS(U4LU2GF, idx)] / U4LU0;
    
prims->vz = auxevol_gfs[IDX4ptS(U4U3GF, idx)] / U4U0;
prims_r->vz = auxevol_gfs[IDX4ptS(U4RU3GF, idx)] / U4RU0;
prims_l->vz = auxevol_gfs[IDX4ptS(U4LU3GF, idx)] / U4LU0;
    
prims->Bx = auxevol_gfs[IDX4ptS(BU0GF, idx)];
prims_r->Bx = auxevol_gfs[IDX4ptS(BRU0GF, idx)];
prims_l->Bx = auxevol_gfs[IDX4ptS(BLU0GF, idx)];
    
prims->By = auxevol_gfs[IDX4ptS(BU1GF, idx)];
prims_r->By = auxevol_gfs[IDX4ptS(BRU1GF, idx)];
prims_l->By = auxevol_gfs[IDX4ptS(BLU1GF, idx)];
    
prims->Bz = auxevol_gfs[IDX4ptS(BU2GF, idx)];
prims_r->Bz = auxevol_gfs[IDX4ptS(BRU2GF, idx)];
prims_l->Bz = auxevol_gfs[IDX4ptS(BLU2GF, idx)];
    
prims->press = auxevol_gfs[IDX4ptS(PGF, idx)];
prims_r->press = auxevol_gfs[IDX4ptS(P_RGF, idx)];
prims_l->press = auxevol_gfs[IDX4ptS(P_LGF, idx)];
    
prims->rho = auxevol_gfs[IDX4ptS(RHOBGF, idx)];
prims_r->rho = auxevol_gfs[IDX4ptS(RHOB_RGF, idx)];
prims_l->rho = auxevol_gfs[IDX4ptS(RHOB_LGF, idx)];
    
metric->lapse = auxevol_gfs[IDX4ptS(ALPHAGF, idx)];
 
metric_face->lapse = auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx)];
    
metric->betax = auxevol_gfs[IDX4ptS(BETAU0GF, idx)];
 
metric_face->betax = auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx)];
    
metric->betay = auxevol_gfs[IDX4ptS(BETAU1GF, idx)];
 
metric_face->betay = auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx)];
    
metric->betaz = auxevol_gfs[IDX4ptS(BETAU2GF, idx)];
 
metric_face->betaz = auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx)];
    
metric->adm_gxy = auxevol_gfs[IDX4ptS(GAMMADD01GF, idx)];
 
metric_face->adm_gxy = auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx)];
    
metric->adm_gxz = auxevol_gfs[IDX4ptS(GAMMADD02GF, idx)];
 
metric_face->adm_gxz = auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx)];
    
metric->adm_gyz = auxevol_gfs[IDX4ptS(GAMMADD12GF, idx)];
 
metric_face->adm_gyz = auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx)];
    
metric->adm_gxx = auxevol_gfs[IDX4ptS(GAMMADD00GF, idx)];
 
metric_face->adm_gxx = auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx)];
    
metric->adm_gyy = auxevol_gfs[IDX4ptS(GAMMADD11GF, idx)];
 
metric_face->adm_gyy = auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx)];
    
metric->adm_gzz = auxevol_gfs[IDX4ptS(GAMMADD22GF, idx)];
 
metric_face->adm_gzz = auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx)];
    
curv->Kxy = auxevol_gfs[IDX4ptS(KDD01GF, idx)];
    
curv->Kxz = auxevol_gfs[IDX4ptS(KDD02GF, idx)];
    
curv->Kyz = auxevol_gfs[IDX4ptS(KDD12GF, idx)];
    
curv->Kxx = auxevol_gfs[IDX4ptS(KDD00GF, idx)];
    
curv->Kyy = auxevol_gfs[IDX4ptS(KDD11GF, idx)];
    
curv->Kzz = auxevol_gfs[IDX4ptS(KDD22GF, idx)];

metric_derivs->lapse[0] = auxevol_gfs[IDX4ptS(ALPHA_DD0GF, idx)];
metric_derivs->lapse[1] = auxevol_gfs[IDX4ptS(ALPHA_DD1GF, idx)];
metric_derivs->lapse[2] = auxevol_gfs[IDX4ptS(ALPHA_DD2GF, idx)];
metric_derivs->betax[0] = auxevol_gfs[IDX4ptS(BETAU_DD00GF, idx)];
metric_derivs->betay[0] = auxevol_gfs[IDX4ptS(BETAU_DD10GF, idx)];
metric_derivs->betaz[0] = auxevol_gfs[IDX4ptS(BETAU_DD20GF, idx)];
metric_derivs->adm_gxx[0] = auxevol_gfs[IDX4ptS(GAMMADD_DD000GF, idx)];
metric_derivs->adm_gxy[0] = auxevol_gfs[IDX4ptS(GAMMADD_DD010GF, idx)];
metric_derivs->adm_gxz[0] = auxevol_gfs[IDX4ptS(GAMMADD_DD020GF, idx)];
metric_derivs->adm_gyy[0] = auxevol_gfs[IDX4ptS(GAMMADD_DD110GF, idx)];
metric_derivs->adm_gyz[0] = auxevol_gfs[IDX4ptS(GAMMADD_DD120GF, idx)];
metric_derivs->adm_gzz[0] = auxevol_gfs[IDX4ptS(GAMMADD_DD220GF, idx)];
metric_derivs->betax[1] = auxevol_gfs[IDX4ptS(BETAU_DD01GF, idx)];
metric_derivs->betay[1] = auxevol_gfs[IDX4ptS(BETAU_DD11GF, idx)];
metric_derivs->betaz[1] = auxevol_gfs[IDX4ptS(BETAU_DD21GF, idx)];
metric_derivs->adm_gxx[1] = auxevol_gfs[IDX4ptS(GAMMADD_DD001GF, idx)];
metric_derivs->adm_gxy[1] = auxevol_gfs[IDX4ptS(GAMMADD_DD011GF, idx)];
metric_derivs->adm_gxz[1] = auxevol_gfs[IDX4ptS(GAMMADD_DD021GF, idx)];
metric_derivs->adm_gyy[1] = auxevol_gfs[IDX4ptS(GAMMADD_DD111GF, idx)];
metric_derivs->adm_gyz[1] = auxevol_gfs[IDX4ptS(GAMMADD_DD121GF, idx)];
metric_derivs->adm_gzz[1] = auxevol_gfs[IDX4ptS(GAMMADD_DD221GF, idx)];
metric_derivs->betax[2] = auxevol_gfs[IDX4ptS(BETAU_DD02GF, idx)];
metric_derivs->betay[2] = auxevol_gfs[IDX4ptS(BETAU_DD12GF, idx)];
metric_derivs->betaz[2] = auxevol_gfs[IDX4ptS(BETAU_DD22GF, idx)];
metric_derivs->adm_gxx[2] = auxevol_gfs[IDX4ptS(GAMMADD_DD002GF, idx)];
metric_derivs->adm_gxy[2] = auxevol_gfs[IDX4ptS(GAMMADD_DD012GF, idx)];
metric_derivs->adm_gxz[2] = auxevol_gfs[IDX4ptS(GAMMADD_DD022GF, idx)];
metric_derivs->adm_gyy[2] = auxevol_gfs[IDX4ptS(GAMMADD_DD112GF, idx)];
metric_derivs->adm_gyz[2] = auxevol_gfs[IDX4ptS(GAMMADD_DD122GF, idx)];
metric_derivs->adm_gzz[2] = auxevol_gfs[IDX4ptS(GAMMADD_DD222GF, idx)];
}
