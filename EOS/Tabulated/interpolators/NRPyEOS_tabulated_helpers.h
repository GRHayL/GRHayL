
/*
 * (c) 2022 Leo Werneck
 *
 * This file contains modified functions from the original
 * helpers.hh file from the Zelmani eosdrivercxx repository.
 * Source: https://bitbucket.org/zelmani/eosdrivercxx
 */

//------------------------------------------
static inline __attribute__((always_inline))
int NRPyEOS_checkbounds(const NRPyEOS_params *restrict eos_params,
                        const double xrho,
                        const double xtemp,
                        const double xye) {

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  if(xrho > eos_params->eos_rhomax) {
    return 105;
  }
  if(xrho < eos_params->eos_rhomin) {
    return 106;
  }
  if(xye > eos_params->eos_yemax) {
    return 101;
  }
  if(xye < eos_params->eos_yemin) {
    // this is probably not pure and should be removed
    fprintf(stderr,"xye: %15.6E eos_yemin: %15.6E\n",xye,eos_params->eos_yemin);
    return 102;
  }
  if(xtemp > eos_params->eos_tempmax) {
    return 103;
  }
  if(xtemp < eos_params->eos_tempmin) {
    return 104;
  }
  return 0;
}
//------------------------------------------
static inline __attribute__((always_inline))
int NRPyEOS_checkbounds_kt0_noTcheck(const NRPyEOS_params *restrict eos_params,
                                     const double xrho,
                                     const double xye) {

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 105 -- rho too high
  // 106 -- rho too low

  if(xrho > eos_params->eos_rhomax) {
    return 105;
  }
  if(xrho < eos_params->eos_rhomin) {
    return 106;
  }
  if(xye > eos_params->eos_yemax) {
    return 101;
  }
  if(xye < eos_params->eos_yemin) {
    return 102;
  }
  return 0;
}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_get_interp_spots(const NRPyEOS_params *restrict eos_params,
                              const double x,
                              const double y,
                              const double z,
                              double *restrict delx,
                              double *restrict dely,
                              double *restrict delz,
                              int *restrict idx) {

  int ix = 1 + (int)( (x - eos_params->logrho[0]  - 1.0e-10) * eos_params->drhoi  );
  int iy = 1 + (int)( (y - eos_params->logtemp[0] - 1.0e-10) * eos_params->dtempi );
  int iz = 1 + (int)( (z - eos_params->yes[0]     - 1.0e-10) * eos_params->dyei  );

  ix = MAX( 1, MIN( ix, eos_params->nrho -1 ) );
  iy = MAX( 1, MIN( iy, eos_params->ntemp-1 ) );
  iz = MAX( 1, MIN( iz, eos_params->nye  -1 ) );

  idx[0] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*iz    ));
  idx[1] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*iz    ));
  idx[2] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*iz    ));
  idx[3] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[4] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*iz    ));
  idx[5] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[6] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));
  idx[7] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));

  // set up aux vars for interpolation
  *delx = eos_params->logrho[ix]  - x;
  *dely = eos_params->logtemp[iy] - y;
  *delz = eos_params->yes[iz]     - z;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_get_interp_spots_linT_low(const NRPyEOS_params *restrict eos_params,
                                       const double x,
                                       const double y,
                                       const double z,
                                       double *restrict delx,
                                       double *restrict dely,
                                       double *restrict delz,
                                       int *restrict idx) {

  int ix = 1 + (int)( (x - eos_params->logrho[0] - 1.0e-10) * eos_params->drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - eos_params->yes[0]    - 1.0e-10) * eos_params->dyei );

  ix = MAX( 1, MIN( ix, eos_params->nrho-1 ) );
  iz = MAX( 1, MIN( iz, eos_params->nye -1 ) );

  idx[0] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[1] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[2] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[3] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[4] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[5] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[6] = NRPyEOS_ntablekeys*(ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));
  idx[7] = NRPyEOS_ntablekeys*((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));

  // set up aux vars for interpolation
  *delx = eos_params->logrho[ix] - x;
  *dely = eos_params->temp1      - y;
  *delz = eos_params->yes[iz]    - z;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_get_interp_spots_linT_low_eps(const NRPyEOS_params *restrict eos_params,
                                           const double x,
                                           const double y,
                                           const double z,
                                           double *restrict delx,
                                           double *restrict dely,
                                           double *restrict delz,
                                           int *restrict idx) {

  int ix = 1 + (int)( (x - eos_params->logrho[0] - 1.0e-10) * eos_params->drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - eos_params->yes[0]    - 1.0e-10) * eos_params->dyei );

  ix = MAX( 1, MIN( ix, eos_params->nrho-1 ) );
  iz = MAX( 1, MIN( iz, eos_params->nye -1 ) );

  idx[0] = (ix     + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[1] = ((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*iz));
  idx[2] = (ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[3] = (ix     + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[4] = ((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*iz));
  idx[5] = ((ix-1) + eos_params->nrho*(iy     + eos_params->ntemp*(iz-1)));
  idx[6] = (ix     + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));
  idx[7] = ((ix-1) + eos_params->nrho*((iy-1) + eos_params->ntemp*(iz-1)));

  // set up aux vars for interpolation
  *delx = eos_params->logrho[ix] - x;
  *dely = eos_params->temp1      - y;
  *delz = eos_params->yes[iz]    - z;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_linterp_one(const NRPyEOS_params *restrict eos_params,
                         const int *restrict idx,
                         const double delx,
                         const double dely,
                         const double delz,
                         double *restrict f,
                         const int iv) {

  // helper variables
  double fh[8], a[8];

  fh[0] = eos_params->alltables[iv+idx[0]];
  fh[1] = eos_params->alltables[iv+idx[1]];
  fh[2] = eos_params->alltables[iv+idx[2]];
  fh[3] = eos_params->alltables[iv+idx[3]];
  fh[4] = eos_params->alltables[iv+idx[4]];
  fh[5] = eos_params->alltables[iv+idx[5]];
  fh[6] = eos_params->alltables[iv+idx[6]];
  fh[7] = eos_params->alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = eos_params->drhoi       * ( fh[1] - fh[0] );
  a[2] = eos_params->dtempi      * ( fh[2] - fh[0] );
  a[3] = eos_params->dyei        * ( fh[3] - fh[0] );
  a[4] = eos_params->drhotempi   * ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = eos_params->drhoyei     * ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = eos_params->dtempyei    * ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = eos_params->drhotempyei * ( fh[7] - fh[0] + fh[1] + fh[2] +
                                     fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0]
     + a[1] * delx
     + a[2] * dely
     + a[3] * delz
     + a[4] * delx * dely
     + a[5] * delx * delz
     + a[6] * dely * delz
     + a[7] * delx * dely * delz;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_linterp_one_linT_low(const NRPyEOS_params *restrict eos_params,
                                  const int *restrict idx,
                                  const double delx,
                                  const double dely,
                                  const double delz,
                                  double *restrict f,
                                  const int iv) {

  // helper variables
  double fh[8], a[8];

  fh[0] = eos_params->alltables[iv+idx[0]];
  fh[1] = eos_params->alltables[iv+idx[1]];
  fh[2] = eos_params->alltables[iv+idx[2]];
  fh[3] = eos_params->alltables[iv+idx[3]];
  fh[4] = eos_params->alltables[iv+idx[4]];
  fh[5] = eos_params->alltables[iv+idx[5]];
  fh[6] = eos_params->alltables[iv+idx[6]];
  fh[7] = eos_params->alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = eos_params->drhoi          * ( fh[1] - fh[0] );
  a[2] = eos_params->dlintempi      * ( fh[2] - fh[0] );
  a[3] = eos_params->dyei           * ( fh[3] - fh[0] );
  a[4] = eos_params->drholintempi   * ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = eos_params->drhoyei        * ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = eos_params->dlintempyei    * ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = eos_params->drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] +
                                        fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0]
     + a[1] * delx
     + a[2] * dely
     + a[3] * delz
     + a[4] * delx * dely
     + a[5] * delx * delz
     + a[6] * dely * delz
     + a[7] * delx * dely * delz;

}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_linterp_one_linT_low_eps(const NRPyEOS_params *restrict eos_params,
                                      const int *restrict idx,
                                      const double delx,
                                      const double dely,
                                      const double delz,
                                      double *restrict f) {

  // helper variables
  double fh[8], a[8];

  fh[0] = eos_params->epstable[idx[0]];
  fh[1] = eos_params->epstable[idx[1]];
  fh[2] = eos_params->epstable[idx[2]];
  fh[3] = eos_params->epstable[idx[3]];
  fh[4] = eos_params->epstable[idx[4]];
  fh[5] = eos_params->epstable[idx[5]];
  fh[6] = eos_params->epstable[idx[6]];
  fh[7] = eos_params->epstable[idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = eos_params->drhoi          * ( fh[1] - fh[0] );
  a[2] = eos_params->dlintempi      * ( fh[2] - fh[0] );
  a[3] = eos_params->dyei           * ( fh[3] - fh[0] );
  a[4] = eos_params->drholintempi   * ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = eos_params->drhoyei        * ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = eos_params->dlintempyei    * ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = eos_params->drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] +
                                        fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0]
     + a[1] * delx
     + a[2] * dely
     + a[3] * delz
     + a[4] * delx * dely
     + a[5] * delx * delz
     + a[6] * dely * delz
     + a[7] * delx * dely * delz;

}
//------------------------------------------
static inline __attribute__((always_inline))
double NRPyEOS_linterp2D(const double *restrict xs,
                         const double *restrict ys,
                         const double *restrict fs,
                         const double x,
                         const double y) {

  //  2     3
  //
  //  0     1
  //
  // first interpolate in x between 0 and 1, 2 and 3
  // then interpolate in y
  // assume rectangular grid

  double dxi = 1./(xs[1]-xs[0]);
  double dyi = 1./(ys[1]-ys[0]); // x*1./y uses faster instructions than x/y
  double t1 = (fs[1]-fs[0])*dxi * (x - xs[0]) + fs[0];
  double t2 = (fs[3]-fs[2])*dxi * (x - xs[0]) + fs[2];

  return (t2 - t1)*dyi * (y-ys[0]) + t1;
}
//------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_bisection(const NRPyEOS_params *restrict eos_params,
                       const double lr,
                       const double lt0,
                       const double ye,
                       const double leps0,
                       const double prec,
                       double *restrict ltout,
                       const int iv,
                       int *restrict keyerrt) {
  // iv is the index of the variable we do the bisection on

  int bcount = 0;
  int maxbcount = 80;
  int itmax = 50;

  const double dlt0p = log(1.1);
  const double dlt0m = log(0.9);
  const double dltp  = log(1.2);
  const double dltm  = log(0.8);

  double leps0_prec = fabs(leps0*prec);

  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = eos_params->logtemp[0];
  double ltmax = eos_params->logtemp[eos_params->ntemp-1];
  double f1,f2,fmid,dlt,ltmid;
  double f1a = 0.0;
  double f2a = 0.0;
  double delx,dely,delz;
  int idx[8];

  // LSMOD (Modification made by Lorenzo Sala)
  // LSMOD: The following lines calculate eps in
  //        f2a = eps(rho,Tmin, Ye) and f1a = eps(rho,Tmax,Ye)
  NRPyEOS_get_interp_spots(eos_params,lr,ltmax,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f1a,iv);
  NRPyEOS_get_interp_spots(eos_params,lr,ltmin,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

  // prepare
  // check if your energy is actually tabulated at this rho and ye.
  // f2a is the energy evaluated at ltmin, so it is the minimum energy tabulated
  // at this rho ad ye.
  // If leps0 <= f2a, then ltout is likely to be the minimum temperature tabulated.
  if(leps0 <= f2a) { // + 1.0E-6
    *ltout = ltmin;
    return;
  }

  /* // If leps0 >= f1a, then ltout is likely to be the maximum temperature tabulated.
     if(leps0 >= f1a) { // + 1.0E-6
     *ltout = ltmax;
     return;
     } */

  // otherwise, proceed finding extrema for applying bisection method.
  lt = lt0;
  lt1 = MIN(lt0 + dlt0p,ltmax);
  lt2 = MAX(lt0 + dlt0m,ltmin);

  NRPyEOS_get_interp_spots(eos_params,lr,lt1,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f1a,iv);

  NRPyEOS_get_interp_spots(eos_params,lr,lt2,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

  f1=f1a-leps0;
  f2=f2a-leps0;

  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
  while(f1*f2 >= 0.0) {
    lt1 = MIN(lt1 + dltp,ltmax);
    lt2 = MAX(lt2 + dltm,ltmin);
    NRPyEOS_get_interp_spots(eos_params,lr,lt1,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f1a,iv);

    NRPyEOS_get_interp_spots(eos_params,lr,lt2,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

    f1=f1a-leps0;
    f2=f2a-leps0;

#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
            bcount,f1,f2,lt1,lt2,f1a,f2a,leps0);
#endif

    bcount++;
    if(bcount >= maxbcount) {
#if DEBUG
      fprintf(stderr,"bcount out of range it %d, lr: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E, ye: %15.6E\n",
              bcount,lr,lt1,lt2,f1a,f2a,leps0,ye);
#endif
      *keyerrt = 667;
      return;
    }
  } // while

  if(f1 < 0.0) {
    lt = lt1;
    dlt = lt2 - lt1;
  } else {
    lt = lt2;
    dlt = lt1 - lt2;
  }

#if DEBUG
  fprintf(stderr,"bisection step 2 it -1, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
          f2,lt,dlt);
  fprintf(stderr,"ltmax: %15.6E\n",ltmax);
#endif

  int it;
  for(it=0;it<itmax;it++) {
    dlt = dlt * 0.5;
    ltmid = lt + dlt;
    NRPyEOS_get_interp_spots(eos_params,lr,ltmid,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&f2a,iv);

    fmid=f2a-leps0;
    if(fmid <= 0.0) lt=ltmid;
#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E f2a: %15.6E lt: %15.6E ltmid: %15.6E dlt: %15.6E\n",
            it,fmid,f2a,lt,ltmid,dlt);
#endif

    if(fabs(leps0-f2a) <= leps0_prec) {
      *ltout = ltmid;
      return;
    }
  } // for it = 0

  *keyerrt = 667;
  return;
} // bisection
  //------------------------------------------
static inline __attribute__((always_inline))
void NRPyEOS_findtemp_from_any( const NRPyEOS_params *restrict eos_params,
                                const int tablevar_key,
                                const double lr,
                                const double lt0,
                                const double ye,
                                const double tablevar_in,
                                const double prec,
                                double *restrict ltout,
                                int *keyerrt ) {

  // local variables
  const int itmax = 200; // use at most 10 iterations, then go to bisection
  double dtablevardlti; // 1 / derivative dlogeps/dlogT
  double ldt;
  double tablevar; // temp vars for eps
  double ltn; // temp vars for temperature
  const double ltmax = eos_params->logtemp[eos_params->ntemp-1]; // max temp
  const double ltmin = eos_params->logtemp[0]; // min temp
  int it = 0;

  // setting up some vars
  *keyerrt  = 0;
  double lt = lt0;

  // step 1: do we already have the right temperature
  int idx[8];
  double delx,dely,delz;
  NRPyEOS_get_interp_spots(eos_params,lr,lt,ye,&delx,&dely,&delz,idx);
  NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&tablevar,tablevar_key);

  // TODO: profile this to see which outcome is more likely
  if(fabs(tablevar-tablevar_in) < prec*fabs(tablevar_in)) {
    *ltout = lt0;
    return;
  }

  double oerr = 1.0e90;
  double fac  = 1.0;
  const int irho = MIN(MAX(1 + (int)(( lr - eos_params->logrho[0] - 1.0e-12) * eos_params->drhoi),1),eos_params->nrho-1);
  const int iye  = MIN(MAX(1 + (int)(( ye - eos_params->yes[0]    - 1.0e-12) * eos_params->dyei ),1),eos_params->nye -1);

  /* ******* if temp low for high density, switch directly to bisection.
     Verifying Newton-Raphson result evaluating the derivative.
     The variable shouldgotobisection will be modified accordingly
     to the value of derivative of eps wrt temp ******* */
  bool shouldgotobisection = false; // LSMOD
  while(it < itmax && shouldgotobisection == false) {
    it++;

    // step 2: check if the two bounding values of the temperature
    //         give eps values that enclose the new eps.
    const int itemp = MIN(MAX(1 + (int)(( lt - eos_params->logtemp[0] - 1.0e-12) * eos_params->dtempi),1),eos_params->ntemp-1);

    double tablevart1, tablevart2;
    // lower temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye-1)));
      fs[0]   = eos_params->alltables[ifs];
      // point 1
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye-1)));
      fs[1]   = eos_params->alltables[ifs];
      // point 2
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye)));
      fs[2]   = eos_params->alltables[ifs];
      // point 3
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp-1) + eos_params->ntemp*(iye)));
      fs[3]   = eos_params->alltables[ifs];

      tablevart1 = NRPyEOS_linterp2D(&eos_params->logrho[irho-1],&eos_params->yes[iye-1], fs, lr, ye);
    }
    // upper temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp) + eos_params->ntemp*(iye-1)));
      fs[0]   = eos_params->alltables[ifs];
      // point 1
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp) + eos_params->ntemp*(iye-1)));
      fs[1]   = eos_params->alltables[ifs];
      // point 2
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho-1 + eos_params->nrho*((itemp) + eos_params->ntemp*(iye)));
      fs[2]   = eos_params->alltables[ifs];
      // point 3
      ifs     = tablevar_key + NRPyEOS_ntablekeys*(irho   + eos_params->nrho*((itemp) + eos_params->ntemp*(iye)));
      fs[3]   = eos_params->alltables[ifs];

      tablevart2 = NRPyEOS_linterp2D(&eos_params->logrho[irho-1],&eos_params->yes[iye-1], fs, lr, ye);
    }

    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if((tablevar_in - tablevart1) * (tablevar_in - tablevart2) <= 0.) {

      *ltout = (eos_params->logtemp[itemp]-eos_params->logtemp[itemp-1]) / (tablevart2 - tablevart1) *
        (tablevar_in - tablevart1) + eos_params->logtemp[itemp-1];

      return;
    }

    // well, then do a Newton-Raphson step
    // first, guess the derivative
    dtablevardlti = (eos_params->logtemp[itemp]-eos_params->logtemp[itemp-1])/(tablevart2-tablevart1);
    ldt = -(tablevar - tablevar_in) * dtablevardlti * fac;

    //LSMOD: too large a dlt means that the energy dependence on the temperature
    //       is weak ==> We'd better try bisection.
    //       Factor 1/12.0 come from tests by LSMOD
    //       This is done in order to limit the "velocity" of T variation
    //       given by Newton-Raphson.
    if(ldt > (ltmax-ltmin) / 12.0 ) shouldgotobisection = true;

    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt = ltn;

    NRPyEOS_get_interp_spots(eos_params,lr,lt,ye,&delx,&dely,&delz,idx);
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&tablevar,tablevar_key);

    // drive the thing into the right direction
    double err = fabs(tablevar-tablevar_in);
    if(oerr < err) fac *= 0.9;
    oerr = err;

    if(err < prec*fabs(tablevar_in)) {
      *ltout = lt;
      return;
    }

  } // while(it < itmax)

    // try bisection
  NRPyEOS_bisection(eos_params,lr,lt0,ye,tablevar_in,prec,ltout,tablevar_key,keyerrt);

  return;
}
