#include "con2prim_header.h"
#include "EOS_hybrid_header.h"
#include "harm_u2p_util.h"

/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble,
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to
    solve the relativistic magnetohydrodynamic equations of motion on a
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model.

    You are morally obligated to cite the following two papers in his/her
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003,
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006,
        Astrophysical Journal, 641, 626.


    Further, we strongly encourage you to obtain the latest version of
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_2d.c:
---------------

    Uses the 2D method:
       -- solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method
       -- can be used (in principle) with a general equation of state.

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want
      to change this aspect of the code so that it still calculates the
      velocity and so that you can floor the densities.  If you want to
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

#ifdef  NEWT_DIM
#undef  NEWT_DIM
#endif
#define NEWT_DIM (2)

// Declarations:
double vsq_calc(const harm_aux_vars_struct *restrict harm_aux, const double W);

int general_newton_raphson( const eos_parameters *restrict eos, harm_aux_vars_struct *restrict harm_aux,
                            double x[], int n,
                            void (*funcd)(const eos_parameters *restrict, harm_aux_vars_struct *restrict,double [], double [], double [], double [][NEWT_DIM], double *, double *, int) );

void func_vsq( const eos_parameters *restrict eos, harm_aux_vars_struct *restrict harm_aux,double [], double [], double [], double [][NEWT_DIM], double *f, double *df, int n);

double x1_of_x0(const harm_aux_vars_struct *restrict harm_aux, const double x0 ) ;

double pressure_W_vsq(const eos_parameters *restrict eos, const double W, const double vsq, const double D) ;
double dpdW_calc_vsq(const eos_parameters *restrict eos, const double W, const double vsq);
double dpdvsq_calc(const eos_parameters *restrict eos, const double W, const double vsq, const double D);
int Utoprim_new_body(const eos_parameters *restrict eos,
                     const double *restrict U,
                     const double gcov[NDIM][NDIM],
                     const double gcon[NDIM][NDIM],
                     double *restrict prim);

/**********************************************************************/
/******************************************************************

  Utoprim_2d():

  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may
     wish to alter the translation as they see fit.  Note that Greek
     indices run 0,1,2,3 and Latin indices run 1,2,3 (spatial only).


             /     rho u^t     \
         U = | T^t_t + rho u^t |  sqrt(-det(g_{\mu\nu}))
             |      T^t_i      |
             \       B^i       /

             /     rho     \
         P = |     uu      |
             | \tilde{u}^i |
             \     B^i     /


   Arguments:
       U[NPR]           = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR]        = primitive variables (guess on input, calculated values on
                                                output if there are no problems);

   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/

/* some mnemonics */
/* for primitive variables */
static const int UU       =1;
static const int UTCON1   =2;
static const int UTCON2   =3;
static const int UTCON3   =4;
static const int BCON1    =5;
static const int BCON2    =6;
static const int BCON3    =7;

/* for conserved variables */
static const int QCOV0    =1;
//static const int QCOV1    =2;
//static const int QCOV2    =3;
//static const int QCOV3    =4;

/* Also for primitive variables */
static const int RHO      =0;
//static const int EPS      =1;
//static const int v1_con   =2;
//static const int v2_con   =3;
//static const int v3_con   =4;
static const int B1_con   =5;
static const int B2_con   =6;
static const int B3_con   =7;
static const int YE       =8;
static const int TEMP     =9;
//static const int PRESS    =10;
static const int WLORENTZ =11;
//static const int ENT      =12;
static const int numprims =13; // rho, eps, v^{x,y,z}, B^{x,y,z}, Ye, T, P, W, S

/* Also for conservative variables */
static const int DD       =0;
static const int S1_cov   =2;
static const int S2_cov   =3;
static const int S3_cov   =4;
static const int TAU      =9;
static const int WS       =10;
static const int numcons  =11; // D, UU, S_{x,y,z}, B^{x,y,z}, DYe, tau, DS

int con2prim_Noble2D( const eos_parameters *restrict eos,
                      const metric_quantities *restrict metric,
                      const conservative_quantities *restrict cons_undens,
                      primitive_quantities *restrict prims,
                      con2prim_diagnostics *restrict diagnostics ) {

  double uu = -cons_undens->tau*metric->lapse - (metric->lapm1)*cons_undens->rho +
    metric->betax*cons_undens->S_x + metric->betay*cons_undens->S_y  + metric->betaz*cons_undens->S_z;


  double new_cons[numcons];
  double new_prims[numprims];
  new_cons[DD] = cons_undens->rho;
  new_cons[UU] = uu - cons_undens->rho;
  new_cons[S1_cov] = cons_undens->S_x;
  new_cons[S2_cov] = cons_undens->S_y;
  new_cons[S3_cov] = cons_undens->S_z;
  new_cons[B1_con] = prims->Bx * ONE_OVER_SQRT_4PI;
  new_cons[B2_con] = prims->By * ONE_OVER_SQRT_4PI;
  new_cons[B3_con] = prims->Bz * ONE_OVER_SQRT_4PI;
  new_cons[YE] = cons_undens->Y_e;
  new_cons[TAU] = MAX(cons_undens->tau, eos->tau_atm);
  new_cons[WS] = cons_undens->entropy;

  int polytropic_index = find_polytropic_K_and_Gamma_index(eos,prims->rho);
  double Gamma_ppoly_tab = eos->Gamma_ppoly_tab[polytropic_index];

  new_prims[RHO] = prims->rho;
  new_prims[UU] = prims->press/(Gamma_ppoly_tab - 1.0);
  new_prims[UTCON1] = 0.0;
  new_prims[UTCON2] = 0.0;
  new_prims[UTCON3] = 0.0;
  new_prims[BCON1] = new_cons[B1_con];
  new_prims[BCON2] = new_cons[B2_con];
  new_prims[BCON3] = new_cons[B3_con];
  new_prims[WLORENTZ] = 1.0;
  if( eos->eos_type == 1 ) {
    new_prims[TEMP] = prims->temp;
    new_prims[YE] = prims->Y_e;
  }

  int retval = Utoprim_new_body(eos, new_cons, metric->g4dn, metric->g4up, new_prims);

  //TODO: missing eps, entropy. Are these not used for this routine?

  prims->rho = new_prims[RHO];
  if( eos->eos_type == 1 ) {
    prims->temp = new_prims[TEMP];
    prims->Y_e = new_prims[YE];
  }

  limit_velocity_and_convert_utilde_to_v(eos, metric, &new_prims[UTCON1], &new_prims[UTCON2],
                                         &new_prims[UTCON3], cons_undens->rho, prims, diagnostics);

  if( eos->eos_type == 0 )
    prims->press = pressure_rho0_u(eos, prims->rho,new_prims[UU]);

  return retval;
}

/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the
        Newton-Raphson routine.

  -- assumes that
             /   rho gamma   \
         U = | alpha T^t_\mu |
             \   alpha B^i   /



             /     rho     \
      prim = |     uu      |
             | \tilde{u}^i |
             \  alpha B^i  /


return:  (i*100 + j)  where
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used)
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence
                   (occurrence of "nan" or "+/-inf" ;

         j = 0 -> success
             1 -> failure: some sort of failure in Newton-Raphson;
             2 -> failure: utsq<0 w/ initial p[] guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

int Utoprim_new_body( const eos_parameters *restrict eos,
                      const double *restrict U,
                      const double gcov[NDIM][NDIM],
                      const double gcon[NDIM][NDIM],
                      double *restrict prim )
{

  double x_2d[NEWT_DIM];
  double Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qtcon[NDIM];
  double rho0,u,p,w,gammasq,gtmp,W_last,W,utsq,vsq;
  int i,j, n, retval, i_increase;

  // Contains Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,gamma,gamma_times_S,ye
  harm_aux_vars_struct harm_aux;

  n = NEWT_DIM ;

  // Assume ok initially:
  retval = 0;

  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  lower_g(Bcon,gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  raise_g(Qcov,gcon,Qcon) ;


  harm_aux.Bsq = 0. ;
  for(i=1;i<4;i++) harm_aux.Bsq += Bcon[i]*Bcov[i] ;

  harm_aux.QdotB = 0. ;
  for(i=0;i<4;i++) harm_aux.QdotB += Qcov[i]*Bcon[i] ;
  harm_aux.QdotBsq = harm_aux.QdotB*harm_aux.QdotB ;

  ncov_calc(gcon,ncov) ;
  // FIXME: The exact form of n^{\mu} can be found
  //        in eq. (2.116) and implementing it
  //        directly is a lot more efficient than
  //        performing n^{\mu} = g^{\mu\nu}n_{nu}
  raise_g(ncov,gcon,ncon);

  harm_aux.Qdotn = Qcon[0]*ncov[0] ;

  harm_aux.Qsq = 0. ;
  for(i=0;i<4;i++) harm_aux.Qsq += Qcov[i]*Qcon[i] ;

  harm_aux.Qtsq = harm_aux.Qsq + harm_aux.Qdotn*harm_aux.Qdotn ;

  harm_aux.D    = U[RHO];

  /* calculate W from last timestep and use for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) {
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  gammasq = 1. + utsq ;
  harm_aux.gamma  = sqrt(gammasq);

  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) .
  rho0 = harm_aux.D / harm_aux.gamma ;
  u = prim[UU];

  // p = 0.0;
  // if( eos.is_Hybrid ) {
    p = pressure_rho0_u(eos, rho0,u);
  // }
  // else if( eos.is_Tabulated ) {
  //   harm_aux.ye            = U[YE]/U[RHO];
  //   harm_aux.gamma_times_S = U[WS];
  //   harm_aux.use_entropy   = false;
  //   harm_aux.T_guess       = prim[TEMP];
  //   double xrho         = rho0;
  //   double xye          = harm_aux.ye;
  //   double xtemp        = harm_aux.T_guess;
  //   double xprs         = 0.0;
  //   double xeps         = 0.0;
  //   double xdepsdT      = 0.0;

  //   // Now compute P and eps from (rho,Ye,T). Note that
  //   // at this point we do not know W, so we do not
  //   // use the entropy in this function call.
  //   get_P_eps_and_depsdT_from_rho_Ye_and_T( eos,xrho,xye,xtemp, &xprs,&xeps,&xdepsdT );
  //   p = xprs;
  //   u = xeps*xrho;
  //   if( xdepsdT < eos.depsdT_threshold ) harm_aux.use_entropy = true;
  // }

  w = rho0 + u + p ;
  W_last = w*gammasq ;

  // Make sure that W is large enough so that v^2 < 1 :
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*harm_aux.Bsq )
            - harm_aux.QdotBsq*(2.*W_last + harm_aux.Bsq) ) <= W_last*W_last*(harm_aux.Qtsq-harm_aux.Bsq*harm_aux.Bsq))
         && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }

  // Calculate W and vsq:
  x_2d[0] =  fabs( W_last );
  x_2d[1] = x1_of_x0( &harm_aux, W_last ) ;

  retval = general_newton_raphson( eos, &harm_aux, x_2d, n, func_vsq) ;

  W = x_2d[0];
  vsq = x_2d[1];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
      return(retval) ;
    }
  }

  // Calculate v^2:
  if( vsq >= 1. ) {
    vsq = 1.-2.e-16;
    //retval = 4;
    //return(retval) ;
  }


  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  harm_aux.gamma = 1./gtmp ;
  rho0 = harm_aux.D * gtmp;

  w = W * (1. - vsq) ;

  // if( eos.is_Hybrid ) {
    p = pressure_rho0_w(eos, rho0,w) ;
    u = w - (rho0 + p) ; // u = rho0 eps, w = rho0 h
    prim[RHO] = rho0 ;
    prim[UU ] = u ;
  // }
  // else {
  //   double xrho  = rho0;
  //   double xye   = harm_aux.ye;
  //   double xtemp = harm_aux.T_guess;
  //   double xent  = harm_aux.gamma_times_S / W;
  //   double xprs  = 0.0;
  //   double xuu   = 0.0;
  //   double xeps  = 0.0;
  //   if( harm_aux.use_entropy ) {
  //     get_P_eps_and_T_from_rho_Ye_and_S( eos,xrho,xye,xent, &xprs,&xeps,&xtemp );
  //   }
  //   else {
  //     xprs  = -0.5*harm_aux.Bsq/(harm_aux.gamma*harm_aux.gamma)+harm_aux.Qdotn+W+harm_aux.Bsq-0.5*harm_aux.QdotBsq/(W*W);;
  //     xuu   = (W-harm_aux.D*harm_aux.gamma-xprs*harm_aux.gamma*harm_aux.gamma)/(harm_aux.D*harm_aux.gamma) * rho0;
  //     xeps  = xuu/xrho;
  //     get_P_S_and_T_from_rho_Ye_and_eps( eos,xrho,xye,xeps, &xprs,&xent,&xtemp );
  //   }

  //   // Update P and T in the prim array
  //   prim[RHO  ] = xrho;
  //   prim[YE   ] = harm_aux.ye;
  //   prim[TEMP ] = MIN(MAX(xtemp,eos.T_atm),eos.T_max);
  //   prim[PRESS] = xprs;
  //   prim[EPS  ] = xeps;
  //   prim[ENT  ] = xent;
  // }

  if( (rho0 <= 0.) || (u <= 0.) ) {
    // User may want to handle this case differently, e.g. do NOT return upon
    // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:

    retval = 5;
    //return(retval) ;
  }

  for(i=1;i<4;i++) Qtcon[i] = Qcon[i] + ncon[i] * harm_aux.Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = harm_aux.gamma/(W+harm_aux.Bsq) * ( Qtcon[i] + harm_aux.QdotB*Bcon[i]/W ) ;

  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  /* done! */
  return(retval) ;

}



/**********************************************************************/
/****************************************************************************
   vsq_calc():

      -- evaluate v^2 (spatial, normalized velocity) from
            W = \gamma^2 w

****************************************************************************/
double vsq_calc(const harm_aux_vars_struct *restrict harm_aux, const double W)
{
  double Wsq = W*W ;
  double Xsq = (harm_aux->Bsq + W) * (harm_aux->Bsq + W);

  return(  ( Wsq * harm_aux->Qtsq  + harm_aux->QdotBsq * (harm_aux->Bsq + 2.*W)) / (Wsq*Xsq) );
}



/********************************************************************

  x1_of_x0():

    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

double x1_of_x0(const harm_aux_vars_struct *restrict harm_aux, const double x0 )
{
  double dv  = 1.e-15;
  double vsq = fabs(vsq_calc(harm_aux,x0)) ; // guaranteed to be positive

  return( ( vsq > 1. ) ? (1.0 - dv) : vsq   );
}


/********************************************************************

  validate_x():

    -- makes sure that x[0,1] have physical values, based upon
       their definitions:

*********************************************************************/

void validate_x(double x[2], const double x0[2] )
{

  double dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];


  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;

}


/************************************************************

  general_newton_raphson():

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
int general_newton_raphson( const eos_parameters *restrict eos, harm_aux_vars_struct *restrict harm_aux,
                            double x[], int n,
                            void (*funcd)(const eos_parameters *restrict, harm_aux_vars_struct *restrict, double [], double [], double [],
                                          double [][NEWT_DIM], double *, double *, int) )
{
  double f, df, dx[NEWT_DIM], x_old[NEWT_DIM];
  double resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  double errx, x_orig[NEWT_DIM];
  int id, i_extra, doing_extra;
  int n_iter;
  int keep_iterating;


  // Initialize various parameters and variables:
  n_iter = 0;
  errx = 1. ;
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < NEWT_DIM ; id++)  x_old[id] = x_orig[id] = x[id] ;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) {

    (*funcd) (eos, harm_aux, x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */


    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    validate_x( x, x_old ) ;


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/

    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0))
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)

  /*  Check for bad untrapped divergences : */
//TODO: error-checking
//  if( (robust_isfinite(f)==0) ||  (robust_isfinite(df)==0) ) {
//    return(2);
//  }


  if( fabs(errx) > MIN_NEWT_TOL){
    //CCTK_VInfo(CCTK_THORNSTRING,"%d %e %e %e %e",n_iter,f,df,errx,MIN_NEWT_TOL);
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}




/**********************************************************************/
/*********************************************************************************
   func_vsq():

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
*********************************************************************************/

void func_vsq(const eos_parameters *restrict eos, harm_aux_vars_struct *restrict harm_aux, double x[], double dx[], double resid[],
              double jac[][NEWT_DIM], double *f, double *df, int n)
{

  double W   = x[0];
  double vsq = x[1];
  double Wsq = W*W;

  double p_tmp, dPdvsq, dPdW;

  // if( eos.is_Hybrid ) {
    p_tmp  = pressure_W_vsq( eos, W, vsq , harm_aux->D);
    dPdW   = dpdW_calc_vsq( eos, W, vsq );
    dPdvsq = dpdvsq_calc( eos, W, vsq, harm_aux->D );
  // }
  // else {
  //   // Here we must compute dPdW and dPdvsq for tabulated EOS.
  //   // We will be using the expressions
  //   double gamma_sq = 1.0/(1.0-vsq);
  //   harm_aux.gamma     = sqrt(gamma_sq);
  //   if( harm_aux.gamma > eos.W_max ) {
  //     harm_aux.gamma = eos.W_max;
  //     gamma_sq       = harm_aux.gamma * harm_aux.gamma;
  //   }

  //   const double rho = MAX(harm_aux.D / harm_aux.gamma,eos.rho_min);
  //   const double xye = harm_aux.ye;
  //   const double h   = fabs(W /(rho*gamma_sq)); // W := rho*h*gamma^{2}
  //   const double ent = harm_aux.gamma_times_S / harm_aux.gamma;
  //   double T         = harm_aux.T_guess;
  //   double prs       = 0.0;
  //   double eps       = 0.0;
  //   double dPdrho    = 0.0;
  //   double dPdeps    = 0.0;

  //   // Now compute the pressure and its derivatives with respect to rho and eps
  //   if( harm_aux.use_entropy ) {
  //     get_P_eps_T_dPdrho_and_dPdeps_from_rho_Ye_and_S(eos,rho,xye,ent, &prs,&eps,&T,&dPdrho,&dPdeps);
  //   }
  //   else {
  //     get_P_eps_T_dPdrho_and_dPdeps_from_rho_Ye_and_h(eos,rho,xye,h, &prs,&eps,&T,&dPdrho,&dPdeps);
  //   }

  //   // Set P
  //   p_tmp = prs;

  //   // Now compute dP/dW
  //   const double dPdeps_o_rho = dPdeps/rho;
  //   dPdW = ( dPdeps_o_rho/(1.+dPdeps_o_rho) )/gamma_sq;

  //   // And finally dP/d(v^{2})
  //   const double dPdvsq_1 = -0.5*harm_aux.D*harm_aux.gamma*dPdrho;
  //   const double dPdvsq_2 = -0.5*(W + prs*gamma_sq)/rho;
  //   dPdvsq = (dPdvsq_1 + dPdeps*dPdvsq_2)/(1+dPdeps_o_rho);
  // }

  // These expressions were calculated using Mathematica, but made into efficient
  // code using Maple.  Since we know the analytic form of the equations, we can
  // explicitly calculate the Newton-Raphson step:

  double t2  = -0.5*harm_aux->Bsq+dPdvsq;
  double t3  = harm_aux->Bsq+W;
  double t4  = t3*t3;
  double t9  = 1/Wsq;
  double t11 = harm_aux->Qtsq-vsq*t4+harm_aux->QdotBsq*(harm_aux->Bsq+2.0*W)*t9;
  double t16 = harm_aux->QdotBsq*t9;
  double t18 = -harm_aux->Qdotn-0.5*harm_aux->Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  double t21 = 1/t3;
  double t23 = 1/W;
  double t24 = t16*t23;
  double t25 = -1.0+dPdW-t24;
  double t35 = t25*t3+(harm_aux->Bsq-2.0*dPdvsq)*(harm_aux->QdotBsq+vsq*Wsq*W)*t9*t23;
  double t36 = 1/t35;
  double t40 = (vsq+t24)*t3;
  dx[0] = -(t2*t11+t4*t18)*t21*t36;
  dx[1] = -(-t25*t11-2.0*t40*t18)*t21*t36;
  //detJ = t3*t35; // <- set but not used...
  jac[0][0] = -2.0*t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;



  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}



/**********************************************************************
 **********************************************************************

 The following routines specify the equation of state.  All routines
  above here should be indpendent of EOS.  If the user wishes
  to use another equation of state, the below functions must be replaced
  by equivalent routines based upon the new EOS.

  **********************************************************************
  **********************************************************************/


/**********************************************************************/
/**********************************************************************
  pressure_W_vsq():

        -- Hybrid single and piecewise polytropic equation of state;
        -- pressure as a function of P_cold, eps_cold, W, vsq, and D:
**********************************************************************/
double pressure_W_vsq(const eos_parameters *restrict eos, const double W, const double vsq, const double D)
{

  // Compute gamma^{-2} = 1 - v^{2} and gamma^{-1}
  double inv_gammasq = 1.0 - vsq;
  double inv_gamma   = sqrt(inv_gammasq);

  // Compute rho_b = D / gamma
  double rho_b = D*inv_gamma;

  // Compute P_cold and eps_cold
  double P_cold, eps_cold;
  compute_P_cold__eps_cold(eos,rho_b, &P_cold, &eps_cold);

  // Compute p = P_{cold} + P_{th}
  return( ( P_cold + (eos->Gamma_th - 1.0)*( W*inv_gammasq - D*inv_gamma*( 1.0 + eps_cold ) ) )/eos->Gamma_th );

}



/**********************************************************************/
/**********************************************************************
  dpdW_calc_vsq():

      -- partial derivative of pressure with respect to W;
**********************************************************************/
double dpdW_calc_vsq(const eos_parameters *restrict eos, const double W, const double vsq)
{

  return( (eos->Gamma_th - 1.0) * (1.0 - vsq) /  eos->Gamma_th  ) ;

}


/**********************************************************************/
/**********************************************************************
  dpdvsq_calc():

      -- partial derivative of pressure with respect to vsq
**********************************************************************/
double dpdvsq_calc(const eos_parameters *restrict eos, const double W, const double vsq, const double D)
{

  // Set gamma and rho
  double gamma = 1.0/sqrt(1.0 - vsq);
  double rho_b = D/gamma;

  // Compute P_cold and eps_cold
  double P_cold, eps_cold;
  compute_P_cold__eps_cold(eos, rho_b, &P_cold, &eps_cold);

  // Set basic polytropic quantities
  int polytropic_index = find_polytropic_K_and_Gamma_index(eos,rho_b);
  double Gamma_ppoly_tab = eos->Gamma_ppoly_tab[polytropic_index];


  /* Now we implement the derivative of P_cold with respect
   * to v^{2}, given by
   *  ----------------------------------------------------
   * | dP_cold/dvsq = gamma^{2 + Gamma_{poly}/2} P_{cold} |
   *  ----------------------------------------------------
   */
  double dPcold_dvsq = P_cold * pow(gamma,2.0 + 0.5*Gamma_ppoly_tab);


  /* Now we implement the derivative of eps_cold with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------------
   * | deps_cold/dvsq = gamma/(D*(Gamma_ppoly_tab-1)) * (dP_cold/dvsq + gamma^{2} P_cold / 2) |
   *  -----------------------------------------------------------------------------------
   */
  double depscold_dvsq = ( gamma/(D*(Gamma_ppoly_tab-1.0)) ) * ( dPcold_dvsq + 0.5*gamma*gamma*P_cold );

  /* Now we implement the derivative of p_hybrid with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------
   * | dp/dvsq = Gamma_th^{-1}( dP_cold/dvsq                                       |
   * |                          + (Gamma_{th}-1)*(-W                               |
   * |                                            + D gamma (1 + eps_cold)/2       |
   * |                                            - (D/gamma) * deps_cold/dvsq) )  |
   *  -----------------------------------------------------------------------------
   */
  return( ( dPcold_dvsq + (eos->Gamma_th-1.0)*( -W + D*gamma*(1+eps_cold)/2.0 - D*depscold_dvsq/gamma ) )/eos->Gamma_th );
}


/******************************************************************************
             END   OF   UTOPRIM_2D.C
******************************************************************************/
