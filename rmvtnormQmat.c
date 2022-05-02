#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

static void compute_std_vec(double *Qmat,int dim,double *sdvec,double *Qdiaginv);
static double compute_conditional_mean(double *mu,double *x,double *Qrowi,double *Qdiaginv,int dim,int skip);
static void mvtnormgibbsprec(double *Qmat,double *mu,double *Hdiaginv,double *sdvec,double *lower,double *upper,int dim,int nburn,int nthin,double *x);
static void mvtnormgibbs_Sparseprec(double *nnzval,int *nnzInRow,int *colInd,double *mu,double *Hdiaginv,double *sdvec,double *lower,\
                                    double *upper,int dim,int niter,double *x);
                                    
static double rtnorm(double mu,double sig,double lb,double ub);                                    
static double rtnorm_ltrunc(double mu,double sig,double lb);
static double rtnorm_rtrunc(double mu,double sig,double ub);
static double rtnorm_left_truncate(double lb);
static double rtnorm_right_truncate(double ub);                                    

SEXP rmvtnormSparseQmat(SEXP Muvec,SEXP pslot,SEXP islot,SEXP xslot,SEXP lb,SEXP ub,SEXP nburn,SEXP nthin){
    int     i,j,DIM,NBURN,NTHIN,*colInd,*nnzInRow,Niter;
    double  *Muvec_ptr,*nnzval,*lb_ptr,*ub_ptr,*sdvec,*Qdiaginv,*out_ptr;
    SEXP    out;

    
    Muvec_ptr=REAL(Muvec);
    nnzval=REAL(xslot);
    colInd=INTEGER(islot);
    nnzInRow=INTEGER(pslot);
  
    DIM=length(lb);
    lb_ptr=REAL(lb);
    ub_ptr=REAL(ub);
    NBURN=asInteger(nburn);
    NTHIN=asInteger(nthin);
    Niter=NBURN+NTHIN;
    
    
    sdvec=(double *) malloc(DIM*sizeof(double));
    Qdiaginv=(double *) malloc(DIM*sizeof(double));
    if (sdvec==NULL || Qdiaginv==NULL){
        Rprintf("In rmvtnormSparseQmat subroutine, CANNOT allocate memory for sdvec or Qdiaginv.");
        exit(EXIT_FAILURE);
    }
    
    
    PROTECT(out=allocVector(REALSXP,DIM));
    out_ptr=REAL(out);
    
    for (i=0;i<DIM;i++){
         /*------Precalculation------*/
         for (j=nnzInRow[i];j<nnzInRow[i+1];j++){
    	      if (colInd[j]==i){
    	     	  Qdiaginv[i]=1.0/nnzval[j];
    	     	  sdvec[i]=sqrt(Qdiaginv[i]);
    	     	  break;	
    	      }
    	 }
    	 
    	 /*----Set initial value----*/
         if (R_FINITE(lb_ptr[i]) && R_FINITE(ub_ptr[i])){
             out_ptr[i]=.5*(lb_ptr[i]+ub_ptr[i]);
         }
         else if(R_FINITE(lb_ptr[i])){
                 out_ptr[i]=lb_ptr[i]+1.0;
         }
         else if(R_FINITE(ub_ptr[i])){
                 out_ptr[i]=ub_ptr[i]-1.0;
         }
         else{
              out_ptr[i]=0.0;
         } 
    }
     
    GetRNGstate();       
    mvtnormgibbs_Sparseprec(nnzval,nnzInRow,colInd,Muvec_ptr,Qdiaginv,sdvec,lb_ptr,ub_ptr,DIM,Niter,out_ptr);
    PutRNGstate();  
      
    free(sdvec);
    free(Qdiaginv);
    UNPROTECT(1);
    
    return(out);
}


static void mvtnormgibbs_Sparseprec(double *nnzval,int *nnzInRow,int *colInd,double *mu,double *Hdiaginv,\
                                    double *sdvec,double *lower,double *upper,int dim,int niter,double *x){
   int    i,j,k,cind;
   double mu_i,sd_i;

   for (j=0;j<niter;j++){
        /*---Loop over all dimension.---*/
        for (i=0;i<dim;i++){
        	 mu_i=0.0;
        	 for (k=nnzInRow[i];k<nnzInRow[i+1];k++){
        	      cind=colInd[k];
        	      mu_i-=nnzval[k]*(x[cind]-mu[cind]);
        	 }
        	 mu_i*=Hdiaginv[i];
             mu_i+=x[i];
             sd_i=sdvec[i];
             x[i]=rtnorm(mu_i,sd_i,lower[i],upper[i]);
        }
   }
}


/*The rtnorm subroutine samples from a normal distribution truncated from below
 *or above (but not both).
 */
static double rtnorm(double mu,double sig,double lb,double ub){
   double x;
   
   if (lb==R_NegInf){
       x=rtnorm_rtrunc(mu,sig,ub);
   }
   else if (ub==R_PosInf){
            x=rtnorm_ltrunc(mu,sig,lb);
   }
   else {
         x=NA_REAL;
   }
   return x;
}


static double rtnorm_ltrunc(double mu,double sig,double lb){
   double lb_star, z, siginv, res;
   
   siginv=1.0/sig;
   lb_star=(lb-mu)*siginv;
   z=rtnorm_left_truncate(lb_star);
   res=mu+sig*z;
   
   return res;
}


static double rtnorm_rtrunc(double mu,double sig,double ub){
   double ub_star, z, siginv, res;
   
   siginv=1.0/sig;
   ub_star=(ub-mu)*siginv;
   z=rtnorm_right_truncate(ub_star);
   res=mu+sig*z;
   
   return res;
}


static double rtnorm_left_truncate(double lb){
    double  u,v,phia,tmp,x,lbsq;

    if (lb<4.0){
        u=runif(0.0,1.0);
        phia=pnorm(lb,0.0,1.0,1,0);
        tmp=phia+(1.0-phia)*u;
        x=qnorm(tmp,0.0,1.0,1,0);
    } 
    else{
         lbsq=lb*lb;
         do {
             u=runif(0.0,1.0);
             tmp=lbsq-2.0*log(u);
             x=sqrt(tmp);
             v=lb*runif(0.0,1.0);
         } while(v>x);
    }
    
    return x;
}


static double rtnorm_right_truncate(double ub){
    double zm, minus_ub, y;
    
    minus_ub=-1.0*ub;
    zm=rtnorm_left_truncate(minus_ub);
    y=-1.0*zm;

    return y;
}




