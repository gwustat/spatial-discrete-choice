#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


static void compute_WsparseVec_Multiply(double *nnzval,int *nnzInRow,int *colInd,int dim,double *vec,double *out);
static void compute_VecWsparse_Multiply(double *nnzval,int *nnzInRow,int *rowInd,int dim,double *vec,double *out);

SEXP compute_trWpowers(SEXP pslot,SEXP islot,SEXP xslot,SEXP nTrunc){
    int     i,j,k,DIM,*colInd,*nnzInRow,maxOrder,Kmax;
    double  *nnzval,*unitvec,*Wunitvec,*tmp,*out_ptr;
    SEXP    out;

    nnzval=REAL(xslot);
    colInd=INTEGER(islot);
    nnzInRow=INTEGER(pslot);
    Kmax=asInteger(nTrunc);
  
    DIM=length(pslot)-1;
    maxOrder=2*Kmax+2;
    
    unitvec=(double *) calloc(DIM,sizeof(double));
    Wunitvec=(double *) malloc(DIM*sizeof(double));
    tmp=(double *) malloc(DIM*sizeof(double));
    if (unitvec==NULL || Wunitvec==NULL || tmp==NULL){
        Rprintf("In compute_trWpowers subroutine, CANNOT allocate sufficient memory.");
        exit(EXIT_FAILURE);
    }
    
    PROTECT(out=allocVector(REALSXP,maxOrder));
    out_ptr=REAL(out);
    
    /*---Initialization.---*/
    for (i=0;i<maxOrder;i++){
         out_ptr[i]=0.0;
    }
    
    /*---Compute trace of W powers.---*/
    for (i=0;i<DIM;i++){
         unitvec[i]=1.0;
         compute_WsparseVec_Multiply(nnzval,nnzInRow,colInd,DIM,unitvec,Wunitvec);
         
         for (j=1;j<maxOrder;j++){
              compute_WsparseVec_Multiply(nnzval,nnzInRow,colInd,DIM,Wunitvec,tmp);
              
              /*Copy tmp to Wunitvec.*/
              for (k=0;k<DIM;k++){
                   Wunitvec[k]=tmp[k];
              }
              out_ptr[j]+=tmp[i];
         }
         unitvec[i]=0.0;
    }
    
    free(tmp);
    free(unitvec);
    free(Wunitvec);
    UNPROTECT(1);
    
    return(out);
}


static void compute_WsparseVec_Multiply(double *nnzval,int *nnzInRow,int *colInd,int dim,double *vec,double *out){
   int    i,j,cind;
   double res; 
   
   /*---Loop over all dimension.---*/
   for (i=0;i<dim;i++){
        res=0.0;
        for (j=nnzInRow[i];j<nnzInRow[i+1];j++){
        	 cind=colInd[j];
        	 res+=nnzval[j]*vec[cind];
        }
        out[i]=res;
   }
}


static void compute_VecWsparse_Multiply(double *nnzval,int *nnzInCol,int *rowInd,int dim,double *vec,double *out){
   int    i,j,rind;
   double res; 
   
   /*---Loop over all dimension.---*/
   for (i=0;i<dim;i++){
        res=0.0;
        for (j=nnzInCol[i];j<nnzInCol[i+1];j++){
        	 rind=rowInd[j];
        	 res+=nnzval[j]*vec[rind];
        }
        out[i]=res;
   }
}

    



    
    
    
    
    
    
