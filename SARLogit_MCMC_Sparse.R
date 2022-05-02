SARLogit_MCMC_Sparse<-function(X,y,Wsparse,mubeta_prior,Qbeta_prior,rhovec,nTrunc,nburn,nthin,nsave){

   n<-nrow(X);
   p<-ncol(X);
   nrho<-length(rhovec);
   tX<-t(X);
   avec<-ifelse(y>0,0,-Inf);
   bvec<-ifelse(y>0,Inf,0);
   
   nu<-7.3;
   sigma2<-(pi^2)*(nu-2.0)/(3.0*nu);
   sigma<-sqrt(sigma2);
   Qmubeta_prior<-Qbeta_prior%*%mubeta_prior;
   
   
   #Compute trace vector of W powers and approximate log determinants
   tWsparse<-t(Wsparse);
   Wpower.trace<-compute_trWpowers(tWsparse,nTrunc);
   logdet<-sapply(1:nrho,function(i){
                         ApproximateLogdetOfA(Wpower.trace,rhovec[i],nTrunc)
                         });
   
   niter<-nburn+nthin*nsave;
   rho<-.5*sample(x=rhovec,size=1);
   
   #Identify matrix as a lgCMatrix class
   In<-sparseMatrix(i=1:n,j=1:n,x=Inf);
   Amat<-In-rho*Wsparse;
   dind<-which(is.infinite(Amat@x));
   offind<-which(!is.infinite(Amat@x));
   Amat@x[dind]<-1.0;   
   
   beta<-mubeta_prior+sqrt(diag(Qbeta_prior))*rnorm(p);
   res<-array(data=NA,dim=c(nsave,p+1));   
   
   lamvec<-rgamma(n=n,shape=nu/2.0,scale=2.0/nu);
   Ldiag<-sparseMatrix(i=1:n,j=1:n,x=Inf);
   Ldind<-which(is.infinite(Ldiag@x));
   
   isave<-1;
   for (i in 1:niter) {
        Ldiag@x[Ldind]<-lamvec;
        Amat@x[offind]<--rho*Wsparse@x;
        ystar<-SARLogit_update_ystar_sparse(Amat,Ldiag,sigma2,X,beta,avec,bvec);
     
        Wystar<-Wsparse%*%ystar;
        rho<-SARLogit_rho_grid_sparse(lamvec,sigma2,X,beta,ystar,Wystar,rhovec,nrho,logdet);
             
        beta<-SARLogit_update_beta_sparse(lamvec,sigma2,p,X,tX,ystar,Wystar,rho,Qbeta_prior,Qmubeta_prior);
        
        lamvec<-SARLogit_update_lamvec(n,nu,ystar,Wystar,X,beta,rho,sigma);

        icount<-i-nburn;
        if (icount>0 && (icount%%nthin==0)){
            res[isave,]<-c(beta,rho);
            isave<-isave+1;
        }
        #cat("iteration= ",i,"\n");
   }  
   
   return(res)
}

SARLogit_pMCMC_Sparse<-function(X,y,Wsparse,mubeta_prior,Qbeta_prior,rhovec,nTrunc,nburn,nthin,nsave,nCore=5){

   if (nCore<1){
       stop("In SARLogit_pMCMC_Sparse subroutine, nCore should be positive!\n");
   }

   n<-nrow(X);
   p<-ncol(X);
   nrho<-length(rhovec);
   tX<-t(X);
   avec<-ifelse(y>0,0,-Inf);
   bvec<-ifelse(y>0,Inf,0);
   
   nu<-7.3;
   sigma2<-(pi^2)*(nu-2.0)/(3.0*nu);
   sigma<-sqrt(sigma2);
   Qmubeta_prior<-Qbeta_prior%*%mubeta_prior;
   
   #Identify matrix as a lgCMatrix class
   In<-sparseMatrix(i=1:n,j=1:n,x=Inf);
   rho.tmp<-.5*sample(x=rhovec,size=1);
   Amat.tmp<-In-rho.tmp*Wsparse;
   dind<-which(is.infinite(Amat.tmp@x));
   offind<-which(!is.infinite(Amat.tmp@x));
   
   
   #Compute trace vector of W powers and approximate log determinants
   tWsparse<-t(Wsparse);
   Wpower.trace<-compute_trWpowers(tWsparse,nTrunc);
   
   if (nsave%%nCore==0){
       avgsave<-nsave/nCore;
       avgrun<-nburn+nthin*avgsave;
       niter<-rep(x=avgrun,times=nCore);
       insave<-rep(x=avgsave,times=nCore);
   }else {
          avgsave<-as.integer(nsave/nCore);
          avgrun<-nburn+nthin*avgsave;
          remsave<-nsave-avgsave*(nCore-1);
          remrun<-nburn+nthin*remsave;
          niter<-c(rep(x=avgrun,times=nCore-1),remrun);
          insave<-c(rep(x=avgsave,times=nCore-1),remsave);
   }
   

   #Call gc() function to force garbage collector to return memory to the system
   #Use invisible function to suppress output of gc()
   rm(tWsparse,Amat.tmp);
   invisible(gc());     
   
   exportfun<-c("SARLogit_update_ystar_sparse","SARLogit_rho_grid_sparse","SARLogit_update_beta_sparse",
                "SARLogit_update_lamvec","rmvtnormSparseQmat","ApproximateAinvXvec");
   
   cl<-makeCluster(nCore);
   registerDoParallel(cl);
   logdet<-foreach(i=1:nrho,.export=c("ApproximateLogdetOfA"),.combine=c) %dopar% {
                   ApproximateLogdetOfA(Wpower.trace,rhovec[i],nTrunc)
   }
   
   pres<-foreach(i=1:nCore,.export=exportfun,.packages=c("Matrix","MASS"),.combine=rbind) %dopar% {
                 set.seed(i);
                 rho<-.5*sample(x=rhovec,size=1);
                 Amat<-In-rho*Wsparse;
                 Amat@x[dind]<-1.0;  
                 Ldiag<-In;
                 Ldind<-which(is.infinite(Ldiag@x));
                
                 beta<-mubeta_prior+sqrt(diag(Qbeta_prior))*rnorm(p);
                 lamvec<-rgamma(n=n,shape=nu/2.0,scale=2.0/nu);
                 res<-array(data=NA,dim=c(insave[i],p+1));   
                
                 isave<-1;
                 for (j in 1:niter[i]) {
                      Ldiag@x[Ldind]<-lamvec;
                      Amat@x[offind]<--rho*Wsparse@x;
                      ystar<-SARLogit_update_ystar_sparse(Amat,Ldiag,sigma2,X,beta,avec,bvec);
     
                      Wystar<-Wsparse%*%ystar;
                      rho<-SARLogit_rho_grid_sparse(lamvec,sigma2,X,beta,ystar,Wystar,rhovec,nrho,logdet);
             
                      beta<-SARLogit_update_beta_sparse(lamvec,sigma2,p,X,tX,ystar,Wystar,rho,Qbeta_prior,Qmubeta_prior);
        
                      lamvec<-SARLogit_update_lamvec(n,nu,ystar,Wystar,X,beta,rho,sigma);

                      icount<-j-nburn;
                      if (icount>0 && (icount%%nthin==0)){
                          res[isave,]<-c(beta,rho);
                          isave<-isave+1;
                      }
                 }
                 res  
   } 
   stopCluster(cl);
     
   return(pres)
}



SARLogit_update_beta_sparse<-function(lamvec,sigma2,p,X,tX,ystar,Wystar,rho,Qbeta_prior,Qmubeta_prior){

    lamdivsig2<-lamvec/sigma2;
    Aystar<-ystar-rho*Wystar;

    tmp<-sweep(X,1,lamdivsig2,'*');
    Qbeta_post<-tX%*%tmp+Qbeta_prior;
    abeta_post<-crossprod(tmp,Aystar)+Qmubeta_prior;
    
    Qbeta_post.cholR<-chol(Qbeta_post);
    beta<-backsolve(Qbeta_post.cholR,backsolve(Qbeta_post.cholR,abeta_post,transpose=TRUE)+rnorm(p));
    
    return(beta)
}


SARLogit_update_ystar_sparse<-function(Amat,Ldiag,sigma2,X,beta,ystar_lb,ystar_ub){

    Xbeta<-X%*%beta;
    muparm<-ApproximateAinvXvec(Amat,Xbeta);
    Qparm<-(1.0/sigma2)*(t(Amat)%*%(Ldiag%*%Amat));
    
    ystar<-rmvtnormSparseQmat(muparm,Qparm,ystar_lb,ystar_ub,nburn=10,nthin=1);
    
    return(ystar)
}


SARLogit_rho_grid_sparse<-function(lamvec,sigma2,X,beta,ystar,Wystar,rhovec,nrho,logdet){
 
   lamdivsig2<-lamvec/sigma2;
   ystarmXbeta<-as.vector(ystar-X%*%beta);
   parta<-sum(Wystar*(lamdivsig2*ystarmXbeta));
   partb<-0.5*sum(lamdivsig2*(Wystar^2));
   rhollike<-logdet+0.5*sum(log(lamdivsig2))+rhovec*parta-rhovec^2*partb;
   
   rhollike.max<-max(rhollike);
   rho.weight.shift<-exp(rhollike-rhollike.max);
   
   rho<-sample(x=rhovec,size=1,prob=rho.weight.shift);
   
   return(rho)
}


SARLogit_update_lamvec<-function(n,nu,ystar,Wystar,X,beta,rho,sigma){

      Aystar<-ystar-rho*Wystar;
      Xbeta<-X%*%beta;
      zvec<-drop(Aystar-Xbeta)/sigma;
      sh<-0.5*(nu+1.0);
      sc<-2.0/(nu+zvec^2);
      lamvec<-sapply(1:n,function(i) {
                              rgamma(n=1,shape=sh,scale=sc[i])
                         });

      return(lamvec)
}
