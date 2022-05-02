Construct_Wlist_KNN_Euclidean<-function(cormat,Kvec){
  
  nloc<-nrow(cormat);
  nKvec<-length(Kvec);
  Wlist<-vector("list",nKvec);
  lindex<-1:nloc;
  
  nCore<-detectCores();
  cl<-makeCluster(nCore);
  registerDoParallel(cl);
  
  Wlist<-foreach(i=1:nKvec) %dopar% {
    kval<-Kvec[i];
    tmp<-NULL;
    for (j in 1:nloc){
         cindex<-lindex[-j];
         dvec<-rowSums(sweep(cormat[cindex,],2,cormat[j,],'-')^2);
         kneighbor<-cindex[order(dvec)[1:kval]];
         tmp<-rbind(tmp,cbind(j,kneighbor));
    }
    #Make neighboring relationship symmetric
    temp<-rbind(tmp,tmp[,c(2,1)]);
    tp<-temp[!duplicated(temp),];
    tpsort<-tp[order(tp[,1],tp[,2]),];
    
    #Prepare spare form of the spatial weights W given kval
    rowid.freqtab<-as.data.frame(table(tpsort[,1]));
    nnz<-dim(rowid.freqtab)[1];
    rowid.nfreq<-rowid.freqtab[,2];
    
    weights<-sapply(1:nnz,function(i) {
                          rep(x=1.0/rowid.nfreq[i],times=rowid.nfreq[i])
    });
    
    
    Ws<-cbind(tpsort,unlist(weights));
    colnames(Ws)<-c("rowid","colid","weights");
    Ws
  }
  stopCluster(cl);   

  return(Wlist)
}


SARProbit_SimData_Sparse<-function(X,Wsparse,beta,rho,sigma2=1.0){
    
    if (!is.matrix(X) || !is.vector(beta) || ncol(X)!=length(beta)){
        stop("In SARProbit_SimData_Sparse.R, X should be a matrix and beta should be a vector of compatible size.\n");
    }
    
    n<-nrow(X);
    Xbeta<-X%*%beta;
    err<-sqrt(sigma2)*rnorm(n);
 
    xvec<-Xbeta+err;
    In<-sparseMatrix(i=1:n,j=1:n,x=1.0);
    Amat<-In-rho*Wsparse;
    tmp<-ApproximateAinvXvec(Amat,xvec);
    ystar<-as.vector(tmp);
    yout<-ifelse(ystar>0.0,1L,0L);
    
    list(yout=yout,ystar=ystar);
}


SARLogit_SimData_Sparse<-function(X,Wsparse,beta,rho){
    
    if (!is.matrix(X) || !is.vector(beta) || ncol(X)!=length(beta)){
        stop("In SARLogit_SimData_Sparse.R, X should be a matrix and beta should be a vector of compatible size.\n");
    }
    
    n<-nrow(X);
    Xbeta<-X%*%beta;
    err<-rlogis(n);
 
    xvec<-Xbeta+err;
    In<-sparseMatrix(i=1:n,j=1:n,x=1.0);
    Amat<-In-rho*Wsparse;
    tmp<-ApproximateAinvXvec(Amat,xvec);
    ystar<-as.vector(tmp);
    yout<-ifelse(ystar>0.0,1L,0L);
    
    list(yout=yout,ystar=ystar);
}


#Note that Qmat is a dgCMatrix object
rmvtnormSparseQmat<-function(Muvec,Qmat,lbvec,ubvec,nburn=50,nthin=1){
   
   #load the C code
   if (!is.loaded("rmvtnormSparseQmat")){
       lib.file<-file.path(paste("rmvtnormQmat",.Platform$dynlib.ext,sep=""));
       dyn.load(lib.file);
       cat(" -Loaded ",lib.file,"\n");
   }
   
   res<-.Call("rmvtnormSparseQmat",as.vector(Muvec),as.integer(Qmat@p),as.integer(Qmat@i),Qmat@x,as.vector(lbvec),as.vector(ubvec),as.integer(nburn),as.integer(nthin));
   return(res);
}


ApproximateAinvXvec<-function(Amat,xvec){
      
   out<-solve(qr(Amat),xvec);

}



ApproximateLogdetOfA<-function(TraceVec,rho,nTrunc){
    oneToQ<-1:nTrunc;
    rhotmp<-rho^oneToQ;
    logdet<--sum(rhotmp*TraceVec[oneToQ]/oneToQ);
    
    return(logdet)
}


compute_trWpowers<-function(tWsparse,nTrunc){
   
   #load the C code
   if (!is.loaded("compute_trWpowers")){
       lib.file<-file.path(paste("compute_trWpowers",.Platform$dynlib.ext,sep=""));
       dyn.load(lib.file);
       cat(" -Loaded ",lib.file,"\n");
   }
   
   res<-.Call("compute_trWpowers",as.integer(tWsparse@p),as.integer(tWsparse@i),tWsparse@x,nTrunc);
   
   return(res)
}






