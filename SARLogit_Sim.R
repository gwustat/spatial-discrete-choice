#**************************************************#
# Step 1) Setup working directory and source file  #
#         search path.                             #
#**************************************************#
library(xtable);
library(Matrix);
library(doParallel);
library(foreach);
library(geosphere);
source("Utility.R");
source("SARLogit_MCMC_Sparse.R");

#*****************************#
# Step 3) Simulate data       #
#*****************************#
beta_true<-c(1.2,-0.9,0.6,-0.7,0.8);
p<-length(beta_true);
sigma2_true<-1.0;
rho_true<-0.5;
nTrunc<-50;
Kvec<-8;

#****************************************#
# Step 4) MCMC and Prior specification   #
#****************************************#
mubeta_prior<-rep(x=0.0,times=p);
Sbeta_prior<-100*diag(p);
Qbeta_prior<-.01*diag(p);
Qmubeta_prior<-Qbeta_prior%*%mubeta_prior;
rho_lb<-0.0;
rho_ub<-0.99;
rhovec<-seq(from=rho_lb,to=rho_ub,length.out=100);
nsave<-40000;
nburn<-2000;
nthin<-1;

nCore<-5; #detectCores();
nvec<-c(2000,4000,6000,8000);
nsim<-length(nvec);
timeres<-array(data=NA,dim=c(nsim,3));
colnames(timeres)<-c("size","MCMC","pMCMC");

for (k in 1:nsim){

     #*************************************************#
     # Step 2) Construct a list of sparse W matrices   #
     #*************************************************#
     n<-nvec[k];

     xcoord<-runif(n);
     ycoord<-runif(n)
     cormat<-cbind(xcoord,ycoord);
     Wslist<-Construct_Wlist_KNN_Euclidean(cormat,Kvec);
     Wlist<-Wslist[[1]];

     #Make a sparse weights matrix out of W
     rowIndex<-Wlist[,1];
     colIndex<-Wlist[,2];
     valnnz<-Wlist[,3];
     Wsparse<-sparseMatrix(i=rowIndex,j=colIndex,x=valnnz,dims=c(n,n));
     
     #*************************************************#
     # Step 3) Simulate predictors                     #
     #*************************************************#
     X<-matrix(data=rnorm(n*p),nrow=n,ncol=p);
     
     dsim<-SARLogit_SimData_Sparse(X,Wsparse,beta_true,rho_true);
     ystar_true<-dsim$ystar;
     y<-dsim$yout;


     
     #*******************************#
     # Step 4) Run MCMC              #
     #*******************************#
     time.begin<-Sys.time();
     res<-SARLogit_MCMC_Sparse(X,y,Wsparse,mubeta_prior,Qbeta_prior,rhovec,nTrunc,nburn,nthin,nsave);
     time.end<-Sys.time();
     time.MCMC<-as.numeric(difftime(time.end,time.begin,units="secs"));
     timeres[k,1:2]<-c(n,time.MCMC);

     #Compute effective sample size
     neff<-coda::effectiveSize(res);   
     names(neff)<-c(paste("beta",1:p,sep=""),"rho");
     cat("MCMC: effective sample size: \n");
     print(neff,digits=2)

     MC.post.stat<-cbind(colMeans(res),apply(res,2,sd),t(apply(res,2,function(x) quantile(x,prob=c(.025,.5,.975)))),c(beta_true,rho_true));
     colnames(MC.post.stat)<-c("post_mean","sd","q.025","q.50","q.975","truth");
     cat("MCMC: post statistic summary: \n");
     print(MC.post.stat,digits=4)
     
     time.begin<-Sys.time();
     sarlogitres<-SARLogit_pMCMC_Sparse(X,y,Wsparse,mubeta_prior,Qbeta_prior,rhovec,nTrunc,nburn,nthin,nsave,nCore);
     time.end<-Sys.time();
     time.sarlogitP<-as.numeric(difftime(time.end,time.begin,units="secs")); 
     timeres[k,3]<-time.sarlogitP;
     save(sarlogitres,time.sarlogitP,file=paste(paste("./Simulation/SARLogitP",n,sep="_"),"MCMC.RData",sep="_"));

     #Compute effective sample size
     neffp<-coda::effectiveSize(sarlogitres);   
     names(neffp)<-c(paste("beta",1:p,sep=""),"rho");
     cat("Parallel MCMC: effective sample size: \n");
     print(neffp,digits=2)
     
     sarlogit.post.stat<-cbind(colMeans(sarlogitres),apply(sarlogitres,2,sd),
                               t(apply(sarlogitres,2,function(x) quantile(x,prob=c(.025,.5,.975)))),c(beta_true,rho_true));
     colnames(sarlogit.post.stat)<-c("post_mean","sd","q.025","q.50","q.975","truth");  
     cat("Parallel MCMC: post statistic summary: \n");
     print(sarlogit.post.stat,digits=4) 
 
     #Call gc() function to force garbage collector to return memory to the system
     #Use invisible function to suppress output of gc()
     invisible(gc());
}




























