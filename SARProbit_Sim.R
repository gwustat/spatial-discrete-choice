#**************************************************#
# Step 1) Setup working directory and source file  #
#         search path.                             #
#**************************************************#
library(Matrix);
library(doParallel);
library(foreach);
library(geosphere);
library(spatialprobit);
source("Utility.R");

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

nCore<-5; 
insave<-nsave/nCore;
nvec<-c(2000,4000,6000,8000);
nsim<-length(nvec);
timeres<-array(data=NA,dim=c(nsim,3));
colnames(timeres)<-c("size","Pkg","Pkg-Parallel");

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
     
     dsim<-SARProbit_SimData_Sparse(X,Wsparse,beta_true,rho_true);
     ystar_true<-dsim$ystar;
     y<-dsim$yout;
     
     #*******************************#
     # Step 4) Run MCMC              #
     #*******************************#
     time.begin<-Sys.time();
     restmp<-sar_probit_mcmc(y,X,Wsparse,nsave,nburn,thinning=nthin,prior=list(rmin=rho_lb,rmax=rho_ub,c=mubeta_prior,T=Sbeta_prior),computeMarginalEffects=FALSE);
     time.end<-Sys.time();
     res<-cbind(restmp$bdraw,restmp$pdraw);
     time.MCMC<-as.numeric(difftime(time.end,time.begin,units="secs"));
     timeres[k,1:2]<-c(n,time.MCMC);
     
     #Compute effective sample size
     neff<-coda::effectiveSize(res);
     names(neff)<-c(paste("beta",1:p,sep=""),"rho");
     
     cat("sarprobit pkg: effective sample size from MCMC run \n");
     print(neff,digits=2)
     
     
     MC.post.stat<-cbind(colMeans(res),apply(res,2,sd),t(apply(res,2,function(x) quantile(x,prob=c(.025,.5,.975)))),c(beta_true,rho_true));
     colnames(MC.post.stat)<-c("post_mean","sd","q.025","q.50","q.975","truth");
     cat("sarprobit pkg: post statistic summary from MCMC run \n");
     print(MC.post.stat,digits=4)

     time.begin<-Sys.time();
     cl<-makeCluster(nCore);
     registerDoParallel(cl);
     sarprobitres<-foreach(i=1:nCore,.export=c('sar_probit_mcmc'),.packages=c('spatialprobit'),.combine=rbind) %dopar% { 
            res<-sar_probit_mcmc(y,X,Wsparse,insave,nburn,thinning=nthin,prior=list(rmin=rho_lb,rmax=rho_ub,c=mubeta_prior,T=Sbeta_prior),computeMarginalEffects=FALSE);
            cbind(res$bdraw,res$pdraw)
     }
     stopCluster(cl);
     time.end<-Sys.time();
     time.sarprobitP<-as.numeric(difftime(time.end,time.begin,units="secs")); 
     timeres[k,3]<-time.sarprobitP;
     save(sarprobitres,time.sarprobitP,file=paste(paste("./Simulation/SARProbitP",n,sep="_"),"MCMC.RData",sep="_"));
     
     #Compute effective sample size
     neffp<-coda::effectiveSize(sarprobitres);    
     names(neffp)<-c(paste("beta",1:p,sep=""),"rho");
     cat("sarprobit pkg: effective sample size from pMCMC run \n");
     print(neffp,digits=2)
     
     sarprobit.post.stat<-cbind(colMeans(sarprobitres),apply(sarprobitres,2,sd),
                                t(apply(sarprobitres,2,function(x) quantile(x,prob=c(.025,.5,.975)))),c(beta_true,rho_true));
     colnames(sarprobit.post.stat)<-c("post_mean","sd","q.025","q.50","q.975","truth");
     print(sarprobit.post.stat,digits=4)
     
     #Call gc() function to force garbage collector to return memory to the system
     #Use invisible function to suppress output of gc()
     invisible(gc());
}




























