##############################################################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
packages <- c("MASS","OnlineDetection")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))



##############################################################################
###parameters setting-----------------------------------
p=20
B=50##the number of batches--------
K=10##the number of batches for one group
nb=200##sample size of every batches
penalty.type="AdapLASSO"
Sigma=diag(p)

Nrep=20

beta_group=matrix(0,p,5)
beta_group[1:5,1]=c(0.5,0.8,0.7,1,2)
beta_group[1:5,2]=c(0.6,0.8,0.7,0.8,2)
beta_group[1:5,3]=c(2.5,0.8,1,1.5,3)
beta_group[1:5,4]=c(3,1,2,1,3)
beta_group[1:5,5]=c(3,0.9,1.9,1,3)

Beta_True=NULL
Beta_Target=NULL
Beta_Oracle=NULL
Beta_Tran=NULL

for (repli in 1:Nrep) {
  print(paste("At ",repli,"-th iteration"))
  ###save for one time-----------------
  betatran=matrix(0,B,p)
  betatrue=matrix(0,B,p)
  betatarget=matrix(0,B,p)
  betaoracle=matrix(0,B,p)

  Xora=NULL
  yora=NULL

  para_index=1

  for (b in 1:B) {
    para_index_origal=para_index
    para_index=(b-1)%/%K+1
    betatrue[b,]=beta_group[,para_index]
    Xt=mvrnorm(nb,mu=rep(0,p),Sigma = Sigma)
    yt=Xt%*%betatrue[b,]+rnorm(nb)
    nt=length(yt)

    ##high-dimension
    betatarget[b,]=ls_penalty(Xt,yt,penalty.type = penalty.type)

    ###Oracle
    if(para_index_origal!=para_index){Xora=NULL;yora=NULL}
    Xora=rbind(Xora,Xt)
    yora=rbind(yora,yt)
    nora=length(yora)

    betaoracle[b,]=ls_penalty(Xora,yora,penalty.type = penalty.type)


    ###Trans--------------------------------------
    if(b==1){
      ##initial estimate--------------------------
      betat_adalasso=betatarget[b,]
      ###main function------------------------------
      beta_s=betat_adalasso
      XsXs=t(Xt)%*%Xt
      Xsys=t(Xt)%*%yt
      ysys=sum(yt^2)
      ns=dim(Xt)[1]
      betatran[b,]=beta_s
    }else{
      res_tran<-onlineTL_ls(XsXs,Xsys,ysys,ns,beta_s,Xt,yt)
      betatran[b,]=res_tran$beta
      ysys=res_tran$ysys
      Xsys=res_tran$Xsys
      XsXs=res_tran$XsXs
      ns=res_tran$ns
      beta_s=res_tran$beta
    }
    ############################################
    print(paste("     ","At ",b,"-th batch"))
  }
  Beta_True[[repli]]=betatrue
  Beta_Target[[repli]]=betatarget
  Beta_Oracle[[repli]]=betaoracle
  Beta_Tran[[repli]]=betatran
}

###analysis##############################################################
###MSE----------------------------------------------
mse_rep_target=matrix(0,Nrep,B)
mse_rep_oracle=matrix(0,Nrep,B)
mse_rep_tran=matrix(0,Nrep,B)

for (repli in 1:Nrep) {
  mse_rep_target[repli,]=apply((Beta_Target[[repli]]-Beta_True[[repli]])^2,1,mean)
  mse_rep_oracle[repli,]=apply((Beta_Oracle[[repli]]-Beta_True[[repli]])^2,1,mean)
  mse_rep_tran[repli,]=apply((Beta_Tran[[repli]]-Beta_True[[repli]])^2,1,mean)
}

mse=matrix(0,B,3)
mse[,1]=apply(mse_rep_target,2,mean)
mse[,2]=apply(mse_rep_tran,2,mean)
mse[,3]=apply(mse_rep_oracle,2,mean)
colnames(mse)=c("Target","Tran","Oracle")


###---plots-------------------------------------
mse_plot=mse*(10^3)
plot(1:B,mse_plot[,"Target"],type="l",ylim=c(0,1.2*max(mse_plot)),
     ylab=expression(MSE %*% 10^3),xlab = "Blocks",col="2",lwd=2)
lines(1:B,mse_plot[,"Tran"],type = "l",col="3",lwd=2)
lines(1:B,mse_plot[,"Oracle"],type = "l",col="4",lwd=2)
legend("topleft",legend = c("Sin","Ours","Ora"),col = c(2,3,4),lty=1,lwd=2)

