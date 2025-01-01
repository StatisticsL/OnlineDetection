#' @title The updating for streaming batch data with TL
#'
#' @description The updating for streaming batch data with TL.
#' @param XsXs The statistic \eqn{X^{T}X}.
#' @param Xsys The statistic \eqn{X^{T}y}.
#' @param ysys The statistic \eqn{y^{T}y}.
#' @param ns The number of source data.
#' @param beta_s The estimate from the source data.
#' @param Xt The covariavte from the target batch.
#' @param yt The response from the target batch.
#' @param penalty.type The penalty type, LASSO or Adaptive LASSO. The default value is AdapLASSO.
#' @param gamma The factor in the Adaptive LASSO. The default value is 1.

#' @export
#' @return
#' \item{beta}{The solution of updating.}
#' \item{XsXs}{The updated statistic \eqn{X^{T}X}.}
#' \item{Xsys}{The updated statistic \eqn{X^{T}y}.}
#' \item{ysys}{The updated statistic \eqn{y^{T}y}.}
#' \item{ns}{The updated number of source data.}

onlineTL_ls<-function(XsXs,Xsys,ysys,ns,beta_s,Xt,yt,penalty.type="AdapLASSO",gamma=1){
  p=dim(Xt)[2]
  nt=dim(Xt)[1]

  bias_t=ls_penalty(Xt,yt-Xt%*%beta_s,penalty.type = penalty.type,gamma = gamma)

  ysys=ysys+2*t(bias_t)%*%Xsys+t(bias_t)%*%XsXs%*%bias_t
  Xsys=Xsys+XsXs%*%bias_t

  H=XsXs/(2*ns)
  Q=-Xsys/ns

  beta_init=ista_HQ(H,Q,lambda = sqrt(log(p)/ns))$beta
  lambdaseq_s=sqrt(log(p)/ns)/abs(beta_init)
  lambdaseq_s[which(lambdaseq_s==Inf)]=1e+06
  beta_s=ista_HQ(H,Q,lambda = lambdaseq_s)$beta


  bias_t=ls_penalty(Xt,yt-Xt%*%beta_s,penalty.type = "LASSO" )

  beta=beta_s+bias_t

  ysys=ysys+2*t(bias_t)%*%Xsys+t(bias_t)%*%XsXs%*%bias_t
  Xsys=Xsys+XsXs%*%bias_t

  XsXs=XsXs+t(Xt)%*%Xt
  Xsys=Xsys+t(Xt)%*%yt
  ysys=ysys+sum(yt^2)

  ns=ns+nt

  out=NULL
  out$beta=beta
  out$XsXs=XsXs
  out$Xsys=Xsys
  out$ysys=ysys
  out$ns=ns
  return(out)
}
