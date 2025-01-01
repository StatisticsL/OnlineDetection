#' @title Solve HQ optimization with ista algorithm
#'
#' @description Solve HQ optimization with ista algorithm.
#' @param H The matrix in the optimization.
#' @param Q The vector in the optimization.
#' @param lambda The penalty factor in the optimization.
#' @param beta_start The initial value for the variable.
#' @param conv_thresh The relative error.

#' @export
#' @return
#' \item {beta}{The solution for the HQ optimization.}
#' \item {costS}{The cost valueS at each step.}


ista_HQ<-function(H,Q,lambda,beta_start=NA,conv_thresh=1.e-8){
  ###adaptivelasso：lambdaseq
  p=length(Q)
  tt=0.1
  if(any(is.na(beta_start))){
    beta=rep(0,p)
  }else{
    beta=beta_start
  }

  cost=cost_HQ(H,Q,beta,lambda)
  costs=cost
  g=gra_HQ(H,Q,beta)

  continue=TRUE
  while(continue){
    z=beta-tt*g
    betanew=thresh_fct(z,lambda*tt)
    costnew=cost_HQ(H,Q,betanew,lambda)
    continue_step=(costnew>cost)
    while(continue_step){
      tt=tt*0.9
      z=beta-tt*g
      betanew=thresh_fct(z,lambda*tt)
      costnew=cost_HQ(H,Q,betanew,lambda)
      continue_step=(costnew>cost)
    }

    continue=((abs(costnew-cost)/abs(cost))>conv_thresh)

    costs=c(costs,costnew)
    cost=costnew
    beta=betanew
    g=gra_HQ(H,Q,beta)
  }

  out=NULL
  out$beta=beta
  out$costs=costs
  return(out)
}
