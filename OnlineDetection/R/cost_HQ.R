#' @title The cost of HQ optimization
#'
#' @description The cost function of HQ optimization.
#' @param H The matrix in the optimization.
#' @param Q The vector in the optimization.
#' @param beta The variable in the optimization.
#' @param lambda The penalty factor in the optimization.

#' @export
#' @return
#' \item{cost}{The cost value.}

cost_HQ<-function(H,Q,beta,lambda){
  cost=t(beta)%*%H%*%beta+t(Q)%*%beta+sum(abs(lambda*beta))
  return(cost)
}
