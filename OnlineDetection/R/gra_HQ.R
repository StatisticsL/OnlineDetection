#' @title The gradient of HQ optimization
#'
#' @description The gradient of HQ optimization.
#' @param H The matrix in the optimization.
#' @param Q The vector in the optimization.
#' @param beta The variable in the optimization.

#' @export
#' @return
#' \item{gra}{The gradient value.}

gra_HQ=function(H,Q,beta){
  gra=2*H%*%beta+Q
  return(gra)
}
