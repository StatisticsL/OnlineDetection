#' @title Thresholding function of LASSO
#'
#' @description Thresholding function of LASSO.
#' @param z The constant in the loss .
#' @param lambda The penalty parameter.

#' @export
#' @return
#' \item{fz}{The solution.}
thresh_fct=function(z,lambda){
  fz=sign(z)*pmax(abs(z)-lambda,0)
  return(fz)
}
