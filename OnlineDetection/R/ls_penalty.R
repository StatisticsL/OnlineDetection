#' @title Solve the penalty optimization for linear model
#'
#' @description Solve the LASSO or Adaptive LASSO for linear model.
#' @param X The covariate in the model.
#' @param y The response variable in the model.
#' @param penalty.type The penalty type, LASSO or Adaptive LASSO. The default value is AdapLASSO.
#' @param gamma The factor in the Adaptive LASSO. The default value is 1.

#' @export
#' @return
#' \item{beta}{The solution of optimization with penalty.}

ls_penalty<-function(X,y,penalty.type="AdapLASSO",gamma=1){
  cv.init<-gcdnet::cv.gcdnet(X,y,method = 'ls',intercept=FALSE)
  beta<-gcdnet::gcdnet(X,y,method = 'ls',lambda=cv.init$lambda.min,intercept = FALSE)$beta
  beta<-beta*(abs(beta)>=cv.init$lambda.min)
  beta<-as.vector(beta)
  if(penalty.type=="LASSO"){
    return(beta)
  }else{
    weight=1/(pmax(abs(beta),1e-6)^(gamma))
    Xs=scale(X,center = FALSE,scale = weight)
    beta<-gcdnet::gcdnet(Xs,y,method = 'ls',lambda=cv.init$lambda.min,intercept = FALSE)$beta
    beta<-beta*(abs(beta)>=cv.init$lambda.min)
    beta<-as.vector(beta)/weight
    return(beta)
  }
}
