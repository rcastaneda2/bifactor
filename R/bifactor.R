
#' title: Square a number
#' 
#' Takes in any numeric value and squares it
#' @param N A number of observations per dataset to generate
#' @param items A number of items per dataset to generate
#' @param reps Number of replications
#' @return Rturns Correlations between thetas, and RMSEA values for EAP scores
#' @depends mirt, mvtnorm
#' @imports mirt, mvtnorm
#' @export 

#devtools::use_package("mirt","mvtnorm") # Defaults to imports


gen<-function(N,items,reps){
  factors<-3
  items=items
  for(i in 1:reps){
    #set.seed(reps)
    #Create the Variance Covariance Matrix for 6 items
    sigma<-diag(factors)
    
    theta=as.matrix(mvtnorm::rmvnorm(n=N,mean=c(rep(1,factors)),sigma=sigma))
    
    
    a1<-matrix(1.5,items)
    a2<-matrix(c(rep(1.5,items/2),rep(NA,items/2)))
    a3<-matrix(c(rep(NA,items/2),rep(1.5,items/2)))
    a<-cbind(a1,a2,a3)
    d<-matrix(rnorm(items))
    data<-mirt::simdata(a,d,N,itemtype='dich',Theta=theta,sigma=sigma)
    
    mod1 <- mirt::bfactor(data, c(rep(1,items/2),rep(2,items/2)), itemtype='2PL')
    coe<-coef(mod1)
    scores<-mirt::fscores(mod1,full.scores=TRUE)
    cor1<-cor(scores[,1],theta[,1])
    cor2<-cor(scores[,2],theta[,2])
    cor3<-cor(scores[,3],theta[,3])
    cat('\n')
    cat("rep ", reps,'\n')
    cat("correlation 1 is ",cor1,'\n')
    cat("correlation 2 is ",cor2,'\n')
    cat("correlation 3 is ",cor3,'\n')	
    
    
    #Calculate RMSE for one replication
    cat("RMSE theta 1 ",mean((theta[,1]-scores[,1])^2),'\n')
    cat("RMSE theta 2 ",mean((theta[,2]-scores[,2])^2),'\n')
    cat("RMSE theta 3 ",mean((theta[,3]-scores[,3])^2),'\n')
    
    
  }
  
}