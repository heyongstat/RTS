#----------------------------------------------------------------------------------------
#  Spatial Multivariate Kendall' tau Matrix
#  Input:
#           x ------ n x p  data matrix (n: time dimension, p: cross-section)
#  Output:
#          TK ------ Sample Spatial Multivariate Kendall' tau Matrix
#----------------------------------------------------------------------------------------

SK<-function(X){
  p <- ncol(X)
  n <- nrow(X)
  TK <- matrix(0,p,p)
  for(i in 1:(n-2)){
    TT <- matrix(rep(X[i,],n-i),n-i,p,byrow = TRUE)-X[(i+1):n,]
    TT <- t(diag(1/diag(TT%*%t(TT)))%*%TT)%*%TT
    TK <- TK+TT
  }
  TT <- X[n-1,]-X[n,]
  TK <- TK+TT%*%t(TT)/sum(TT^2)
  TK <- 2/(n*(n-1))*TK
  return(TK)
}

#----------------------------------------------------------------------------------------
#  RTS Method
#  Input:
#           x ------ n x p  data matrix (n: time dimension, p: cross-section)
#           r ------ number of factors
#  Output:
#        Fhat ------ n x r  the estimated factor scores 
#        Lhat ------ p x r  the estimated factor loadings 
#----------------------------------------------------------------------------------------

RTS <- function(X,r){
    p <- ncol(X)
    n <- nrow(X)
    Khat <- SK(X)
    Lhat <- sqrt(p)*as.matrix(eigen(Khat)$vectors[,1:r]) 
    Fhat <- matrix(0,n,r)
    for (i in 1:n){
	  Fhat[i,]=lm(X[i,]~Lhat-1)$coefficients
    }
    return(list(Fhat=Fhat,Lhat=Lhat))
}

#----------------------------------------------------------------------------------------
#  PCA Method
#  Input:
#           x ------ n x p  data matrix (n: time dimension, p: cross-section)
#           r ------ number of factors
#  Output:
#        Fhat ------ n x r  the estimated factor scores 
#        Lhat ------ p x r  the estimated factor loadings 
#----------------------------------------------------------------------------------------


PCA <- function(X,r){
   p <- ncol(X)
   n <- nrow(X)
   Fhat <- sqrt(n)*as.matrix(eigen(X%*%t(X)/n)$vectors[,1:r])
   Lhat <- t(X)%*%Fhat/n
   return(list(Fhat=Fhat,Lhat=Lhat))
}
