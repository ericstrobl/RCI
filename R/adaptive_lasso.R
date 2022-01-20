adaptive_lasso <- function(X,B_ols){
  
  c = ncol(X)
  Bp = matrix(0,c,c)
  for (i in 2:c){        
      Yp = X[,i] 
      Xp = X[,1:(i-1),drop=FALSE]
      Bp[1:(i-1),i] = adalasso(Xp,Yp,B_ols[1:(i-1),i]) #regress on parents
  }
  
  return(Bp)
  
}

adalasso <- function(X,Y,B_ols){
  
    w = 1 / (abs(B_ols)^1)
    if (length(w)>1){
      X = X %*% diag(1/w)
      s = 2
    } else{
      X = X * (1/w)
      X = cbind(999,X)
      s = 3
    }
    cvfit = cv.glmnet(X,Y,intercept=FALSE)
    B_lasso = c(as.matrix(coef(cvfit, s = "lambda.min")))
    return(B_lasso[s:length(B_lasso)])
}


