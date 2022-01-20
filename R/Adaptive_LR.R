Adaptive_LR <-function(X,Y){
  
  B = adaptive_log_lasso1(X,Y)
  
  return(list(scores = X %*% diag(B), order = 1:ncol(X)))
}

adaptive_log_lasso1 <- function(X,Y){
  
  B_ols = glm.fit(cbind(X,1),Y,family=binomial())$coefficients[1:ncol(X)]
  B = ada_log_lasso1(X,Y,B_ols);
  
  return(B)
}

ada_log_lasso1 <- function(X,Y,B_ols){
  
  w = 1 / (abs(B_ols)^1)
  if (length(w)>1){
    X = X %*% diag(1/w)
    s = 2
  } else{
    X = X * (1/w)
    X = cbind(999,X)
    s = 3
  }
  cvfit = cv.glmnet(X,Y,family="binomial",type.measure = "class")
  B_lasso = c(as.matrix(coef(cvfit, s = "lambda.min")))
  return(B_lasso[s:length(B_lasso)])
}