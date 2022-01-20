ttest_scores <- function(X,Y){
  # checked 1/1/22
  
  X0 = X[Y==0,,drop=FALSE]
  X1 = X[Y==1,,drop=FALSE]
  
  ms = colMeans(X0) - colMeans(X1)
  
  n0 = nrow(X0)
  n1 = nrow(X1)
  
  var0 = apply(X0,2,var)
  var1 = apply(X1,2,var)

  t = ms/ sqrt(var0/n0 + var1/n1)
  
  delta = t(matrix(abs(t),ncol(X),nrow(X)))
  
  return( list(scores=delta,order=1:ncol(X)) )
  
  
  
}