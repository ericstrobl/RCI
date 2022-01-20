RCI_noY <- function(X,Y,alpha=0.2){
  # checked 1/1/22
  
  time = proc.time()
  outL = DirectLiNGAM_fast(X)
  # outL = DirectLiNGAM(X)
  time = (proc.time() - time)[3]
  K = outL$K;
  if (length(K)==0){
    delta = matrix(0,nrow(X),ncol(X))
    return(list(delta=delta,K=K))
  }
  E = outL$X[,order(K),drop=FALSE]
  K = sort(K)
  
  ## patient-specific root causes
  beta = glm.fit(cbind(E,1),Y,family=binomial())$coefficients[1:length(K)]
  if (length(beta)==1){
    delta = E * beta
  } else{
    delta = E %*% diag(beta)
  }
  
  return(list(scores=delta,order=K,E=E, time = time))
}