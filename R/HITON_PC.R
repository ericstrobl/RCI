HITON_PC <- function(X,Y,alpha=0.01) 
{
  # checked 1/1/22
  XY = cbind(X,Y)
  
  pc = learn.nbr(as.data.frame(XY),'Y',method='si.hiton.pc',alpha=alpha)
  pc = as.numeric(gsub('V', '', pc))
  
  if (length(pc)==0){
    return(list(order=pc,scores=NULL))
  }
  
  beta = glm.fit(cbind(X[,pc,drop=FALSE],1),Y,family=binomial())$coefficients[seq_len(length(pc))]
  if (length(beta)==1){
    delta = X[,pc,drop=FALSE] * beta
  } else{
    delta = X[,pc,drop=FALSE] %*% diag(beta)
  }
  
  return(list(order=pc,scores=delta))
}