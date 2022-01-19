sample_DAG_LogOdds <- function(B,En){
  
  # sample from marginals
  Xn = En
  for (i in 1:ncol(B)){
    pa = which(B[,i] != 0)
    if (length(pa)>0){
      # print(pa)
      # print(ncol(En))
      Xn[,i] = Xn[,pa,drop=FALSE] %*% B[pa,i,drop=FALSE] + En[,i,drop=FALSE]
    } else{
      Xn[,i] = En[,i]
    }
  }
  
  return(Xn)
}