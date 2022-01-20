FindRoot <- function(X,U,Cov){
  
  r = length(U) #4
  
  if (r==1){ #1
    return(U) #2
  }
  
  p = ncol(Cov)
  # M = matrix(0,p,p)
  S = rep(0,r)
  for (i in 1:(r-1)){
    # V = U[U>i]
    for (j in (i+1):r){

      # M[i,j] = Compare2(X,i,j,Cov)
      # M[j,i] = -M[i,j]
      
      score = Compare2(X,U[i],U[j],Cov)
      S[i] = S[i] + min(0,score)^2
      S[j] = S[j] + min(0,-score)^2
    }
  }

  # S = rowSums(pmin(M[U,U],0)^2)
  
  # print(S)
  
  root = U[S==min(S)][1]
  
  return(root) #output
  
}