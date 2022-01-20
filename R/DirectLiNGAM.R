DirectLiNGAM <- function(X){
  X = as.matrix(X)
  
  K = c() #1 
  U = 1:ncol(X) #2
  X = normalizeData(X) #3
  Cov = cov(X) #4
  
  repeat{ #11
    
    root = FindRoot(X,U,Cov) #6
    
    K = c(K,root) #7
    U = U[-which(U==root)] #8
    
    X = UpdateData(X,U,Cov,root) #9
    if (length(U)==1){
      K = c(K,U)
      break
    }
    Cov = UpdateCovMat(U,Cov,root) #10
  }
  
  return(K) #output
}