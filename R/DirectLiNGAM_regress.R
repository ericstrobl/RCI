DirectLiNGAM_regress <- function(X){
  
  time = proc.time()
  
  K = c() #1 
  U = 1:ncol(X) #2
  X = normalizeData(X) #3
  Cov = cov(X) #4
  
  X0 = X
  Cov0 = Cov
  
  repeat{ #11
    
    root = FindRoot(X,U,Cov) #6
    
    K = c(K,root) #7
    U = U[-which(U==root)] #8
    
    if (length(U)==0){
      break
    }
    X = UpdateData(X,U,Cov,root) #9
    Cov = UpdateCovMat(U,Cov,root) #10
  }
  
  time = (proc.time()-time)[3]
  
  B_ols = OLS_LiNGAM(X0[,K],Cov0[K,K])
  B_lasso = adaptive_lasso(X0[,K],B_ols)
  oK = order(K)
  
  return(list(K = K[oK],B = B_lasso[oK,oK], E = X, X=X0, time = time)) #output
}