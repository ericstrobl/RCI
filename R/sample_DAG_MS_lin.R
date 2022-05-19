sample_DAG_MS_lin <- function(X,Y,G,E,outm=NULL){
  
  # new X to be filled in
  Xn = matrix(0,nrow(X),ncol(X)+1)
  Xn[,1:ncol(X)] = E[,1:ncol(X)]
  
  if (is.null(outm)){
    outm = vector(mode = "list", length = (ncol(X)+1))
    mstart = TRUE
  } else{
    mstart = FALSE
  }
  
  means = Xn
  done=which(colSums(G)==0)
  stop=0;
  while (stop==0){
    for (s in done){
      ch=which(G[s,]==1) 
      for (c in ch){
        if (c %in% done){
          next
        }
        pa=which(G[,c]==1) 
        
        h=intersect(pa,done)
        if (setequal(h,pa)){ # if all parents already done
          if (length(h)>0){
            if (c != (ncol(X)+1)){
              if (mstart){
                outm[[c]] = PredictOut_lin(X[,h,drop=FALSE],X[,c,drop=FALSE],Xn[,h,drop=FALSE])
                Mean = outm[[c]]$Mean
              } else{
                Mean = Xn[,h,drop=FALSE] %*% outm[[c]]$modMean
              }
              Xn[,c] = E[,c,drop=FALSE] + Mean
            } else{
              if (mstart){
                outm[[c]] = PredictOut_lin(X[,h,drop=FALSE],as.matrix(Y),Xn[,h,drop=FALSE],binary_tar = TRUE)
                Mean = outm[[c]]$Mean
              } else{
                Mean = cbind(1,Xn[,h,drop=FALSE]) %*% outm[[c]]$modMean # because log odds
                
              }
              Xn[,c] = Mean
            }
          } else{
            if (c != (ncol(X)+1)){
              Xn[,c] = E[,c]
            } else{
              Xn[,c] = Y
            }
          }
          
          done=unique(c(done, c))
        }
      }
    }
    
    if (length(done) == (ncol(X)+1)){
      stop=1;
    }
  }
  
  
  return(list(X = Xn,E=E,outm=outm))
}

PredictOut_lin <- function(X,Y,X_test,binary_tar = FALSE){
  
  Y = as.matrix(Y)
  X = as.matrix(X)
  X_test = as.matrix(X_test)
  
  if (binary_tar){
    mod = lm.fit(cbind(1,X),c(Y))
    Res = cbind(1,X_test) %*% as.matrix(mod$coefficients)
    
    return(list(Mean = Res, modMean = as.matrix(mod$coefficients)))
    
  } else{
    ## MEAN
    mod = lm.fit(X,Y)
    Res = X_test %*% as.matrix(mod$coefficients)
    
    return(list(Mean = Res, modMean = as.matrix(mod$coefficients)))
    
  }
}
