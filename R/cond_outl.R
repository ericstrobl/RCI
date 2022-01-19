cond_outl <- function(X,Y){
  
  listL = DirectLiNGAM_regress(X)
  B = listL$B
  K = listL$K
  X = listL$X
  B = cbind(B,adaptive_log_lasso(X,Y))
  B = rbind(B,0)
  G = (B!=0)
  
  Anc = c()
  c = ncol(X)
  CO = matrix(0,nrow(X),c)
  for (i in 1:c){
    if (isAnc(G,i,c+1)){
        Anc = c(Anc, i)
        
        pa = which(G[,i])
        num = (X[,i,drop=FALSE] - X[,pa,drop=FALSE] %*% B[pa,i,drop=FALSE])
        CO[,i] = abs(num)/sd(num)
    }
  }
  
  return(list(order=Anc,scores=CO[,Anc], listL = listL))
  
}

adaptive_log_lasso <- function(X,Y){
  
  B_ols = glm.fit(cbind(X,1),Y,family=binomial())$coefficients[1:ncol(X)]
  B = ada_log_lasso(X,Y,B_ols);
  
  return(B)
  
}

ada_log_lasso <- function(X,Y,B_ols){
  
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

isAnc <- function(graph,a,b,visited=rep(FALSE,nrow(graph)))
{
  
  if (a %in% b){
    return(TRUE)
  }
  
  visited[a] = TRUE;
  
  adj = which(graph[a,] & !visited);
  
  out=FALSE;
  for (j in adj){
    out=isAnc(graph, j, b, visited);
    if(out==TRUE){
      break;
    }
  }
  
  return(out)
  
}


