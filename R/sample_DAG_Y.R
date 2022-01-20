sample_DAG_Y <- function(nsamps, DAG){
  
  G = DAG$graph
  r = nrow(G)
  
  Y = sample(which(rowSums(G)==0 & colSums(G)>0),1) #terminal vertex with at least one parent
  
  err=matrix(0,nsamps,r)
  ix = sample(1:3,r,replace=TRUE)
  # ix = sample(1,r,replace=TRUE)
  for (i in 1:r){
    if (ix[i] == 1){
      err[,i]=matrix(rt(nsamps,df=5),nsamps)
    } else if (ix[i]==2){
      err[,i]=matrix(runif(nsamps,-1,1),nsamps)
    } else if (ix[i]==3){
      err[,i]=matrix(rchisq(nsamps,df=3)-3,nsamps)
    }
  }
  err[,Y]=0 #error for diagnosis is zero
  data=err
  
  done=which(colSums(G)==0) # which variables do not have parents
  stop=FALSE;
  while (!stop){
    for (s in done){
      ch=setdiff(which(G[s,]==1),done)
      for (c in ch){
        pa=which(G[,c]==1)
        
        if (all(pa %in% done)){
          # print("parents")
          # print(pa)
          # print("children")
          # print(c)
          data[,c]=(data[,pa,drop=FALSE]%*%DAG$weights[pa,c])+err[,c]
          done=unique(c(done, c))
          break
        }
      }
    }
    
    if (length(done)==r){
      stop=TRUE;
    }
  }
  
  Y0 = data[,Y]
  pY = logistic(data[,Y])
  for (i in 1:nsamps){
    data[i,Y] = rbinom(n=1,size=1,prob=pY[i]) 
  }
  
  return(list(data=data,Y=Y,E=err,Y0=Y0))
}


logistic <- function(X){
  
  return(1/(1+exp(-X)))
}