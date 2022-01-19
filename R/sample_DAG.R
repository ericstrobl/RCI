sample_DAG <- function(nsamps, DAG){
  
  G = DAG$graph
  r = nrow(G)

  # err=matrix(rt(nsamps*r,df=5),nsamps,r)
  # data=matrix(rt(nsamps*r,df=5),nsamps,r)
  
  err=matrix(rnorm(nsamps*r),nsamps,r)
  data=matrix(rnorm(nsamps*r),nsamps,r)
  
  done=which(colSums(G)==0)
  stop=0;
  while (stop==0){
    for (s in done){
      ch=which(G[s,]==1)
      for (c in ch){
        pa=which(G[,c]==1)
        pa=sort(pa)
  
        h=intersect(pa,done)
  
        data[,c]=(data[,h,drop=FALSE]%*%DAG$weights[h,c])+err[,c]
        done=unique(c(done, c))
      }
    }
  
    if (length(done)==r){
      stop=1;
    }
  }
  
 return(data)
}

