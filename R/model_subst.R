model_subst <- function(X,Y,listL){
  
  # listL = DirectLiNGAM_regress(X) ## do this in cond_outl instead
  B = listL$B
  K = listL$K
  X = listL$X
  E = cbind(listL$E,0) # error term for D is zero
  B = cbind(B,adaptive_log_lasso(X,Y))
  B = rbind(B,0)
  G = (B!=0)
  
  # plot(as(G,"graphNEL"))
  
  i0 = which(Y==0)
  i1 = which(Y==1)
  n0 = length(i0)
  E0 = E[i0,]
  
  Anc = c()
  score = c()
  c = ncol(X)
  L0 = mean(sample_DAG_LogOdds(B,E0)[,c+1]) # score for healthy
  for (i in 1:c){ # compute marginal contributions
    if (isAnc(G,i,c+1)){
      Anc = c(Anc, i)

      En = E0
      En[,i]= sample(E[i1,i],n0,replace=TRUE) # bootstrap substitute the errors for diseased
      
      Ln = sample_DAG_LogOdds(B,En) # sample new DAG
      
      score = c(score,mean(Ln[,c+1]) - L0)
    }
  }
  
  return(list(order=Anc,scores=score)) # scores aligned with Anc
  
  
}
