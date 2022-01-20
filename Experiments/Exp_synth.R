ns = c(100,1000,10000)
ps = c(10,50,100)
np = length(ps)
reps = 100

Gs = lapply(1:reps, function(.) lapply(1:length(ns),function(.) vector("list",np)))
RCI_res = Gs
RCI_noY_res = Gs
ICA_res = Gs
CO_res = Gs
MS_res = Gs
HPC_res = Gs
TT_res = Gs
LR_res = Gs

for (i in 1:reps){
  print(i)
  for (n in 1:length(ns)){
    for (p in 1:length(ps)){
      
      ### GENERATE DAG AND SAMPLE FROM IT
      G = generate_DAG(ps[p],en=2)
      Gs[[i]][[n]][[p]] = G
      X = sample_DAG_Y(ns[n],G)
      # print(X$Y)
      # plot(as(G$graph,"graphNEL"))
      
      nW = apply(X$E,2,sd)
      X$data[,-X$Y] = normalizeData(X$data[,-X$Y])
      Gs[[i]][[n]][[p]]$Y = X$Y
      
      
      ### COMPUTE GROUND TRUTH
      weights_norm = G$weights
      for (j1 in 1:ps[p]){
        weights_norm[j1,] = G$weights[j1,]*nW[j1] #first step is from error terms
      }
      
      beta_star = rep(0,ps[p])
      betas1 = weights_norm
      for (j in 1:(ps[p]-1)){
        beta_star = beta_star + betas1[,X$Y]
        betas1 = betas1 %*% G$weights
      }
      
      Gs[[i]][[n]][[p]]$betas = beta_star #beta_star with standardization
      
      Et = sweep(X$E[,-X$Y], 2, apply(X$E[,-X$Y],2,sd), FUN = '/')
      truth = Et %*% diag(beta_star[-X$Y])
      
      ### RUN ALGORITHMS
      ptm <- proc.time()
      out = RCI(X$data[,-X$Y],X$data[,X$Y])
      RCI_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3]
      RCI_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      RCI_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)
      RCI_res[[i]][[n]][[p]]$time_LiNGAM = out$time

      ptm <- proc.time()
      out = RCI_noY(X$data[,-X$Y],X$data[,X$Y])
      RCI_noY_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3]
      RCI_noY_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      RCI_noY_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)
      RCI_noY_res[[i]][[n]][[p]]$time_LiNGAM_fast = out$time

      out = LiNGAM_times(X$data[,-X$Y],X$data[,X$Y])
      RCI_res[[i]][[n]][[p]]$time_LiNGAM_slow_Y = out$time1
      # RCI_res[[i]][[n]][[p]]$time_LiNGAM_fast = out$time2
      
      ptm <- proc.time()
      out = ICA_predict(X$data[,-X$Y],X$data[,X$Y])
      ICA_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3]
      ICA_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      beta = glm.fit(cbind(out$E,1),X$data[,X$Y],family=binomial())$coefficients[1:ncol(out$E)] # compute scores using LR
      if (length(beta)==1){
        out$scores = out$E * beta
      } else{
        out$scores = out$E %*% diag(beta)
      }
      ICA_res[[i]][[n]][[p]]$rank_overlap_LR = eval_scores(truth,out)
      ICA_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)
      ICA_res[[i]][[n]][[p]]$time_ICA = out$time
      
      ptm <- proc.time()
      out = cond_outl(X$data[,-X$Y],X$data[,X$Y])
      CO_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3]
      CO_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      CO_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)
      CO_res[[i]][[n]][[p]]$time_ICA = out$listL$time

      MS_res[[i]][[n]][[p]]$time = out$listL$time
      ptm <- proc.time()
      out = model_subst(X$data[,-X$Y],X$data[,X$Y],out$listL)
      MS_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3] + MS_res[[i]][[n]][[p]]$time
      MS_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      MS_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)

      ptm <- proc.time()
      out = HITON_PC(X$data[,-X$Y],X$data[,X$Y])
      HPC_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3]
      HPC_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      HPC_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)

      ptm <- proc.time()
      out = ttest_scores(X$data[,-X$Y],X$data[,X$Y])
      TT_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3]
      TT_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      TT_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)

      ptm <- proc.time()
      out = Adaptive_LR(X$data[,-X$Y],X$data[,X$Y])
      LR_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3]
      LR_res[[i]][[n]][[p]]$rank_overlap = eval_scores(truth,out)
      LR_res[[i]][[n]][[p]]$L2 = eval_L2(truth,out)


      ## SAVE RESULTS
      save(file="Results_synthetic.RData",Gs,RCI_res,RCI_noY_res,ICA_res,
           CO_res,MS_res,HPC_res,TT_res,LR_res)
      
    }
  }
}

####
## RBO RESULTS

RBO_RCI = matrix(0,length(ns),length(ps))
RBO_noY_RCI = matrix(0,length(ns),length(ps))
RBO_ICA = RBO_RCI
RBO_CO = RBO_RCI
RBO_MS = RBO_RCI
RBO_HPC = RBO_RCI
RBO_TT = RBO_RCI
RBO_LR = RBO_RCI

for (i in 1:reps){
  for (n in 1:length(ns)){
    for (p in 1:length(ps)){
      
      RBO_RCI[n,p] = RBO_RCI[n,p] + RCI_res[[i]][[n]][[p]]$rank_overlap
      RBO_noY_RCI[n,p] = RBO_noY_RCI[n,p] + RCI_noY_res[[i]][[n]][[p]]$rank_overlap
      RBO_ICA[n,p] = RBO_ICA[n,p] + ICA_res[[i]][[n]][[p]]$rank_overlap
      RBO_CO[n,p] = RBO_CO[n,p] + CO_res[[i]][[n]][[p]]$rank_overlap
      RBO_MS[n,p] = RBO_MS[n,p] + MS_res[[i]][[n]][[p]]$rank_overlap
      RBO_HPC[n,p] = RBO_HPC[n,p] + HPC_res[[i]][[n]][[p]]$rank_overlap
      RBO_TT[n,p] = RBO_TT[n,p] + TT_res[[i]][[n]][[p]]$rank_overlap
      RBO_LR[n,p] = RBO_LR[n,p] + LR_res[[i]][[n]][[p]]$rank_overlap
      
      
    }
  }
}

print(RBO_RCI/reps)
print(RBO_noY_RCI/reps)
print(RBO_ICA/reps)
print(RBO_CO/reps)
print(RBO_MS/reps)
print(RBO_HPC/reps)
print(RBO_TT/reps)
print(RBO_LR/reps)


######
## MSE RESULTS

MSE_RCI = matrix(0,length(ns),length(ps))
MSE_noY_RCI = MSE_RCI
MSE_ICA = MSE_RCI
MSE_CO = MSE_RCI
MSE_MS = MSE_RCI
MSE_HPC = MSE_RCI
MSE_TT = MSE_RCI
MSE_LR = MSE_RCI

for (i in 1:100){
  for (n in 1:length(ns)){
    for (p in 1:length(ps)){
      
      MSE_RCI[n,p] = MSE_RCI[n,p] + RCI_res[[i]][[n]][[p]]$L2
      MSE_noY_RCI[n,p] = MSE_noY_RCI[n,p] + RCI_noY_res[[i]][[n]][[p]]$L2
      MSE_ICA[n,p] = MSE_ICA[n,p] + ICA_res[[i]][[n]][[p]]$L2
      MSE_CO[n,p] = MSE_CO[n,p] + CO_res[[i]][[n]][[p]]$L2
      MSE_MS[n,p] = MSE_MS[n,p] + MS_res[[i]][[n]][[p]]$L2
      MSE_HPC[n,p] = MSE_HPC[n,p] + HPC_res[[i]][[n]][[p]]$L2
      MSE_TT[n,p] = MSE_TT[n,p] + TT_res[[i]][[n]][[p]]$L2
      MSE_LR[n,p] = MSE_LR[n,p] + LR_res[[i]][[n]][[p]]$L2
      
      
    }
  }
}

print(MSE_RCI/reps)
print(rowMeans(MSE_RCI/reps))
print(rowMeans(MSE_noY_RCI/reps))
print(MSE_ICA/reps)
print(MSE_CO/reps)
print(MSE_MS/reps)
print(MSE_HPC/reps)
print(MSE_TT/reps)
print(MSE_LR/reps)


######
## TIMING RESULTS

Time_local_plus = matrix(0,length(ns),length(ps))
Time_plus = Time_local_plus
Time_local = Time_local_plus
Time_original = Time_local_plus

for (i in 1:100){
  for (n in 1:length(ns)){
    for (p in 1:length(ps)){
      
      Time_local_plus[n,p] = Time_local_plus[n,p] + RCI_res[[i]][[n]][[p]]$time_LiNGAM
      Time_plus[n,p] = Time_plus[n,p] + RCI_noY_res[[i]][[n]][[p]]$time_LiNGAM_fast
      Time_local[n,p] = Time_local[n,p] + RCI_res[[i]][[n]][[p]]$time_LiNGAM_slow_Y
      Time_original[n,p] = Time_original[n,p] + CO_res[[i]][[n]][[p]]$time_ICA
      
      
    }
  }
}

print(Time_local_plus/reps)
print(Time_local/reps)
print(Time_plus/reps)
print(Time_original/reps)


