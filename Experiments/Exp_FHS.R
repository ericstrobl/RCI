data = read.csv('FHS.csv') #requires prior authorization from https://biolincc.nhlbi.nih.gov/studies/framcohort/
data = data[,-c(5,13,21)]
data = as.matrix(data[complete.cases(data),])
BMI = 25
Y = 0+(data[,5]>=BMI)
X = data[,-c(5,12,19,9,16)]

Xo = X
X = normalizeData(X)

reps = 1000

RCI_res = vector("list",reps)
RCI_noY_res = RCI_res
ICA_res = RCI_res
CO_res = RCI_res
MS_res = RCI_res
HPC_res = RCI_res
TT_res = RCI_res
LR_res = RCI_res

r= length(Y)
c = ncol(X)
for (i in 1:reps){
  print(i)
  
  is = sample(r,r,replace=TRUE)
  Xt = X[is,]
  Yt = Y[is]
  
  id = sample(1:c,c,replace=FALSE)
  Xt = Xt[,id]
  
  ### RUN ALGORITHMS
  ptm <- proc.time()
  out = RCI(Xt,Yt)
  RCI_res[[i]]$time = (proc.time() - ptm)[3]
  RCI_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)
  RCI_res[[i]]$time_LiNGAM = out$time

  ptm <- proc.time()
  out = RCI_noY(Xt,Yt)
  RCI_noY_res[[i]]$time = (proc.time() - ptm)[3]
  RCI_noY_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)
  RCI_noY_res[[i]]$time_LiNGAM_fast = out$time

  out = LiNGAM_times(Xt,Yt)
  RCI_res[[i]]$time_LiNGAM_slow_Y = out$time1
  # RCI_res[[i]]$time_LiNGAM_fast = out$time2
  
  ptm <- proc.time()
  out = ICA_predict(Xt,Yt)
  ICA_res[[i]]$time = (proc.time() - ptm)[3]
  ICA_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)
  beta = glm.fit(cbind(out$E,1),Yt,family=binomial())$coefficients[1:ncol(out$E)] # compute scores using LR
  if (length(beta)==1){
    out$scores = out$E * beta
  } else{
    out$scores = out$E %*% diag(beta)
  }
  ICA_res[[i]]$rank_overlap_LR = eval_scores_fhs(out,X[is,],id)
  ICA_res[[i]]$time_ICA = out$time
  
  ptm <- proc.time()
  out = cond_outl(Xt,Yt)
  CO_res[[i]]$time = (proc.time() - ptm)[3]
  CO_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)
  CO_res[[i]]$time_ICA = out$listL$time

  MS_res[[i]]$time = out$listL$time
  ptm <- proc.time()
  out = model_subst(Xt,Yt,out$listL)
  MS_res[[i]]$time = (proc.time() - ptm)[3] + MS_res[[i]]$time
  MS_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)

  ptm <- proc.time()
  Xm = as.matrix(Xt); colnames(Xm) = NULL;
  out = HITON_PC(Xm,Yt)
  HPC_res[[i]]$time = (proc.time() - ptm)[3]
  HPC_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)

  ptm <- proc.time()
  out = ttest_scores(Xt,Yt)
  TT_res[[i]]$time = (proc.time() - ptm)[3]
  TT_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)

  ptm <- proc.time()
  out = Adaptive_LR(Xt,Yt)
  LR_res[[i]]$time = (proc.time() - ptm)[3]
  LR_res[[i]]$rank_overlap = eval_scores_fhs(out,X[is,],id)
  
  
  ### SAVE RESULTS
  save(file="Res_FHS.RData",Gs,RCI_res,RCI_noY_res,ICA_res,
       CO_res,MS_res,HPC_res,TT_res,LR_res)
  
}


### RBO
RBO = matrix(0,reps,8)
for (i in 1:reps){
  RBO[i,1] = RCI_res[[i]]$rank_overlap
  RBO[i,2] = RCI_noY_res[[i]]$rank_overlap
  RBO[i,3] = ICA_res[[i]]$rank_overlap
  RBO[i,4] = CO_res[[i]]$rank_overlap
  RBO[i,5] = MS_res[[i]]$rank_overlap
  RBO[i,6] = HPC_res[[i]]$rank_overlap
  RBO[i,7] = TT_res[[i]]$rank_overlap
  RBO[i,8] = LR_res[[i]]$rank_overlap
}

print(colMeans(RBO))


### TIMING
Time = matrix(0,reps,4)
for (i in 1:reps){
  Time[i,1] = RCI_res[[i]]$time_LiNGAM
  Time[i,2] = RCI_noY_res[[i]]$time_LiNGAM_fast
  Time[i,3] = RCI_res[[i]]$time_LiNGAM_slow_Y
  Time[i,4] = CO_res[[i]]$time_ICA
}

print(colMeans(Time[,4]/Time))
