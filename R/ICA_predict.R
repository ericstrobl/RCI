ICA_predict <- function(X,Y){
  c = ncol(X)
  time = proc.time()
  mm = icafast(X, c)
  time = (proc.time() - time)[3]
  
  id = HungarianSolver(t(1/abs(mm$M)))$pairs[,2]
  E = mm$S[,order(id)]
  
  # mm$W = t(mm$W)
  # id = c()
  # for (i in 1:c){
  #   Mi = abs(mm$W[,i]) # what is the variable in X that contributes most to source S[,i]
  #   id  = c(id,which(Mi == max(Mi))[1])
  # }
  # E = mm$S[,order(id)]
  
  imp = randomForest(E, y=factor(Y), localImp=TRUE)$localImp
  # imp = randomForest(E, y=factor(Y), importance=TRUE)$importance
  # print(imp)
  
  return(list(scores = t(imp), order = 1:c, E = E, time = time))
  
}


