eval_L2 <- function(true_scores,output){
  
  est_scores = matrix(0,nrow(true_scores),ncol(true_scores))
  est_scores[,output$order]= output$scores
  
  return(mean( (true_scores - est_scores)^2 ))
}