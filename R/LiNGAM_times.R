LiNGAM_times <- function(X,Y,alpha=0.2){
  
  time1 = proc.time()
  outL = DirectLiNGAM_slow_Y(X,Y)
  time1 = (proc.time() - time1)[3]
  
  # time2 = proc.time()
  # outL = DirectLiNGAM_fast1(X)
  # time2 = (proc.time() - time2)[3]
  
  return(list(time1 = time1))
  
}