generate_DAG <- function(p,en){
  
  N = p*p - p;
  samplesB = rbinom(N/2,1, en/(p-1) );
  
  DAG = list()
  graph = matrix(0,p,p)
  graph[upper.tri(graph, diag=FALSE)] <- samplesB;
  
  DAG = list()
  DAG$graph = graph
  
  weights = matrix((0.75*runif(p^2)+0.25)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  DAG$weights = weights*DAG$graph
  
  ord = sample(1:p,p,replace=FALSE) # permute order
  DAG$graph = DAG$graph[ord,ord]
  DAG$weights = DAG$weights[ord,ord]
  
  return(DAG)
}