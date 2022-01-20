OLS_LiNGAM <- function(X,Cov){
  
  listQL = QL(inverse.sqrt(Cov));
  L = listQL$L
  
  c = ncol(X)
  L = L/(as.matrix(diag(L)) %*% matrix(1,1,c))
  
  B = diag(c)-L
  
  return(t(B))
  
}

QR <- function(A){
  tmp <- qr(A)
  Q <- qr.Q(tmp)
  R <- qr.R(tmp)
  sign_flips <- sign(diag(R))
  sign_flips[sign_flips == 0] <- 1
  Q <- Q %*% diag(sign_flips)
  R <- diag(sign_flips) %*% R
  return(list(Q=Q, R=R))
}

QL <- function(A){
  B <- A[,ncol(A):1]
  tmp <- QR(B)
  Q <- tmp$Q
  R <- tmp$R
  Q <- Q[,ncol(Q):1]
  L <- R[nrow(R):1,ncol(R):1]
  return(list(Q=Q, L=L))
}

inverse.sqrt <- function(S){
  ei <- eigen(S)
  V <- ei$vectors
  id = which(ei$values>1E-10)
  res <- V[,id] %*% diag(1 / sqrt(ei$values[id])) %*% t(V[,id])
  return(res)
}