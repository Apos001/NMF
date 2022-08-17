
# Creates the S matrix for a clustering lambda, which is an n-dimensional numeric vector, and k is the number of clusters
# n is the number of objects
my_S = function(lambda,k){
  n = length(lambda)
  H = matrix(0, nrow=n, ncol=k)
  for(i in 1:n) H[i,lambda[i]]  = 1
  S = H %*% t(H)
  return(S)
}

# assume that you have at hand 3 clusterings, so you create all 3 S matrices
S_1 = my_S(lambda_1,k)
S_2 = my_S(lambda_2,k)
S_3 = my_S(lambda_3,k)

# then you create the combined S matrix
S_data = (S_1+S_2+S_3)/3

S = S_data

# this function performs the NMF ensemble method given the combined S matrix, the number of objects n, and desired
# number of clusters k
my_NMF = function(k,n,S){
 
  random_clustering = kmeans(S, centers=k)$cluster
  H = matrix(0, nrow=n, ncol=k)
  for(i in 1:n) H[i,random_clustering[i]]  = 1
  
  h_t = t(H) %*% H
  D = h_t
  
  for(j in 1:k)  h_t[j,j] = h_t[j,j]^(-1/2) 
  Q = H %*% h_t
  
  for(i in 1:50){
    Q_new = Q * sqrt(  (S %*% Q %*% D)  / (Q %*% t(Q) %*% S %*% Q %*% D) )
    
    D_new = D * sqrt( (t(Q) %*% S %*% Q) / (t(Q) %*% Q %*% D %*% t(Q) %*% Q))
    D_new = diag(diag(D_new))
    
    Q = Q_new
    D = D_new
  }
  S = Q_new %*% t(D_new) %*% t(Q_new)
  return(list(Q=Q, D=D, S=S))
}

# run the NMF
res = my_NMF(k,n,S_data)
Q = res$Q

label = rep(NA,n)
for(i in 1:n){
  label[i] = which.max(Q[i,] != 0)
}
produced_clusters = label # this is the ensemble clustering



