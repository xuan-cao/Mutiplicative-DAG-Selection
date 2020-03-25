rm(list=ls())
source("function_multi_new.R")
library(reshape2)
library(ggplot2)
# generate true data
n <- 150
p <- 3*n
c = 2.001
tau = sqrt(n)
alpha1 = 0.05
alpha2 = 0.05
spars <- 3
lambda1 <- 0.05
lambda2 <- p^c



true.L <- matrix(1,p,p)
true.L[upper.tri(true.L)] <- 0
for(t in 1:(p-1))
{
  true.L[(t+1):p,t] <- true.L[(t+1):p,t]*rbinom(p-t,1,min(spars/(p-t),1))
}
nonzero <- as.vector(1*colSums(true.L))
for(w in 1:(p-1)) {
  if(nonzero[w] == 1) true.L[p,w] = 1
}
true.Adj <- true.L
diag(true.L) <- 1
diag(true.Adj) <- 0

true.Sigma <- solve(true.L%*%t(true.L))



X <- mvrnorm(n,rep(0, p),Sigma=true.Sigma)
S <- 1/n*t(X)%*%X
S1 <- S + 1/(n*tau^2)*diag(p)


epsilon <- seq(0.2,0.2,by=0.01) # can tune the epsilon value
epsilonmax <- list(maxAdj = list(), maxlog = list())
time = Sys.time()
for (u in 1:length(epsilon)) {
  print(epsilon[u])
  ch <- chol(solve(S+unlist(epsilon[u])*diag(p))) 
  dd <- diag(ch) 
  L.epsilon <- t(ch/dd)
  #################################################################
  thres <- seq(0.1,0.5,by=0.1)
  nonzero.adj <- rep(0,p)
  logscore <- rep(0, length(thres))
  for(k in 1:length(thres)){
    cat('thres:', thres[k],'\n')
    L.epsilon1 <- L.epsilon
    for(i in 2:p){
      for(j in 1:(i-1)){
        if(abs(L.epsilon[i,j]) < thres[k]) L.epsilon1[i,j] <- 0}}
    nonzero.adj <- as.vector(1*colSums(L.epsilon1!=0))
    for(m in 1:(p-1)) {
      if(nonzero.adj[m] == 1) L.epsilon1[p,m] = 1
    }
    Adj.epsilon <- matrix(0,p,p)
    Adj.epsilon[which(L.epsilon1 !=0)] <- 1
    diag(Adj.epsilon) <- 0
    logscore[k] <- logposterior(Adj=Adj.epsilon,tau,alpha1, alpha2,lambda1, lambda2,c)
    print(logscore[k])
    
  }
  
  thres_opt <- thres[which.max(logscore)]


}


for(i in 2:p){
  for(j in 1:(i-1)){
    if(abs(L.epsilon[i,j]) < thres_opt) L.epsilon[i,j] <- 0}}
diag(L.epsilon) <- 0
L.epsilon <- 1*(L.epsilon != 0)


####Evaluate the model performance
print(evaluation.dag(true.Adj,L.epsilon))





#######run the SSS or greedy search (optional)

#greedDAG <- greedy(Adj0 = L.epsilon, iter.max = 5)
#logposterior(Adj=L.epsilon,tau,alpha1, alpha2,lambda1, lambda2,c)
#StepList <- SSS(Adj0 = L.epsilon,adaptive=TRUE,gamma=1,X1=1,X3=2,stop.num=2,iter.max=2)
#print(evaluation.dag(true.Adj,StepList$Graph[[which.max(StepList$logScore)]]))
