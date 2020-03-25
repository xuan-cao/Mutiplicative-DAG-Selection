source("multi_graph_laplace.R")
library(MASS)

GetDiagSub <- function(ind,S1){
  return(S1[ind,ind])
}
GetVecSub <- function(ind,vec){
  return(vec[ind])
}
  SpecialDet <- function(x){
    if(length(x)==0) return(1)
    if(length(x)==1) return(as.numeric(x))
    return(det(x))
  }

logposterior <- function(Adj, tau = sqrt(n), alpha1 = 0.05, alpha2 = 0.05,lambda1 = 0.05, lambda2 = p^2.001, c=2.001) {
  graph_term <- int_multi(Adj,lambda1, lambda2,c)
  pa.set <- apply(Adj,2,function(x){which(x!=0)})
  pa.count <- unlist(lapply(pa.set, length))
  pa.count <- pa.count[-p]
  pa.mat.det <- rep(0,p-1)
  pa.con <- rep(0,p-1)
  logz <- rep(0,p-1)
  for(i in 1:(p-1)){
    if(length(pa.set[[i]]) == 0) logz[i] = 0
    else{
      pa.mat.det[i] = SpecialDet(GetDiagSub(pa.set[[i]], S1))
      pa.con[i] = S1[i,i] - as.numeric(1*t(GetVecSub(as.vector(pa.set[[i]]), S1[,i]))%*%solve(GetDiagSub(as.vector(pa.set[[i]]), S1))%*%GetVecSub(pa.set[[i]], S1[,i]))
      logz[i] = (-n/2-alpha1)*log(pa.con[i]*n/2-1/(2*tau^2) + alpha2) - 1/2*log(pa.mat.det[i]) - pa.count[i]*log(sqrt(n*tau^2))
    }
  }
  return(sum(logz) + graph_term)

}

evaluation.dag <- function(Adj1, Adj2){
  true.index <- which(Adj1==1)
  false.index <- which(Adj1==0)
  positive.index <- which(Adj2==1)
  negative.index <- which(Adj2==0)
  
  TP <- length(intersect(true.index,positive.index))
  FP <- length(intersect(false.index,positive.index))
  FN <- length(intersect(true.index,negative.index))
  TN <- length(intersect(false.index,negative.index))
  
  MCC.denom <- sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  if(MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if((TN+FP)==0) MCC <- 1
  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1
  
  return(list(Precision=Precision,Recall=Recall,Sensitivity=Sensitivity,Specific=Specific,MCC=MCC,TP=TP,FP=FP,TN=TN,FN=FN))
}


DAGLasso <- function(Y, lambda, maxitr=100, tol=1e-4){
  require(lassoshooting)  
  p = ncol(Y)
  n = nrow(Y)
  S = (t(Y) %*% Y)/n
  T = diag(p)
  D = rep(1, p)
  itr_log = eps_log = NULL
  for (k in 2:p){
    nuk_old = nuk_new = c(rep(0, k-1), 1)
    r = 0
    km1_ = 1:(k-1)
    repeat {
      r = r + 1
      nuk_new[k] = D[k]
      output = lassoshooting(XtX= S[km1_,km1_,drop=FALSE], 
                             Xty=-S[km1_,k], maxit = 100, lambda = 0.5*nuk_new[k]*lambda)
      nuk_new[km1_] = output$coefficients
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
        T[k, km1_] = nuk_new[km1_]
        eps_log = c(eps_log, maxdiff)
        itr_log = c(itr_log, r)
        break
      } else {
        nuk_old = nuk_new
      }
    }
  }
  return(T)
}




DAGLassoseq <- function(Y, lambda.seq, maxitr=100, tol=1e-4){
  require(lassoshooting)  
  p = ncol(Y)
  n = nrow(Y)
  S = (t(Y) %*% Y)/n
  T = diag(p)
  D = rep(1, p)
  itr_log = eps_log = NULL
  for (k in 2:p){
    nuk_old = nuk_new = c(rep(0, k-1), 1)
    r = 0
    km1_ = 1:(k-1)
    repeat {
      r = r + 1
      nuk_new[k] = D[k]
      output = lassoshooting(XtX= S[km1_,km1_,drop=FALSE], 
                             Xty=-S[km1_,k], maxit = 100, lambda=0.5*nuk_new[k]*lambda.seq[k])
      nuk_new[km1_] = output$coefficients
      
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
        T[k, km1_] = nuk_new[km1_]
        eps_log = c(eps_log, maxdiff)
        itr_log = c(itr_log, r)
        break
      } else {
        nuk_old = nuk_new
      }
    }
  }
  Adj <- matrix(0,p,p)
  Adj[which(T!=0)] <- 1
  for(i in 1:p) Adj[i,i] = 0
  return(Adj)
}






#### SimpleGraphMove: the function to generate new graphs by adding/deleting on edge of the given graph
#### Adj: the adjacency matrix defined as before.
#### n: the number of new graphs to generate
#### checklist: a list of adjacency matrices that one wants to avoid generating. By default NULL

SimpleGraphMove <- function(Adj,n=1,checklist=NULL){
  if(!is.null(checklist)){
    result.list <- list()
    total.edge.index <- which(lower.tri(Adj))
    total.edge.number <- length(total.edge.index)
    if(total.edge.number<(n+length(checklist))) stop(paste("Cannot generate",n,"new graphs!"))
    perm.edge.index <- sample(total.edge.index,size=total.edge.number)
    checklist.edge.list <- lapply(checklist,function(x)return(which(x==1)))
    iter <- 1
    for(k in 1:total.edge.number){
      tmp.Adj <- Adj
      tmp.Adj[perm.edge.index[k]] <- 1-tmp.Adj[perm.edge.index[k]]
      tmp.edge.index <- which(tmp.Adj==1)
      if(sum(unlist(lapply(checklist.edge.list,function(x)return(setequal(x,tmp.edge.index)))))==0)
      {
        result.list[[iter]] <- tmp.Adj
        iter <- iter+1
      }
      if(iter>n) break
    }
    return(result.list)
  }
  result.list <- list()
  total.edge.index <- which(lower.tri(Adj))
  total.edge.number <- length(total.edge.index)
  if(total.edge.number<n) stop(paste("Cannot generate",n,"new graphs!"))
  perm.edge.index <- sample(total.edge.index,size=total.edge.number)
  iter <- 1
  for(k in 1:total.edge.number){
    tmp.Adj <- Adj
    tmp.Adj[perm.edge.index[k]] <- 1-tmp.Adj[perm.edge.index[k]]
    result.list[[iter]] <- tmp.Adj
    iter <- iter+1
    if(iter>n) break
  }
  return(result.list)
}



#### SSS : the stochastic shotegun search with DAG Wishart prior. The logic is shown in Jones 2005.
### x, U: defined as before
### Adj0, the starting graph adjacency matrix
### gamma: the annealing parameter in SSS of Jones 2005. It is set to be 1 by default
### X1, X3 are parameters defined the same as in Jones 2005. X1 is the graphs to explore in each step, X2 is the number of graphs to choose from for new steps (we set X2=X1 as suggested in the paper) and X3 is the number of final retained list.
### stop.num: the stopping rule. If the top X3 list has not been changed in the recent stop.num of iteration, then the algorithm stops.
### iter.max: the maximum number of iterations, if stop.num is not achieved.

greedy <- function(Adj0,iter.max=5){
  Adj.ini <- Adj0
  for(i in 1:iter.max){
    logscore.ini <- logposterior(Adj.ini, tau, alpha1, alpha2,lambda1, lambda2, c)
    Adj.remove <- Adj.ini
    Adj.remove.vector <- as.vector(Adj.remove)
    Adj.remove.vector[sample(which(Adj.remove.vector!=0),1)] <- 0
    Adj.remove <- matrix(Adj.remove.vector,p,p)
    if(logposterior(Adj.remove, tau, alpha1, alpha2,lambda1, lambda2, c)> logscore.ini) 
      Adj.ini <- Adj.remove
    Adj.add <- Adj.ini
    Adj.add[upper.tri(Adj.add)] <- 1
    diag(Adj.add) <- 1
    Adj.add.vector <- as.vector(Adj.add)
    Adj.add.vector[sample(which(Adj.add.vector==0),1)] <- 1
    Adj.add <- matrix(Adj.add.vector,p,p)
    Adj.add[upper.tri(Adj.add)] <- 0
    diag(Adj.add) <- 0
    #print("hh")
    #print(Adj.add)
    if(logposterior(Adj.add, tau, alpha1, alpha2,c) > logscore.ini) Adj.ini <- Adj.add
  }
  return(Adj.ini)
}




SSS <- function(Adj0,adaptive=TRUE,gamma=1,X1,X3,stop.num,iter.max){
  ## basic checking
  #if(X2>X1) stop("X2 cannot be larger than X1!")
  
  ## object setup
  final.list <- list(Graph=list(),logScore=rep(0,X3),MinIndex=0,Length=0)
  iter.num <- 0
  nochange.num <- 0
  pre.MinScore <- NULL
  
  
  
  ini.score <- logposterior(Adj0, tau, alpha1, alpha2,lambda1, lambda2, c)
  # ini.score <- logDAGMargScore(x,Adj0,U=diag(p),c=10,d=d)
  final.list$Graph[[1]] <- Adj0
  final.list$logScore[1] <- ini.score
  final.list$MinIndex <- 1
  final.list$Length <- 1
  current.Adj <- Adj0
  pre.MinScore <- ini.score
  DAG.checklist <- final.list$Graph
  
  while((nochange.num < stop.num) && (iter.num < iter.max)){
    ## Generate new graphs
    if(X3==0) DAG.checklist <- NULL
    tmp.graphs <- SimpleGraphMove(current.Adj,n=X1,checklist = DAG.checklist)
    tmp.score <- unlist(lapply(tmp.graphs,function(Adj){
      score <- logposterior(Adj, tau, alpha1, alpha2,lambda1, lambda2, c)
      return(score)
    }))
    ## update candidate list
    for(k in 1:X1){
      final.list <- UpdateGraphList(final.list,tmp.graphs[[k]],tmp.score[k],max.length=X3)
    }
    DAG.checklist <- final.list$Graph
    ## select the graph to start next search
    new.index <- ScoreSampling(tmp.score,gamma)
    current.Adj <- tmp.graphs[[new.index]]
    iter.num <- iter.num + 1
    current.MinScore <- final.list$logScore[final.list$MinIndex]
    nochange.num <- nochange.num + 1
    if(pre.MinScore != current.MinScore) nochange.num <- 0
    pre.MinScore <- current.MinScore
  }
  final.list$iter <- iter.num
  final.list$nochange <- nochange.num
  return(final.list)
}


########################
#### UpdateGraphList: function to update candidate list according to current evaluation. For internal use only.
### glist: the candidate list to update
### Adj: adjacency matrix to consider
### score: log marginal score of the graph
### max.length: the limit for the list
UpdateGraphList <- function(glist,Adj,score,max.length){
  if(glist$Length <  max.length){
    glist$Graph[[glist$Length+1]] <- Adj
    glist$logScore[glist$Length+1] <- score
    glist$Length <- glist$Length+1
    glist$MinIndex <- which.min(glist$logScore)
    #print("reach")
    return(glist)
  }
  if(glist$logScore[glist$MinIndex]<score){
    glist$Graph[[glist$MinIndex]] <- Adj
    glist$logScore[[glist$MinIndex]] <- score
    glist$MinIndex <- which.min(glist$logScore)
    return(glist)
  }
  return(glist)
}


#### ScoreSampling : internal use - function to sample the index to next graph according to the rule in Jones 2005. Special treatment for numerical stability is used.
### logScore: the sequence of log scores
### beta: annealing parameter

log.sum.exp <- function(x){
  maxx <- max(x)
  return(log(sum(exp(x-maxx)))+maxx)
}

ScoreSampling <- function(logScore,gamma=1){
  lscore.vec <- gamma*logScore
  lscore.vec <- lscore.vec - (min(lscore.vec))
  log.score.sum <- log.sum.exp(lscore.vec)
  p.weight <- exp(lscore.vec-log.score.sum)
  p.weight <- p.weight/sum(p.weight)
  select.index <- sample(1:length(logScore),size=1,prob=p.weight)
  return(select.index)
}



