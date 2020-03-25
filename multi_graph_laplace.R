int_multi <- function(Adj, lambda1, lambda2,c){
  log_l <- function(Adj, lambda1, lambda2,c, u){
    omega <- exp(u)/(1+exp(u))
    new_adj <- Matrix::forceSymmetric(Adj,uplo="L")
    pa.set <- apply(new_adj,2,function(x){which(x!=0)})
    pa.count <- unlist(lapply(pa.set, length))
    lulu <- (1 - Adj) * log(1 - omega %*% t(omega))
    lulu <- lulu[lower.tri(lulu)]
    sum(lulu)
    return(sum((lambda1 + pa.count)*u - (lambda1 + pa.count + 1)*log(1 + exp(u)) + (lambda2 - 1)*log(1-omega))) + lulu - p*log(beta(lambda1, lambda2))}
  
  hessian_matrix <- function(Adj, lambda1, lambda2,c, u){
    omega <- exp(u)/(1+exp(u))
    new_adj <- Matrix::forceSymmetric(Adj,uplo="L")
    pa.set <- apply(new_adj,2,function(x){which(x!=0)})
    pa.count <- unlist(lapply(pa.set, length))
    di <- omega / (1 + exp(u))
    didi <- di * (1 - exp(u)) / (1 + exp(u))
    lulu <- - (1 - Adj) * (di %*% t(di)) / ((1 - omega %*% t(omega))^2)
    lala <- - (lambda1 + pa.count + 1) * di - (lambda2 - 1) * (didi / (1 - omega) + (di / (1 - omega))^2)
    momo <- (1 - new_adj) * ((didi %*% t(omega)) / (1 - omega %*% t(omega)) + ((di %*% t(omega)) / (1 - omega %*% t(omega)))^2)
    diag(momo) <- 0
    diag(lulu) <- lala - rowSums(as.matrix(momo))
    return(lulu)
  }
  initial <- rep(p^(-c), p)
  wrapper <- function(u1) -log_l(Adj, lambda1, lambda2,c, u1)
  o = nlm(wrapper, initial, steptol = 1e-3, iterlim = 5)
  f_x = o$minimum 
  u0 = o$estimate
  ccc = 0
  H = sum(log(abs(eigen(hessian_matrix(Adj, lambda1, lambda2,c, u0))$values)))
  int = ccc-0.5*H-f_x
  return(int)
}