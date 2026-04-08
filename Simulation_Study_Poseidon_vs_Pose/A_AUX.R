library(tidyverse)
library(Poser)

# Functions ---------------------------------------------------------------

make_mix <- function(ns,mus,sigs){
  l <- as.list(c(sapply(1:length(ns), function(y) rnorm(ns[y],mus[y],sigs[y]))))
  do.call(c,l)
}

plot_elbo_poseidon <- function(XXX,RUNS=250){
  par(mfrow=c(1,2))
  plot(XXX[[1]]$Elbo,type="l")
  for(i in 1:RUNS){
    lines(XXX[[i]]$Elbo,col=i)
  }
  plot(diff(XXX[[1]]$Elbo),type="l",ylim=c(-5,4))
  for(i in 1:RUNS){
    lines(diff(XXX[[i]]$Elbo),col=i)
    if(!all(diff(XXX[[i]]$Elbo)>0)) cat(paste(i,"\n"))
  }
  abline(h=0)
  par(mfrow=c(1,1))
}


f <- function(i,Y){
  set.seed(14332*i)

  par <- maotai::kmeanspp(t(Y),k=K)

  RHOi <- from_partition_to_probs(partition = par,
                                  K = K,
                                  perturb = .5)
  N <- nrow(Y)
  XI_ikl <- array(NA,c(N,K,L))

  for(k in 1:K){
    par_row <- maotai::kmeanspp(as.matrix((Y[,par==k])),k=5)
    XI_ikl[,k,] = from_partition_to_probs(par_row, K=L, perturb = .2)
  }

  VV =   Poser::SEr_fiSAN(
    Y = Y,
    L = L,
    K = K,
    hyperpar = hyperpar,
    epsilon = 1e-8,
    RHO_init = RHOi,
    XI_init = XI_ikl,
    nsim = 1500,verbose = F)
  VV
}


g <- function(i,YY){

  set.seed(3121*i)

  par <- maotai::kmeanspp(t(do.call(rbind,YY)),k=K)

  RHOi <- Poser::from_partition_to_probs(partition = par,
                                  K = K,
                                  perturb = .5)
  XIi = list()

  col_cl <- par

  XXX <- list()

  for(tt in 1:length(YY)){
    N <- nrow(YY[[tt]])
    XI_ikl <- array(NA,c(N,K,L))


    for(k in 1:K){
      par_row <- maotai::kmeanspp(as.matrix(YY[[tt]][,par==k]),k=5)
      XI_ikl[,k,] = Poser::from_partition_to_probs(par_row, K=L, perturb = .2)
    }
    XXX[[tt]] = XI_ikl
  }

  VV =  POSEIDON::POSEIDON_onlyfiSAN(Y_list = YY,
                                 RHOi = RHOi,
                                 XIi = XXX,
                                 L = L,
                                 K = K,verbose = F,
                                 epsilon = 1e-8,
                                 nsim = 1500,
                                 hyperpar = hyperpar)
  VV
}
