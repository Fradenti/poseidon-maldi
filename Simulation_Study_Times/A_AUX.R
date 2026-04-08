library(tidyverse)
library(Poser)

# Functions ---------------------------------------------------------------

make_mix <- function(ns,mus,sigs){
  l <- as.list(c(sapply(1:length(ns), function(y) rnorm(ns[y],mus[y],sigs[y]))))
  do.call(c,l)
}

g_time_zeroepsilon_10iterations <- function(i,YY){

  set.seed(3121*i)

  par <- maotai::kmeanspp(t(do.call(rbind,YY)),k=K)

  RHOi <- Poser::from_partition_to_probs(partition = par,
                                  K = K,
                                  perturb = .5)

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

  VV =  POSEIDON::POSEIDON_onlyfiSAN(
                                 Y_list = YY,
                                 RHOi = RHOi,
                                 XIi = XXX,
                                 L = L,
                                 K = K,
                                 verbose = FALSE,
                                 epsilon = 0,
                                 nsim = 10,
                                 hyperpar = hyperpar)
  VV
}

# Single datasets -------------------------------------------------------------------------
fisan <- function(i,Y){
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
    epsilon = 0,
    RHO_init = RHOi,
    XI_init = XI_ikl,
    nsim = 10,verbose = F)
  VV
}


f_cam <- function(i,Y){
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

  VV =   Poser::SEr_CAM(
    Y = Y,
    L = L,
    K = K,
    hyperpar = hyperpar,
    epsilon = 0,
    RHO_init = RHOi,
    XI_init = XI_ikl,
    nsim = 10,
    verbose = F)
  VV
}

