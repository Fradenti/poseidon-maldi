# using kmeans(nstart = 20) instead of kmeanspp
f <- function(i,Y,K,L,hyperpar, eps1= .05, eps2= .2){
  set.seed(7325*i)

  # km_out <- maotai::kmeanspp(t(Y),k=K)
  par <- kmeans(t(Y),centers=K, nstart = 20,iter.max = 5000,
                algorithm = "Hartigan-Wong")$cluster
  RHOi <- Poser::from_partition_to_probs(par,
                                         K = K,
                                         perturb = eps1)
  p <- nrow(Y)
  XI_ikl <- array(NA,c(p,K,L))
  for(k in 1:K){
    par_row <- kmeans(Y[,par==k, drop=FALSE],centers=L, nstart = 25)$cluster
    XI_ikl[,k,] = Poser::from_partition_to_probs(par_row, K=L, perturb = eps2)
  }

  VV = POSEr_fiSAN(
    Y = Y,
    L = L,
    K = K,
    NN_inds = nnList_cpp,
    itemp = 1,
    estimate_itemp = FALSE,
    # itemp_range = seq(0,2,by=.01),
    epsilon = 1e-8,
    hyperpar = hyperpar,
    RHO_init = RHOi,
    XI_init = XI_ikl,
    nsim = 2500,verbose = F,
    replicates_potts = 5)
  VV
}

# specify starting_point
g <- function(i,Y,K,L,hyperpar, starting_point, eps1= .05, eps2= .2){
  set.seed(7325*i)

  if(!is.null(starting_point)){
    par = starting_point[[1]]
  } else {
    par <- kmeans(t(Y),centers=K, nstart = 20,
                  algorithm = "Hartigan-Wong")$cluster
  }
  K <- length(unique(par))
  RHOi <- Poser::from_partition_to_probs(par,
                                         K = K,
                                         perturb = eps1)
  p <- nrow(Y)
  XI_ikl <- array(NA,c(p,K,L))
  for(k in 1:K){
    if(!is.null(starting_point)){
      par_row <- starting_point[[2]][[k]]
    } else{
      par_row <- kmeans(Y[,par==k, drop=FALSE],centers=L, nstart = 25)$cluster
    }
    XI_ikl[,k,] = Poser::from_partition_to_probs(par_row, K=L, perturb = eps2)
  }

  VV = POSEr_fiSAN(
    Y = Y,
    L = L,
    K = K,
    NN_inds = nnList_cpp,
    itemp = 1,
    estimate_itemp = FALSE,
    # itemp_range = seq(0,2,by=.01),
    epsilon = 1e-7,
    hyperpar = hyperpar,
    RHO_init = RHOi,
    XI_init = XI_ikl,
    nsim = 1500,verbose = F,
    replicates_potts = 5)
  VV
}

mapper1 <- function(Y){
  z = Y/max(Y)
  z - mean(z)
}
mapper1_ext <- function(data,Y){
  data = data/max(Y)
  z = Y/max(Y)
  data - mean(z)
}
