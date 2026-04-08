#' Title
#'
#' @param Y ...
#' @param L ...
#' @param K ...
#' @param seed ...
#' @param hyperpar ...
#' @param warmstart ...
#' @param BY ...
#'
#' @return ...
#' @export
#'
#' @examples
#'  1+1
random_init_poseidon = function(Y, # must be a list
                            L,
                            K,
                            seed = NULL,
                            hyperpar,
                            warmstart,
                            BY = .1){
  if(!is.null(seed)) set.seed(seed)
  Tr = length(Y)
  J = ncol(Y[[1]])
  MKCD0 = MKCD_init = array(NA,c(L,4,Tr))
  a  = hyperpar$a
  b  = hyperpar$b
  m0 = hyperpar$m0
  k0 = hyperpar$k0
  c0 = hyperpar$c0
  d0 = hyperpar$d0
  for(t in 1:Tr){
  if(warmstart){
    ml  <- stats::kmeans(c(Y[[t]]),centers = L)$centers
  }else{
    ml  <- stats::runif(L, min(c(Y[[t]])),max(c(Y[[t]])))
  }
  kl  <- stats::rgamma(L,hyperpar$k0,1)
  al  <- stats::rgamma(L,hyperpar$c0,1)
  bl  <- stats::rgamma(L,hyperpar$d0,1)

  MKCD_init[,,t] = cbind(ml,kl,al,bl)
  MKCD0temp = cbind(c(rep(m0,L)),rep(k0,L),rep(c0,L),rep(d0,L))
  # MKCD0temp[1,] = c(0,100,100,1)
  MKCD0[,,t] = MKCD0temp
  }

  if(!is.null(seed)) set.seed(seed)
  a = rep(a,K)
  b = array(b,c(L,K,Tr))

  if(K<J){
    log.RHO_jk <- log(matrix(stats::rbeta(J*K,1,1),J,K))
    Z2         <- apply(log.RHO_jk, c(1), function(x) matrixStats::logSumExp(x))
    R     <- exp(sapply(1:K, function(qq) log.RHO_jk[,qq]-Z2,simplify = "matrix"))
  }else if(K==J){
    R     <- diag(J)
  }else{
    R     <- cbind(diag(J), matrix(0,J,K-J))
  }
  XXX = list()
  for(t in 1:Tr){
    N = nrow(Y[[t]])
    XI_ikl = array(pnorm(rnorm(N*K*L,0,10),0,10),c(N,K,L))
    XI_ikl = apply(XI_ikl,c(1,2),function(x)x/sum(x))
    XXX[[t]] = aperm(XI_ikl, c(2,3,1))
  }
  return(list(a=a,b=b,
              X=XXX,R=R,
              MKCD0 = MKCD0, MKCD_init = MKCD_init))
}
















#' Title
#'
#' @param Y ...
#' @param L ...
#' @param K ...
#' @param seed ...
#' @param hyperpar ...
#' @param warmstart ...
#' @param BY ...
#'
#' @return ...
#' @export
#'
#' @examples
#'  1+1
random_init_poseidon_v2 = function(Y, # must be a list
                                L,
                                K,
                                seed = NULL,
                                hyperpar,
                                warmstart,
                                BY = .1){
  if(!is.null(seed)) set.seed(seed)
  Tr = length(Y)
  J = ncol(Y[[1]])
  MKCD0 = MKCD_init = array(NA,c(L,4,Tr))
  a  = hyperpar$a
  b  = hyperpar$b
  m0 = hyperpar$m0
  k0 = hyperpar$k0
  c0 = hyperpar$c0
  d0 = hyperpar$d0
  for(t in 1:Tr){
    if(warmstart){
      ml  <- stats::kmeans(c(Y[[t]]),centers = L)$centers
    }else{
      ml  <- stats::runif(L, min(c(Y[[t]])),max(c(Y[[t]])))
    }
    kl  <- stats::rgamma(L,hyperpar$k0,1)
    al  <- stats::rgamma(L,hyperpar$c0,1)
    bl  <- stats::rgamma(L,hyperpar$d0,1)

    MKCD_init[,,t] = cbind(ml,kl,al,bl)
    MKCD0temp = cbind(c(rep(m0,L)),rep(k0,L),rep(c0,L),rep(d0,L))
    # MKCD0temp[1,] = c(0,100,100,1)
    MKCD0[,,t] = MKCD0temp
  }

  if(!is.null(seed)) set.seed(seed)
  a = rep(a,K)
  b = array(b,c(L,K,Tr))

  if(K<J){
    log.RHO_jk <- log(matrix(stats::rbeta(J*K,1,1),J,K))
    Z2         <- apply(log.RHO_jk, c(1), function(x) matrixStats::logSumExp(x))
    R     <- exp(sapply(1:K, function(qq) log.RHO_jk[,qq]-Z2,simplify = "matrix"))
  }else if(K==J){
    R     <- diag(J)
  }else{
    R     <- cbind(diag(J), matrix(0,J,K-J))
  }

  col_cl <- apply(R,1,which.max )

  XXX = list()

  for(t in 1:Tr){
    N = nrow(Y[[t]])
    XI_ikl <- array(NA,c(N,K,L))

    for(k in 1:K){
      inds <- which(col_cl==k)
      if(length(inds) == 0 ){
        XI_ikl[,k,] = matrix(1/L,N,L)
      }else{
        cl = e1071::cmeans(as.matrix(Y[[t]][,inds],nrow=N),L)
        XI_ikl[,k,] = cl$membership
      }
    }
    XXX[[t]] = XI_ikl
  }
  return(list(a=a,b=b,
              X=XXX,R=R,
              MKCD0 = MKCD0, MKCD_init = MKCD_init))
}

