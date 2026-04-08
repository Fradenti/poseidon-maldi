#' Title
#'
#' @param Y
#' @param L
#' @param K
#' @param seed
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
random_init_pose_RHO = function(Y,
                            L,
                            K,
                            seed = NULL,
                            hyperpar,
                            BY = .1,
                            p = .2,
                            perturb = .5){


  if(!is.null(seed)) set.seed(seed)

  J = ncol(Y)
  N = nrow(Y)

  # ----------------------------------------------------------------------------
  a = rep(hyperpar$a,K)
  b = matrix(hyperpar$b,L,K)
  # ----------------------------------------------------------------------------

  if(is.null(warmMAT)){
    # No warmstart, uniform allocation -----------------------------------------
    if(K<J){
      log.RHO_jk <- log(matrix(rbeta(J*K,1,1),J,K))
      Z2         <- apply(log.RHO_jk, c(1), function(x) matrixStats::logSumExp(x))
      R     <- exp(sapply(1:K, function(qq) log.RHO_jk[,qq]-Z2,simplify = "matrix"))
    }else if(K==J){
      R     <- diag(J)
    }else{
      R     <- cbind(diag(J), matrix(0,J,K-J))
    }
    XI_ikl = array(pnorm(rnorm(N*K*L,0,10),0,10),c(N,K,L))
    XI_ikl = apply(XI_ikl,c(1,2),function(x)x/sum(x))
    X = aperm(XI_ikl, c(2,3,1))
  }else{
    # Warmstart, perturb the warmstart allocations -----------------------------
    R = apply(warmMAT$warmR + runif(K*ncol(Y),-perturb,perturb),
              2,
              function(x) abs(x)/sum(abs(x)))

    for(k in 1:K){
      QQ = warmMAT$warmX[,k,] +runif(L*nrow(Y),-perturb,perturb)
      warmMAT$warmX[,k,] = apply(QQ,1,function(x) abs(x)/sum(abs(x)))
    }
    X = warmMAT$warmX
  }

  return(list(a=a,b=b,X=X,R=R))
}



#' Title
#'
#' @param Y
#' @param L
#' @param K
#' @param seed
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
random_init_pose_v2 = function(Y,
                            L,
                            K,
                            seed = NULL,
                            warmMAT,
                            hyperpar,
                            BY = .1,
                            p = .2,
                            perturb = .5){
  if(!is.null(seed)) set.seed(seed)

  J = ncol(Y)
  N = nrow(Y)

  # ----------------------------------------------------------------------------
  a = rep(hyperpar$a,K)
  b = matrix(hyperpar$b,L,K)
  # ----------------------------------------------------------------------------

  if(is.null(warmMAT)){
    # No warmstart, uniform allocation -----------------------------------------
  #  if(K<J){
      log.RHO_jk <- log(matrix(rbeta(J*K,1,1),J,K))
      Z2         <- apply(log.RHO_jk, c(1), function(x) matrixStats::logSumExp(x))
      R     <- exp(sapply(1:K, function(qq) log.RHO_jk[,qq]-Z2,simplify = "matrix"))
    # }else if(K==J){
    #   R     <- diag(J)
    # }else{
    #   R     <- cbind(diag(J), matrix(0,J,K-J))
    # }
    col_cl <- apply(R,1,which.max )
    XI_ikl <- array(NA,c(N,K,L))

    for(k in 1:K){
      inds <- which(col_cl==k)
      if(length(inds) == 0 ){
        XI_ikl[,k,] = matrix(1/L,N,L)
      }else{
      cl = e1071::cmeans(as.matrix(Y[,inds],nrow=N),L)
      XI_ikl[,k,] = cl$membership
      }
    }

    X = XI_ikl
  }else{
    # Warmstart, perturb the warmstart allocations -----------------------------
    R = apply(warmMAT$warmR + runif(K*ncol(Y),-perturb,perturb),
              2,
              function(x) abs(x)/sum(abs(x)))

    for(k in 1:K){
      QQ = warmMAT$warmX[,k,] +runif(L*nrow(Y),-perturb,perturb)
      warmMAT$warmX[,k,] = apply(QQ,1,function(x) abs(x)/sum(abs(x)))
    }
    X = warmMAT$warmX
  }
  # Prior parameters on atoms --------------------------------------------------
  MKCD0 = cbind(c(rep(hyperpar$m0,L)),
                rep(hyperpar$k0,L),
                rep(hyperpar$c0,L),
                rep(hyperpar$d0,L))

  return(list(a=a,b=b,X=X,R=R, MKCD0 = MKCD0))
}



#' Title
#'
#' @param Y
#' @param L
#' @param K
#'
#' @return
#' @export
#'
#' @examples
warm_clust_pose <- function(Y,L,K){

  seqq = seq(min(Y)-1,max(Y)+1,length.out = 100)
  Q <- parallel::mclapply(1:ncol(Y), function(y)
    ecdf(Y[,y])(seqq),mc.cores = 5)
  M <-  do.call(rbind,(Q))

  Dend = hclust(as.dist(Rfast::Dist(M)))
  # distance matrix between ecdfs
  QQ = model.matrix(~-1+as.factor(cutree(Dend,K)))
  R = t(apply(QQ,1,function(x) abs(x)/sum(abs(x))))
  colnames(R) = NULL
  col_cl_init = apply(R,1,which.max)
  # within the distributional clusters, assign row clusters
  XI_ikl = array(NA,c(nrow(Y),K,L))
  for(kk in 1:K){
    s = kmeans(Y[,col_cl_init==kk],centers = L)$cl
    QQ = model.matrix(~-1+as.factor(s))
    X = t(apply(QQ,1,function(x) abs(x)/sum(abs(x))))
    XI_ikl[,kk,] = X
  }
  X=XI_ikl
  return(list(warmR=R,warmX=X))
}


#' @export
warm_clust_pose_double_kmeans <- function(Y,L,K){

  kcol = kmeans(t(Y),K)

  # distance matrix between ecdfs
  QQ = model.matrix(~-1+as.factor(kcol$cluster))
  R = t(apply(QQ,1,function(x) abs(x)/sum(abs(x))))
  colnames(R) = NULL
  col_cl_init = apply(R,1,which.max)
  # within the distributional clusters, assign row clusters
  XI_ikl = array(NA,c(nrow(Y),K,L))
  for(kk in 1:K){
    s = kmeans(Y[,col_cl_init==kk],centers = L)$cl
    QQ = model.matrix(~-1+as.factor(s))
    X = t(apply(QQ,1,function(x) abs(x)/sum(abs(x))))
    XI_ikl[,kk,] = X
  }
  X=XI_ikl
  return(list(warmR=R,warmX=X))
}




#' @export
#'
#' @examples
warm_clust_pose2 <- function(Y,L0,L,K){

  ce <- cl <- Y

  for(j in 1:ncol(Y)){
    KM <- kmeans(Y[,j], L0)
    cl[,j] <- KM$clu
    ce[,j] <- KM$centers[KM$clu]
  }

  DD = Rfast::Dist(t(ce))
  Dend = hclust(as.dist(DD))
  # distance matrix between ecdfs
  QQ = model.matrix(~-1+as.factor(cutree(Dend,K)))
  R = t(apply(QQ,1,function(x) abs(x)/sum(abs(x))))
  colnames(R) = NULL
  col_cl_init = apply(R,1,which.max)

  # within the distributional clusters, assign row clusters
  XI_ikl = array(NA,c(nrow(Y),K,L))
  for(kk in 1:K){
    s = kmeans(Y[,col_cl_init==kk],centers = L)$cl
    QQ = model.matrix(~-1+as.factor(s))
    X = t(apply(QQ,1,function(x) abs(x)/sum(abs(x))))
    XI_ikl[,kk,] = X
  }
  X=XI_ikl
  return(list(warmR=R,warmX=X))
}



#' Title
#'
#' @param partition
#' @param row_part
#' @param K
#' @param L
#' @param perturb
#'
#' @return
#' @export
#'
#' @examples
from_partition_to_probs <- function(partition,
                                    K, perturb=0){

  cp  <- factor(partition,levels = 1:K)
  RHO <- model.matrix(~-1+cp)
  if(perturb>0){
    RHOu <- RHO + runif(prod(dim(RHO)),0,perturb)
    RHO  <- t(apply(RHOu,1,function(x) x/sum(x)))
    }
  return(RHO = RHO)
}
