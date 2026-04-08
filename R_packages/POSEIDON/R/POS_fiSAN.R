#' Title
#'
#' @param L ...
#' @param K ...
#' @param NN_inds ...
#' @param nsim ...
#' @param replicates_potts ...
#' @param epsilon ...
#' @param seed ...
#' @param hyperpar ...
#' @param warmstart ...
#' @param Y_list ...
#' @param itemp ...
#' @param itemp_grid ...
#' @param compute_elbo_every ...
#' @param IN ...
#' @param verbose ...
#' @param estimate_itemp ...
#'
#' @return ...
#' @export
#'
#' @examples
#' 1+1
POSEIDON_fiSAN <- function(Y_list,
                     L,
                     K,
                     NN_inds,
                     XIi, RHOi,
                     itemp = 0,
                     itemp_grid = seq(0,2,by=.1),
                     nsim = 50,
                     replicates_potts = 5,
                     epsilon = 1e-4,
                     hyperpar,
                     seed = NULL,
                     verbose = TRUE,
                     estimate_itemp = FALSE,
                     compute_elbo_every = 1
){

  cols <- unlist(lapply(Y_list, function(x) ncol(x)))
  if(!all(cols/min(cols)==1)){
    stop("We need same number of columns across datasets")
  }

      t1 <- Sys.time()
      if( is.null(hyperpar$s1) | is.null(hyperpar$s1)){
        stop("You need to specify either s1 or s2 in the hyperpar object")
      }

      Tr <- length(Y_list)
      b  <- hyperpar$b
      m0 <- hyperpar$m0
      k0 <- hyperpar$k0
      c0 <- hyperpar$c0
      d0 <- hyperpar$d0
      a  <- rep(hyperpar$a,K)
      b  <- array(hyperpar$b,c(L,K,Tr))

      MKCD0 = array(NA,c(L,4,Tr))
      for(t in 1:Tr){
        MKCD0[,,t] <- cbind(c(rep(m0,L)),rep(k0,L),rep(c0,L),rep(d0,L))
      }

      res <- POSEIDON_cpp_fiSAN(Y_list = Y_list,
                               NN_inds = NN_inds,
                               itemp = itemp,
                               itemp_grid = itemp_grid,
                               L = L,
                               K = K,
                               D = length(Y_list),
                               B0 = b,
                               mkcd_0 = MKCD0,
                               RHO_init = RHOi,
                               XI_init = XIi,
                               nsim = nsim,
                               conc_hyper = c(hyperpar$s1,hyperpar$s2),
                               compute_elbo_every = compute_elbo_every,
                               epsilon = epsilon,
                               replicates_potts = replicates_potts,
                               verbose = verbose,
                               estimate_itemp = estimate_itemp)
      t2 <- Sys.time()
      res[["model"]] <- "fiSAN"
      res[["CEE"]] <- compute_elbo_every
      res[["time"]] <- t2-t1
      structure(res,class = c("POSEIDON", class(res)))


}
