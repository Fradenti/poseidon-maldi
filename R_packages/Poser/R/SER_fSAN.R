#' SER: only finite mixtures
#'
#' @param Y the dataset to bicluster, N x J
#' @param L the observational clusters (rows) upper bound
#' @param K the distributional clusters (cols) upper bound
#' @param nsim max number of iterations that is considered
#' @param epsilon threshold for stoping rule
#' @param compute_elbo_every how many iterations between elbo evaluations?
#' @param verbose should I print the evolution?
#' @param hyperpar a list with hyperparamaters for the model (see examples)
#'
#' @return
#' @export
#'
#' @examples
SEr_fSAN <- function(Y,
                        L,
                        K,
                        nsim = 50,
                        epsilon = 1e-4,
                        hyperpar,
                        compute_elbo_every = 1,
                        RHO_init,
                        XI_init,
                        verbose = TRUE){
  # ----------------------------------------------------------------------------
  MKCD0 = as.matrix(cbind(c(rep(hyperpar$m0,L)),
                rep(hyperpar$k0,L),
                rep(hyperpar$c0,L),
                rep(hyperpar$d0,L)))

  t1 <- Sys.time()
  res <- SE_fSAN(Y = Y,
                 L = L, K = K,
                 A0 = rep(hyperpar$a,K),
                 B0 = matrix(hyperpar$b,L,K),
                 mkcd_0 = MKCD0,
                 RHO_init = RHO_init,
                 XI_init = XI_init,
                 nsim = nsim,
                 compute_elbo_every = compute_elbo_every,
                 epsilon = epsilon,
                 verbose = verbose)
  t2 <- Sys.time()
  res[["model"]] <- "fiSAN"
  res[["CEE"]] <- compute_elbo_every
  res[["time"]] <- t2-t1
  structure(res,class = c("poser", class(res)))

}
