#' POSER: potts model + CAM
#'
#' @param Y the dataset to bicluster, N x J
#' @param L the observational clusters (rows) upper bound
#' @param K the distributional clusters (cols) upper bound
#' @param NN_inds a list of J elements. Each element contains a vector with C++ indexes (numbers starting from 0) indicating the neighbors
#' @param nsim max number of iterations that is considered
#' @param replicates_potts number of sweeps of responsibilities updates at each iteration of the algorithm
#' @param epsilon threshold for stoping rule
#' @param compute_elbo_every how many iterations between elbo evaluations?
#' @param itemp,estimate_itemp,itemp_range starting point (if estimate_itemp = TRUE) or fixed value (if estimate_itemp = FALSE) for the
#' inverse temperature parameter. If estimated, itemp_range, a vector containing min and max, is needed
#' @param IN a list with INitial conditions
#' @param verbose should I print the evolution?
#' @param hyperpar a list with hyperparamaters for the model (see examples)
#'
#' @return
#' @export
#'
#' @examples
POSEr_CAM <- function(Y,
                        L,
                        K,
                        NN_inds,
                        itemp = 0.25,
                        estimate_itemp = F,
                        itemp_range = c(0,5),
                        nsim = 50,
                        replicates_potts = 5,
                        epsilon = 1e-4,
                        hyperpar,
                        compute_elbo_every = 1,
                        RHO_init,
                        XI_init,
                        verbose = TRUE){
  # ----------------------------------------------------------------------------
  MKCD0 = cbind(c(rep(hyperpar$m0,L)),
                  rep(hyperpar$k0,L),
                  rep(hyperpar$c0,L),
                  rep(hyperpar$d0,L))
  t1 <- Sys.time()
  res <- POSE_CAM(Y = Y,
                    NN_inds = NN_inds,
                    itemp = itemp,
                    estimate_itemp = estimate_itemp,
                    itemp_range = itemp_range,
                    L = L, K = K,
                    conc_hyper = c(hyperpar$s1,hyperpar$s2,
                                   hyperpar$s3,hyperpar$s4),
                    mkcd_0 = MKCD0,
                    RHO_init = RHO_init,
                    XI_init  = XI_init,
                    nsim = nsim,
                    compute_elbo_every = compute_elbo_every,
                    epsilon = epsilon,
                    replicates_potts = replicates_potts,
                    verbose = verbose)
  t2 <- Sys.time()
  res[["model"]] <- "CAM"
  res[["CEE"]] <- compute_elbo_every
  res[["time"]] <- t2-t1
  structure(res,class = c("poser", class(res)))

}
