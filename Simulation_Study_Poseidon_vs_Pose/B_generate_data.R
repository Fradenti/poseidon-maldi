source("A_Simulation_1/A_AUX.R")
# Generate data -----------------------------------------------------------

D0list <- list()
D1list <- list()
D2list <- list()
D3list <- list()
D4list <- list()
CClist <- list()
RClist <- list()

for(i in 1:100){

  set.seed(1221*i)
  # reference dataset
  s_d1 <- 1
  NCOL <- 10
  NOBS <- 25

  D0 <- cbind(
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(4,0),sigs = c(s_d1,s_d1))),
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(0,0),sigs = c(s_d1,s_d1))))

  # high noise
  s_d1 <- 1
  D1 <- cbind(
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(.25,0),sigs = c(s_d1,s_d1))),
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(0,  0),sigs = c(s_d1,s_d1))))

  s_d1 <- 1
  D2 <- cbind(
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(.5,0),sigs = c(s_d1,s_d1))),
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(0, 0),sigs = c(s_d1,s_d1))))

  s_d1 <- 1
  D3 <- cbind(
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(.75,0),sigs = c(s_d1,s_d1))),
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(0,  0),sigs = c(s_d1,s_d1))))

  # low noise
  s_d1 <- 1
  D4 <- cbind(
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(1,0),sigs = c(s_d1,s_d1))),
    replicate(NCOL,make_mix(ns = c(NOBS,NOBS),mus = c(0,0),sigs = c(s_d1,s_d1))))

  TRUTH_BICL = matrix(0,NOBS*2,NCOL*2)
  TRUTH_BICL[1:NOBS,1:NCOL] = 1

  D0list[[i]] <- D0
  D1list[[i]] <- D1
  D2list[[i]] <- D2
  D3list[[i]] <- D3
  D4list[[i]] <- D4
  CClist[[i]] <- TRUTH_CC  <- rep(0:1,each=NCOL)
  RClist[[i]] <- TRUTH_BICL

}

saveRDS(list(
  D0list = D0list,
  D1list = D1list,
  D2list = D2list,
  D3list = D3list,
  D4list = D4list,
  CClist = CClist,
  RClist = RClist),"Simulation_Study_Poseidon_vs_Pose/Output/RDS/data_and_true_clustering.RDS")

cat(paste("\nDONE!"))
