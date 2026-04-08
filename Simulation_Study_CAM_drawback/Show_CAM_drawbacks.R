library(Poser)

# Hyperparameters ---------------------------------------------------------
List_cpp_NN <- NULL
hyperpar = list(
  a  = 1e-4,
  b  = 1e-4,
  m0 = 0,
  k0 = .01,
  c0 = 3,
  d0 = 2,
  s1 = 1,
  s2 = 1,
  s3 = 1,
  s4 = 1)
make_mix <- function(ns,mus,sigs){
  l <- as.list(c(sapply(1:length(ns), function(y) rnorm(ns[y],mus[y],sigs[y]))))
  do.call(c,l)
}
f_CAM <- function(i,Y){
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

  VV =   SEr_CAM(
    Y = Y,
    L = L,
    K = K,
    hyperpar = hyperpar,
    epsilon = 1e-8,
    RHO_init = RHOi,
    XI_init = XI_ikl,
    nsim = 500,verbose = F)
  VV
}
f_fiSAN <- function(i,Y){
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

  VV =   SEr_fiSAN(
    Y = Y,
    L = L,
    K = K,
    hyperpar = hyperpar,
    epsilon = 1e-8,
    RHO_init = RHOi,
    XI_init = XI_ikl,
    nsim = 500,verbose = F)
  VV
}



# 10 obs per column -------------------------------------------------------


set.seed(252525)
s_d1 <- .25
D1 <- cbind(
  replicate(10,make_mix(ns = c(5,5),mus = c(0,0),sigs = c(s_d1,s_d1))),
  #  replicate(4,make_mix(c(15,15),c(0,0),c(s_d1,s_d1))),
  replicate(10,make_mix(c(2,8),c(0,-3),c(s_d1,s_d1))),
  replicate(10,make_mix(c(2,8),c(0,3),c(s_d1,s_d1))))
saveRDS(D1, "Simulation_Study_CAM_drawback/Output/RDS/D_10obs.RDS")
image(D1)
K <- 20
L <- 15
res_fisan <- res_cam <- list()

for(t in 1:100){
  res_cam[[t]]   <- f_CAM(i = 2*t,D1)
  res_fisan[[t]] <- f_fiSAN(i = 2*t,D1)
  cat(t)
}
saveRDS(res_cam,"Simulation_Study_CAM_drawback/Output/RDS/run_CAM_10obs.RDS")
saveRDS(res_fisan,"Simulation_Study_CAM_drawback/Output/RDS/run_fiSAN_10obs.RDS")

best_cam <- extract_best_run_from_list(res_list = res_cam)
best_fisan <- extract_best_run_from_list(res_list = res_fisan)

# draw_dens_pose(best_cam,k = 1)
# draw_dens_pose(best_cam,k = 2)
# draw_dens_pose(best_cam,k = 3)
# image(best_cam$RHO)
# draw_dens_pose(best_fisan,hist_breaks =25,k =1)
# draw_dens_pose(best_fisan,hist_breaks =25,k =2)
# draw_dens_pose(best_fisan,hist_breaks =25,k =3)
# image(best_fisan$RHO)



# 30 obs per column -------------------------------------------------------------------------
set.seed(252525)
s_d1 <- .25
D1 <- cbind(
  replicate(10,make_mix(ns = c(15,15),mus = c(0,0),sigs = c(s_d1,s_d1))),
  #  replicate(4,make_mix(c(15,15),c(0,0),c(s_d1,s_d1))),
  replicate(10,make_mix(c(12,18),c(0,-3),c(s_d1,s_d1))),
  replicate(10,make_mix(c(12,18),c(0,3),c(s_d1,s_d1))))
saveRDS(D1, "Simulation_Study_CAM_drawback/Output/RDS/D_30obs.RDS")
image(D1)
K <- 20
L <- 15
res_fisan <- res_cam <- list()

for(t in 1:100){
  res_cam[[t]]   <- f_CAM(i = 2*t,D1)
  res_fisan[[t]] <- f_fiSAN(i = 2*t,D1)
  cat(t)
}
saveRDS(res_cam,"Simulation_Study_CAM_drawback/Output/RDS/run_CAM_30obs.RDS")
saveRDS(res_fisan,"Simulation_Study_CAM_drawback/Output/RDS/run_fiSAN_30obs.RDS")

best_cam <- extract_best_run_from_list(res_list = res_cam)
best_fisan <- extract_best_run_from_list(res_list = res_fisan)


# draw_dens_pose(best_cam,k = 1,hist_breaks =25)
# draw_dens_pose(best_cam,k = 2,hist_breaks =25)
# draw_dens_pose(best_cam,k = 3,hist_breaks =25)
# image(best_cam$RHO)
# draw_dens_pose(best_fisan,hist_breaks =25,k =1)
# draw_dens_pose(best_fisan,hist_breaks =25,k =2)
# draw_dens_pose(best_fisan,hist_breaks =25,k =3)
# image(best_fisan$RHO)
