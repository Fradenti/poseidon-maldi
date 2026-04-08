source("Simulation_VaryingSpatial_competitors/A_Aux_Function_1.R")
source("Simulation_VaryingSpatial_competitors/A_Aux_Function_2.R")
source("Simulation_VaryingSpatial_competitors/B_Hyperpar.R")

library(Poser)
library(mcclust)
library(raster)

# -------------------------------------------------------------------------
# Simulation Script with i_sim input and result saving
# -------------------------------------------------------------------------
# i_sim <- 1

# Accept i_sim from environment (or default to 1)
# i_sim <- as.integer(Sys.getenv("i_sim", unset = "1"))
# i_str <- as.integer(Sys.getenv("i_str", unset = "1"))
# cat("Running simulation number:", i_sim, " with i_str ", i_str,"\n")

for(i_str in 1:4){
str <- strength[i_str]

for(i_sim in 1:40){

  full_data = read.csv(file = paste0("Simulation_VaryingSpatial_competitors/data/simulated_dataC/simulated_data_str",str,"_",i_sim,".csv"))
  metadata_cols = c('cell_id', 'x', 'y', 'true_cluster', 'sample_id')
  expression_cols = setdiff(colnames(full_data),metadata_cols)
  Xt = apply(full_data[,expression_cols], 2, as.numeric)
  X <- t(Xt)
  obs = full_data[,metadata_cols]
  obs$true_cluster <- as.numeric(obs$true_cluster)

  X1 <- mapper1(X)

  # this is what we'd have here
  n_side_x <- max(obs$x)
  n_side_y <- max(obs$y)
  rast <- raster::raster(matrix(1,nrow = n_side_y,ncol = n_side_x))
  ADJ  <- raster::adjacent(rast,
                           cells = as.numeric(substring(obs$cell,5)),
                           directions=ifelse(ROOK_vs_QUEEN,4,8),sorted=T)
  nnList = list()
  for(i in 1:raster::ncell(rast)){
    nnList[[i]] = ADJ[which(ADJ[,1]==i),2]
  }
  nnList_cpp = lapply(nnList,function(x)x-1)


  bicl_data = read.csv(file = paste0("Simulation_VaryingSpatial_competitors/data/simulated_dataC/simulated_data_bicl_str",str,"_",i_sim,".csv"))
  mean_cols = 1:(ncol(bicl_data)/2)
  truebicl_cols = (ncol(bicl_data)/2) + 1:(ncol(bicl_data)/2)
  true_mean = apply(bicl_data[,mean_cols], 2, as.numeric)
  true_bicl = apply(bicl_data[,truebicl_cols], 2, as.factor)
  true_mean = t(true_mean)
  true_bicl = t(true_bicl)

  tm <- Sys.time()

  Res <- list()
  for(i in 1:nrep){
    cat(i, " ")
    Res[[i]] <- f_sp(i = i,Y = X1,K = K,L = L,hyperpar = hyperpar1)
  }
  cat("\n")
  indmax = which.max(unlist(lapply(Res, function(x)max(x$Elbo_val) )))
  res1 <- Res[[indmax]]

  colcl1 = apply(res1$RHO,1,which.max)
  tmp = Poser::estimate_clust_pose(res1)

  est_colcl <- colcl1
  est_bicl <- tmp$bicl
  ARI_colcl <- mcclust::arandi(obs$true_cluster,est_colcl)
  VI_colcl <- mcclust::vi.dist(obs$true_cluster,est_colcl)
  RMSE_bicl <- sqrt(mean((mapper1_ext(true_mean,X) - est_bicl)^2))
  ARI_bicl <- mcclust::arandi(as.factor(c(est_bicl)), as.factor(c(true_bicl)))
  VI_bicl <- mcclust::vi.dist(as.factor(c(est_bicl)), as.factor(c(true_bicl)))


  tot_time <- Sys.time() - tm

  save_list <- c(save_list, "res1",
                 "ARI_colcl", "VI_colcl","RMSE_bicl", "ARI_bicl", "VI_bicl",
                 "est_colcl", "est_bicl", "tot_time")
  save_file <- paste0("Simulation_VaryingSpatial_competitors/results/RData/Pose/POSE_result_str",str,"_", i_sim, ".RData")
  save(list = save_list, file = save_file)
}

}


