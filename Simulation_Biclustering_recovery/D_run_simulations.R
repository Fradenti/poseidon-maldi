library(Cardinal)
library(superheat)
library(POSEIDON)
library(Poser)
library(sparseBC)
library(cluster)
source("Simulation_Biclustering_recovery/A_Aux_Function_1.R")
source("Simulation_Biclustering_recovery/A_Aux_Function_2.R")
source("Simulation_Biclustering_recovery/A_Aux_Function_3.R")
source("Simulation_Biclustering_recovery/B_Hyperpar.R")

# -------------------------------------------------------------------------
nrep  <- 50
n_sim <- 50

for(i_sim in 1:n_sim){
  tm <- Sys.time()
  X_data <- readRDS(paste0("Simulation_Biclustering_recovery/OUTPUT/RDS/DATASETS/data_bicl_",i_sim,".RDS"))
  X1 <- X_data$X1
  nnList_cpp <- X_data$nnList_cpp
  cl_pix <- X_data$cl_pix
  data1 <- X_data$data1
  true_mean <- X_data$true_mean
  true_bicl <- X_data$true_bicl

  for(i_met in 1:4){

    if(i_met == 1){ # POSER
    cat("--- Poser ---\n")
    # -------------------------------------------------------------------------
   t1 <- Sys.time()
   Res <- list()
   for(i in 1:nrep){
     cat(i, " ")
     Res[[i]] <- f(i,X1,K,L,hyperpar1)
   }
   cat("\n")
   indmax = which.max(unlist(lapply(Res, function(x)max(x$Elbo_val) )))
   res1 <- Res[[indmax]]

   colcl1 <- apply(res1$RHO,1,which.max)

   tmp <- Poser::estimate_clust_pose(res1)

   est_colcl <- colcl1
   est_bicl <- tmp$bicl
   ARI_colcl <- salso::ARI(cl_pix,colcl1)
   RMSE_bicl <- sqrt(mean((mapper1_ext(true_mean,data1) - tmp$bicl)^2))
   ARI_bicl <- salso::ARI(as.factor(c(tmp$bicl)), as.factor(c(true_bicl)))

  } else if(i_met == 2){ # sparseBC with K_true

      cat("--- sparseBC ---\n")
      t1 <- Sys.time()

       # -------------------------------------------------------------------------

       K_true <- length(unique(cl_pix))
       out_sparseBC <- sparseBC(x = t(X1),k = K_true,r = L_upper, lambda = 10)

       ## postprocessing to get simpler biclustering
       set.seed(123*i_sim)
       data_tmp = as.matrix(c(out_sparseBC$Mus),ncol=1)
       gap_stat <- clusGap(data_tmp, FUN = kmeans, nstart = 25, K.max = 25, B = 50)
       plot(gap_stat)
       optimal_k <- with(gap_stat, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method="firstmax"))
       print(optimal_k)
       final_kmeans <- kmeans(data_tmp, centers = optimal_k, nstart = 25)

       Clustered_mus_matrix <- out_sparseBC$mus
       for(i in seq_along(c(out_sparseBC$Mus))){
         mu = c(out_sparseBC$Mus)[i]
         tmp_ind <- which(out_sparseBC$mus == mu, arr.ind = TRUE)
         Clustered_mus_matrix[tmp_ind] <- final_kmeans$centers[final_kmeans$cluster[i]]
       }
       Clustered_mus_matrix = t(Clustered_mus_matrix)

       est_colcl <- out_sparseBC$Cs
       est_bicl <- Clustered_mus_matrix
       ARI_colcl <- salso::ARI(cl_pix,out_sparseBC$Cs)
       RMSE_bicl <- sqrt(mean((mapper1_ext(true_mean,data1) - Clustered_mus_matrix)^2))
       ARI_bicl <- salso::ARI(as.factor(c(Clustered_mus_matrix)), as.factor(c(true_bicl)))

      } else if(i_met == 3){ # sparseBC.choosekr

        cat("--- sparseBC chooserkr ---\n")
        t1 <- Sys.time()

      # -------------------------------------------------------------------------

      tmp <- sparseBC.choosekr_patched(x = t(X1),
                               k = seq(2,K_upper,by=2),
                               r = Ls_choose,
                               lambda = 10,
                               percent = .1,
                               trace= TRUE)
      k <- tmp$estimated_kr[1,1]; k
      l <- tmp$estimated_kr[1,2]; l
      out_sparseBC <- sparseBC(x = t(X1),k = k,r = l, lambda = 10)

      ## postprocessing to get simpler biclustering
      set.seed(123*i_sim)
      data_tmp = as.matrix(c(out_sparseBC$Mus),ncol=1)
      gap_stat <- clusGap(data_tmp, FUN = kmeans, nstart = 25, K.max = L, B = 50)
      plot(gap_stat)
      optimal_k <- with(gap_stat, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method="firstmax"))
      print(optimal_k)
      final_kmeans <- kmeans(data_tmp, centers = optimal_k, nstart = 25)

      Clustered_mus_matrix <- out_sparseBC$mus
      for(i in seq_along(c(out_sparseBC$Mus))){
        mu = c(out_sparseBC$Mus)[i]
        tmp_ind <- which(out_sparseBC$mus == mu, arr.ind = TRUE)
        Clustered_mus_matrix[tmp_ind] <- final_kmeans$centers[final_kmeans$cluster[i]]
      }
      Clustered_mus_matrix = t(Clustered_mus_matrix)

      est_colcl <- out_sparseBC$Cs
      est_bicl <- Clustered_mus_matrix
      ARI_colcl <- salso::ARI(cl_pix,out_sparseBC$Cs)
      RMSE_bicl <- sqrt(mean((mapper1_ext(true_mean,data1) - Clustered_mus_matrix)^2))
      ARI_bicl <- salso::ARI(as.factor(c(Clustered_mus_matrix)), as.factor(c(true_bicl)))

    } else if(i_met == 4){ # double k-means (smart)
      cat("--- double k-means ---\n")
      # -------------------------------------------------------------------------
      t1 <- Sys.time()

      gap_stat <- clusGap(t(X1), FUN = kmeans, nstart = 25, K.max = K_upper, B = 50)
      optimal_k <- with(gap_stat, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method="firstmax"))
      final_kmeans <- kmeans(t(X1), centers = optimal_k, nstart = 25)

      col_avg <- c()
      for(k in seq_len(optimal_k)) {
        tmp_ind <- which(final_kmeans$cluster == k, arr.ind = TRUE)
        col_avg <- c(col_avg, rowMeans(X1[,tmp_ind,drop=FALSE]))
      }

      data_tmp <- as.matrix(col_avg,ncol=1)
      gap_stat_rows <- clusGap(data_tmp, FUN = kmeans, nstart = 25, K.max = L, B = 50)
      optimal_l <- with(gap_stat_rows, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method="firstmax"))
      final_kmeans_rows <- kmeans(data_tmp, centers = optimal_l, nstart = 25)

      bicl_mu_matrix <- X1*0
      for(k in seq_len(optimal_k)) {
        tmp_ind <- which(final_kmeans$cluster == k, arr.ind = TRUE)
        row_index <- (k-1)*nrow(X1) + 1:nrow(X1)
        bicl_mu_matrix[,tmp_ind] <- final_kmeans_rows$centers[final_kmeans_rows$cluster[row_index]]
      }

      est_colcl <- final_kmeans$cluster
      est_bicl <- bicl_mu_matrix
      ARI_colcl <- salso::ARI(cl_pix,final_kmeans$cluster)
      RMSE_bicl <- sqrt(mean((mapper1_ext(true_mean,data1) - bicl_mu_matrix)^2))
      ARI_bicl <- salso::ARI(as.factor(c(bicl_mu_matrix)), as.factor(c(true_bicl)))
    }

  tot_time <- Sys.time() - t1

  save_list <- list("ARI_colcl" = ARI_colcl,
                    "RMSE_bicl" = RMSE_bicl,
                    "ARI_bicl" = ARI_bicl,
                    "est_colcl" = est_colcl,
                    "est_bicl" = est_bicl,
                    "tot_time"=tot_time)

  save_file <- paste0("Simulation_Biclustering_recovery/OUTPUT/RDS/simuC_met", i_met,"_sim", i_sim, ".RDS")
  saveRDS(save_list, file = save_file)
  }
}
# -------------------------------------------------------------------------
