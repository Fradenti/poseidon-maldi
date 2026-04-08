library(tidyverse)
library(raster)

source("Simulation_VaryingSpatial_competitors/A_Aux_Function_1.R")
source("Simulation_VaryingSpatial_competitors/A_Aux_Function_2.R")
source("Simulation_VaryingSpatial_competitors/B_Hyperpar.R")

# -------------------------------------------------------------------------

n_side_x <- 20
n_side_y <- 20
n_pix = n_pixel = n_side_x*n_side_y

# coords = expand.grid(x = 1:n_side, y = 1:n_side)

rast <- raster::raster(matrix(1,nrow = n_side_y,ncol = n_side_x))
coords <- as.data.frame(raster::rowColFromCell(rast, 1:raster::ncell(rast)))
coords$cell <- cellFromRowCol(
  rast,
  coords$row,
  coords$col
)
colnames(coords) <- c("y", "x", "cell")
coords <- coords[,c("cell", "x", "y")]


dist_dahl <- as.matrix(stats::dist(coords[,c("x", "y")], method = ifelse(ROOK_vs_QUEEN,"manhattan","maximum")))

## this code generates multiple (Nsample) partitions! one with their own permutation
Nsample <- 10

for(i_str in seq_along(strength)){
  set.seed(1234)
  perm_matrix <- sapply(1:Nsample, function (x) sample(n_pix)) # one permutation for each column (n x Nsample)
  str <- strength[i_str]
  alpha_DP <- 2; lambda_dahl <- exp(-dist_dahl*str)
  tmp <- apply(perm_matrix, 2, function(perm) sample_dahl(perm,lambda_dahl,alpha_DP) )
  parts <- tmp[1:n_pix,] # if you need the log prior probabilities, they're stored in tmp[n_pix+1,]

  for(i_sim in 1:Nsample){
    cl_pix <- parts[,i_sim]
    K_pix <- length(unique(cl_pix))

    # coords %>% mutate(cluster = cl_pix) %>%
    #   ggplot() + geom_tile(aes(x=x,y=y,fill = as.factor(cluster)), color = "gray") +
    #   scale_fill_brewer(palette = "Set3", name = "Cluster") + xlab("")+ ylab("") +
    #   theme_minimal() + theme(axis.text = element_blank(),
    #                           axis.ticks = element_blank())

    # -------------------------------------------------------------------------

    ### ONE DATASET
    set.seed(20000 + i_sim)

    m_peak <- 1.5
    n_rows <- 500
    data <- matrix(NA, nrow = n_rows, ncol = n_pix)
    true_bicl <- matrix("noise", nrow = n_rows, ncol = n_pix)
    true_mean <- matrix(0, nrow = n_rows, ncol = n_pix)
    n_peaks <- round(n_rows * 0.02/2)
    mult <- 1
    sd <- 1

    starting_point <- list()
    starting_point[[1]] <- cl_pix
    starting_point[[2]] <- list()

    peaks_list <- list()

    for(k in 1:K_pix){
      cl <- which(cl_pix == k)
      z_rows <- rep("noise",n_rows)
      # each cluster has a different pattern, except for the first one (which is just noise, no peaks)
      if(k > 1){
        peaks <- sort(sample(which((1:n_rows)%%50==5), size = n_peaks))
        peaks_list[[k-1]] <- peaks

        while(all.equal(length(peaks_list), length(unique(peaks_list))) != TRUE){
          peaks <- sort(sample(which((1:n_rows)%%50==5), size = n_peaks))
          peaks_list[[k-1]] <- peaks
        }

        med <- c(peaks-1, peaks+2)
        peaks <- c(peaks, peaks + 1)
        noise <- setdiff(1:n_rows, c(peaks,med))
        z_rows[peaks] <- "peaks"
        z_rows[med] <- "med"

        true_bicl[z_rows == "peaks",cl] <- "peaks"
        true_bicl[z_rows == "med",cl] <- "med"

        true_mean[z_rows == "peaks",cl] <- m_peak/mult
        true_mean[z_rows == "med",cl] <- m_peak/mult/2
      }
      z_rows_num <- as.numeric(factor(z_rows, levels = c("noise", "peaks", "med")))
      starting_point[[2]][[k]] <- z_rows_num

      data[,cl] <- true_mean[,cl]
      data[,cl] <- data[,cl] + rnorm(n_rows*length(cl), sd = sd/mult)
    }

    if(all.equal(length(peaks_list), length(unique(peaks_list))) != TRUE){
      stop("-------- peaks are not unique across clusters! --------\n")
    }



    ### export
    paste0("pix_",coords$cell) -> colnames(data)
    paste0("mz_",1:n_rows) -> rownames(data)

    X   <- as.matrix(t(data))     # potentially replace with Matrix::Matrix(data, sparse = TRUE)
    var <- data.frame(feature = colnames(X))
    obs <- data.frame(
      cell_id = colnames(data),
      x = as.integer(coords$x),
      y = as.integer(coords$y),
      true_cluster = cl_pix,
      sample_id = paste0("sim_", i_sim)
    )

    write.csv(cbind(obs,X),
              paste0("Simulation_VaryingSpatial_competitors/data/simulated_dataC/simulated_data_str",str,"_",i_sim,".csv"),
              row.names = FALSE)

    biclust_data <- t(rbind(true_mean,true_bicl))
    write.csv(biclust_data,
              paste0("Simulation_VaryingSpatial_competitors/data/simulated_dataC/simulated_data_bicl_str",str,"_",i_sim,".csv"),
              row.names = FALSE)
  }
}



