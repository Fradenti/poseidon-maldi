# same as new_semireal_simu but trying to have more similar peaks
library(Cardinal)
library(superheat)
library(POSEIDON)
library(Poser)
library(sparseBC)
library(cluster)
source("Simulation_Biclustering_recovery/A_Aux_Function_1.R")
source("Simulation_Biclustering_recovery/A_Aux_Function_2.R")
source("Simulation_Biclustering_recovery/B_Hyperpar.R")


DATA <- list()
# -------------------------------------------------------------------------
# - Generete Data
# -------------------------------------------------------------------------
for(i_sim in 1:50){

    set.seed(10000 + i_sim)
    n_side <- 20
    n_pix <- n_pixel <-  n_side^2

    coords <- expand.grid(x = 1:n_side, y = 1:n_side)
    coords_extended <- cbind(coords,
                            index = 1:n_pixel,
                            square1 = (coords$x-1) %/% 5 + 1,
                            square2 = (coords$y-1) %/% 5 + 1,
                            square3 = (coords$x-1) %/% 10 + 1,
                            square4 = (coords$y-1) %/% 10 + 1)

    coords_extended$cluster <-
      as.factor(paste0(coords_extended$square3, "_", coords_extended$square4))
    coords_extended$cluster_num <- as.numeric(coords_extended$cluster)
    coords_extended$cluster_tmp <-
      as.factor(paste0(coords_extended$square1, "_", coords_extended$square2))
    coords_extended$cluster <- ifelse(coords_extended$cluster == "1_1",
                                      coords_extended$cluster_tmp, coords_extended$cluster)

    coords_extended <- coords_extended[,c("x","y","index","cluster")]


    K_pix <- length(unique(coords_extended$cluster))
    cl_pix <- coords_extended$cluster

    rast <- raster::raster(matrix(1,n_side,n_side))
    ADJ  <- raster::adjacent(rast,
                             cells = 1:raster::ncell(rast),
                             directions=8,sorted=T)
    nnList = list()
    for(i in 1:raster::ncell(rast)){
      nnList[[i]] = ADJ[which(ADJ[,1]==i),2]
    }
    nnList_cpp = lapply(nnList,function(x)x-1)

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

    for(k in 1:K_pix){
      cl <- which(cl_pix == k)
      z_rows <- rep("noise",n_rows)
      # each cluster has a different pattern, except for the first one
      if(k > 1){
        peaks <- sort(sample(which((1:n_rows)%%50==5), size = n_peaks))
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

    data1      <- data
    true_bicl1 <- true_bicl
    X1         <- mapper1(data1)
    DATA   <- list(X1 = X1,
                   data1 = data1,
                   nnList_cpp = nnList_cpp,
                   true_mean = true_mean,
                   true_bicl = true_bicl,
                   cl_pix = cl_pix)
  saveRDS(DATA, paste0("Simulation_Biclustering_recovery/OUTPUT/RDS/DATASETS/data_bicl_",i_sim,".RDS"))
}
