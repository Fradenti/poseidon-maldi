source("Simulation_Study_Poseidon/A_AUX.R")
source("Simulation_Study_Poseidon/C_Hyperpar.R")

s_d1 <- 1
NCOL <- c(100,250,500)
NOBS <- c(25,50,100,250) # these get doubled
num_data <- 4
DATASETS <- DD <- readRDS("Simulation_Study_Times/Output/RDS/All_datasets.RDS")

RESULTS <-  matrix(list(),length(NCOL), length(NOBS))

k <- 1
l <- list()
for(i in 1:length(NCOL)){
  for(j in 1:length(NOBS)){
    Res <-  list()
      for(s in 1:50){ # starting points
        Res[[s]] <- f_cam(i = i,
                          Y = DATASETS[[i,j]][1:k][[1]]
        )
        cat(paste("starting point:", s, "-- J:", i,"-- N:",j,"-- T:",k,"\n"))
      }
      times <- unlist(lapply(Res, function(x) as.numeric(x$time, units = "secs")))
      RESULTS[[i,j]][[k]] <- times
    }
}

saveRDS(RESULTS,"Simulation_Study_Times/Output/RDS/Times_singledata_cam.RDS")


RESULTS <-  matrix(list(),length(NCOL), length(NOBS))

k <- 1
l <- list()
for(i in 1:length(NCOL)){
  for(j in 1:length(NOBS)){
    Res <-  list()
    for(s in 1:50){ # starting points
      Res[[s]] <- fisan(i = i,
                        Y = DATASETS[[i,j]][1:k][[1]]
      )
      cat(paste("starting point:", s, "-- J:", i,"-- N:",j,"-- T:",k,"\n"))
    }
    times <- unlist(lapply(Res, function(x) as.numeric(x$time, units = "secs")))
    RESULTS[[i,j]][[k]] <- times
  }
}

saveRDS(RESULTS,"Simulation_Study_Times/Output/RDS/Times_singledata_fisan.RDS")
