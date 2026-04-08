source("Simulation_Study_Poseidon/A_AUX.R")
source("Simulation_Study_Poseidon/C_Hyperpar.R")

s_d1 <- 1
NCOL <- c(100,250,500)
NOBS <- c(25,50,100,250) # these get doubled
num_data <- 4
DATASETS <- DD <- readRDS("Simulation_Study_Times/Output/RDS/All_datasets.RDS")

RESULTS <-  matrix(list(),length(NCOL), length(NOBS))

for(k in 1:num_data){
l = list()
  for(i in 1:length(NCOL)){
  for(j in 1:length(NOBS)){
  Res <-  list()
  for(s in 1:50){ # starting points
   Res[[s]] <- g_time_zeroepsilon_10iterations(i = i,
        YY = DATASETS[[i,j]][1:k]
      )
  cat(paste("starting point:", s, "-- J:", i,"-- N:",j,"-- T:",k,"\n"))
   }
  times <- unlist(lapply(Res, function(x) as.numeric(x$time, units = "secs")))
  RESULTS[[i,j]][[k]] <- times
  }
  }
}
saveRDS(RESULTS,"Simulation_Study_Times/Output/RDS/Times.RDS")


RESULTS <-  matrix(list(),length(NCOL), length(NOBS))

for(k in 1:num_data){
l = list()
  for(i in 1:length(NCOL)){
  for(j in 1:length(NOBS)){
  Res <-  list()
  for(s in 1:50){ # starting points
   Res[[s]] <- g_time_zeroepsilon_10iterations(i = i,
        YY = DATASETS[[i,j]][1:k]
      )
  cat(paste("starting point:", s, "-- J:", i,"-- N:",j,"-- T:",k,"\n"))
   }
  times <- unlist(lapply(Res, function(x) as.numeric(x$time, units = "secs")))
  RESULTS[[i,j]][[k]] <- times
  }
  }
}
saveRDS(RESULTS,"Simulation_Study_Times/Output/RDS/Times.RDS")





