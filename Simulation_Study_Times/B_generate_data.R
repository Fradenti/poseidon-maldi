source("Simulation_Study_Times/A_AUX.R")
# Generate data -----------------------------------------------------------

s_d1 <- 1
NCOL <- c(100,250,500)
NOBS <- c(25,50,100,250) # these get doubled
num_data <- 4
DATASETS <- matrix(list(),length(NCOL), length(NOBS))


for(i in 1:length(NCOL)){
  for(j in 1:length(NOBS)){
    data_list <- list()
    for(k in 1:num_data){

      data_list[[k]] <- cbind(replicate(NCOL[i],
                            make_mix(ns = c(NOBS[j],NOBS[j]),mus = c(4,0),sigs = c(1,1))))

    }
    DATASETS[[i,j]] <- data_list
  }
}


saveRDS(DATASETS,"Simulation_Study_Times/Output/RDS/All_datasets.RDS")

