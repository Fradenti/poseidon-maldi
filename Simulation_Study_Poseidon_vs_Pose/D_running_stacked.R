source("Simulation_Study_Poseidon_vs_Pose/A_AUX.R")
source("Simulation_Study_Poseidon_vs_Pose/C_Hyperpar.R")

DD <- readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/data_and_true_clustering.RDS")

RESULTS_stacked <- list()
for(qq in 1:4){
  l = list()
  for(j in 1:100){ # replications
    Res = list()
    for(i in 1:50){   # starting points
      Res[[i]] <- f(i,
                    Y = rbind(DD[[1]][[j]],DD[[qq+1]][[j]])
      )
    }

    indmax = which.max(unlist(lapply(Res, function(x)max(x$Elbo_val) )))
    l[[j]] <- Res[[indmax]]
    cat(paste("done with iteration",j, "of dataset",qq, "\n"))
  }
  RESULTS_stacked[[qq]] <- l
}
saveRDS(RESULTS_stacked,"Simulation_Study_Poseidon_vs_Pose/Output/RDS/posers_stacked_dec25.RDS")
