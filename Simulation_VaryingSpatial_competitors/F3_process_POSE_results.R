library(tidyverse)
library(mcclust)

source("Simulation_VaryingSpatial_competitors/A_Aux_Function_1.R")
source("Simulation_VaryingSpatial_competitors/A_Aux_Function_2.R")
source("Simulation_VaryingSpatial_competitors/B_Hyperpar.R")

Nsample <- 40

data.frame(i_sim = rep(1:Nsample, length(strength)),
           i_str = rep(seq_along(strength),each = Nsample),
           str = NA,
           ARI = NA,
           VI = NA,
           ARI_bicl = NA,
           VI_bicl = NA,
           RMSE_bicl = NA) -> results_summary

for(i_str in seq_along(strength)){
  str = strength[i_str]
  for(i_sim in 1:Nsample){
    results <- load(paste0("Simulation_VaryingSpatial_competitors/results/RData/Pose/POSE_result_str",str,"_",i_sim,".RData"))

    index <- which((results_summary$i_sim == i_sim) & (results_summary$i_str == i_str))
    results_summary[index,"str"] <- str

    results_summary[index,"ARI"] <- ARI_colcl
    results_summary[index,"VI"] <- VI_colcl
    results_summary[index,"ARI_bicl"] <- ARI_bicl
    results_summary[index,"VI_bicl"] <- VI_bicl
    results_summary[index,"RMSE_bicl"] <- RMSE_bicl
  }
}

plot(results_summary$i_str, results_summary$ARI)
plot(results_summary$i_str, results_summary$VI)

plot(results_summary$i_str, results_summary$ARI_bicl)
plot(results_summary$i_str, results_summary$VI_bicl)
plot(results_summary$i_str, results_summary$RMSE_bicl)


save(list = c("results_summary"),
     file = "Simulation_VaryingSpatial_competitors/results/POSE_results_summary.RData")
