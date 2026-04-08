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
           VI = NA) -> results_summary

for(i_str in seq_along(strength)){
  str = strength[i_str]
  for(i_sim in 1:Nsample){
    results <- load(paste0("Simulation_VaryingSpatial_competitors/results/RData/BayesSpace/BayesSpace_result_str",str,"_",i_sim,".RData"))

    index <- which((results_summary$i_sim == i_sim) & (results_summary$i_str == i_str))
    results_summary[index,"str"] <- str

    results_summary[index,"ARI"] <- ARI_colcl
    results_summary[index,"VI"] <- VI_colcl
  }
}

plot(results_summary$i_str, results_summary$ARI)
plot(results_summary$i_str, results_summary$VI)

save(list = c("results_summary"),
     file = "Simulation_VaryingSpatial_competitors/results/BayesSpace_results_summary.RData")
