library(tidyverse)
library(mcclust)

source("Simulation_VaryingSpatial_competitors/A_Aux_Function_1.R")
source("Simulation_VaryingSpatial_competitors/A_Aux_Function_2.R")
source("Simulation_VaryingSpatial_competitors/B_Hyperpar.R")

Nsample <- 40

data.frame(i_sim = rep(1:Nsample, length(strength)),
           i_str = rep(seq_along(strength),each = Nsample),
           str = NA,
           ARI1 = NA, ARI2 = NA, ARI3 = NA,ARI4 = NA,
           VI1 = NA, VI2 = NA, VI3 = NA,VI4 = NA ) -> results_summary

for(i_str in seq_along(strength)){
  str = strength[i_str]
  for(i_sim in 1:Nsample){
    results <- read.csv(paste0("Simulation_VaryingSpatial_competitors/results/RData/UTAG/utag_results_str",str,"_",i_sim,".csv"))

    index <- which((results_summary$i_sim == i_sim) & (results_summary$i_str == i_str))
    results_summary[index,"str"] <- str

    z1 <- results$utag_cluster_02
    z2 <- results$utag_cluster_03
    z3 <- results$utag_cluster_04
    z4 <- results$utag_cluster_05


    results_summary[index,"ARI1"] <- mcclust::arandi(results$true_cluster, z1)
    results_summary[index,"ARI2"] <- mcclust::arandi(results$true_cluster, z2)
    results_summary[index,"ARI3"] <- mcclust::arandi(results$true_cluster, z3)
    results_summary[index,"ARI4"] <- mcclust::arandi(results$true_cluster, z4)
    results_summary[index,"VI1"] <- mcclust::vi.dist(results$true_cluster, z1)
    results_summary[index,"VI2"] <- mcclust::vi.dist(results$true_cluster, z2)
    results_summary[index,"VI3"] <- mcclust::vi.dist(results$true_cluster, z3)
    results_summary[index,"VI4"] <- mcclust::vi.dist(results$true_cluster, z4)
  }
}

# boxplot(results_summary$ARI1,results_summary$ARI2,results_summary$ARI3,results_summary$ARI4)
# boxplot(results_summary$VI1,results_summary$VI2,results_summary$VI3,results_summary$VI4)
#
# results_summary$best_ARI <- apply(results_summary[,3+1:4],1,max)
# results_summary$best_VI <- apply(results_summary[,3+4+1:4],1,min)
#
# plot(results_summary$i_str, results_summary$best_ARI)
# plot(results_summary$i_str, results_summary$best_VI)

save(list = c("results_summary", "results"), # note: we pass one of the results objects to see what resolutions parameters it has been constructed with
     file = "Simulation_VaryingSpatial_competitors/results/UTAG_results_summary_str.RData")
