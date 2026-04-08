library(Cardinal)
library(superheat)
library(POSEIDON)
library(Poser)
library(sparseBC)
library(cluster)
library(tidyverse)
source("Simulation_Biclustering_recovery/A_Aux_Function_1.R")
source("Simulation_Biclustering_recovery/A_Aux_Function_2.R")
source("Simulation_Biclustering_recovery/A_Aux_Function_3.R")
source("Simulation_Biclustering_recovery/B_Hyperpar.R")

# -------------------------------------------------------------------------
nrep  <- 50
n_sim <- 50

RES_ari_colcl <- matrix(NA,n_sim,4)
RES_ari_bicl  <- matrix(NA,n_sim,4)
RES_rmse_bicl <- matrix(NA,n_sim,4)

colnames(RES_ari_colcl) <- c("Poser", "sparseBC\n(true K)", "sparseBC\n(choose KR)", "Double\nk-means")
colnames(RES_ari_bicl) <- c("Poser", "sparseBC\n(true K)", "sparseBC\n(choose KR)",  "Double\nk-means")
colnames(RES_rmse_bicl) <- c("Poser", "sparseBC\n(true K)", "sparseBC\n(choose KR)", "Double\nk-means")
for(i_sim in 1:n_sim){
  for(i_met in 1:4){
    res <-  readRDS(paste0("Simulation_Biclustering_recovery/OUTPUT/RDS/simuC_met",i_met,"_sim",i_sim,".RDS"))

    RES_ari_colcl[i_sim,i_met] <- res$ARI_colcl
    RES_rmse_bicl[i_sim,i_met] <- res$RMSE_bicl
    RES_ari_bicl[i_sim,i_met]  <- res$ARI_bicl
    }
}


RES_ari_colcl <- as_tibble(RES_ari_colcl) %>% mutate(type = "ARI (column clustering)")
RES_ari_bicl <- as_tibble(RES_ari_bicl) %>% mutate(type =   "ARI (biclustering)")
RES_rmse_bicl <- as_tibble(RES_rmse_bicl) %>% mutate(type = "RMSE (biclustering)")

RES <- RES_ari_bicl %>% bind_rows(RES_ari_colcl,RES_rmse_bicl) %>% reshape2::melt() %>%
  mutate(type = factor(type,levels = c("ARI (column clustering)","ARI (biclustering)", "RMSE (biclustering)")))

RES
q <- ggplot(RES)+
  theme_bw()+
  ylab("")+xlab("Method")+
  geom_boxplot(aes(x=variable, y=value))+
  facet_wrap(~type,scale = "free_y")+
  theme(text = element_text(size=16))


q+ggview::canvas(h=5,w=15)


ggsave("Simulation_Biclustering_recovery/results_bicl_recovery.pdf",plot = q,width = 15,height = 5)
ggsave("Simulation_Biclustering_recovery/results_bicl_recovery.png",plot = q,width = 15,height = 5)
ggsave("Simulation_Biclustering_recovery/results_bicl_recovery.eps",plot = q,width = 15,height = 5)
