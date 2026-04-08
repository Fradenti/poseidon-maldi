library(tidyverse)
theme_set(theme_bw())

s_d1 <- 1
NCOL <- c(100,250,500)
NOBS <- c(25,50,100,250) # these get doubled
num_data <- 4
DATASETS <- DD <- readRDS("Simulation_Study_Times/Output/RDS/All_datasets.RDS")
Rescam = readRDS("Simulation_Study_Times/Output/RDS/Times_singledata_cam.RDS")
Resfis = readRDS("Simulation_Study_Times/Output/RDS/Times_singledata_fisan.RDS")



d <- c()
for(i in 1:length(NCOL)){ # replications
  for(j in 1:length(NOBS)){
    q = cbind(J = NCOL[i], N = NOBS[j]*2, time =  unlist(Rescam[[i,j]]))
    d = rbind(d,q)
  }
}
d_cam <- d %>% as.data.frame() %>% mutate(model="Common Atoms")
d <- c()
for(i in 1:length(NCOL)){ # replications
  for(j in 1:length(NOBS)){
    q = cbind(J = NCOL[i], N = NOBS[j]*2,  time =  unlist(Resfis[[i,j]]))
    d = rbind(d,q)
  }
}
d_fis <- d %>% as.data.frame() %>% mutate(model="Shared Atoms")


d_all <- rbind(d_cam,d_fis)



Q <- d_all %>%
  as_tibble() %>%   mutate(J = paste("J =",J),
                           N = N) %>%
  ggplot(aes(x = factor(N), y = time/10,col=model,group-model)) +
  geom_boxplot(
    aes(
      #      col  = factor(V3),
      #      group = interaction(V2, V3)
    )
  ) +
  facet_wrap( ~ J)+
  scale_y_log10()+
  scale_color_manual(values = c("gray","black"))+
  theme(text = element_text(size=15),legend.position = "bottom")+
  ylab("Average time per iteration\nin seconds (log10 scale)")+xlab("N")

Q+ggview::canvas(h=5,w=10)

ggsave("Simulation_Study_Times/Output/PLOT/boxplots_seconds_single.png",h=5,w=10)
ggsave("Simulation_Study_Times/Output/PLOT/boxplots_seconds_single.pdf",h=5,w=10)
ggsave("Simulation_Study_Times/Output/PLOT/boxplots_seconds_single.eps",h=5,w=10)
