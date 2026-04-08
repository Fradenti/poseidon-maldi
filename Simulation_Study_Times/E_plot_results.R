library(tidyverse)
theme_set(theme_bw())

s_d1 <- 1
NCOL <- c(100,250,500)
NOBS <- c(25,50,100,250) # these get doubled
num_data <- 4
DATASETS <- DD <- readRDS("Simulation_Study_Times/Output/RDS/All_datasets.RDS")
RESULTS <- readRDS("Simulation_Study_Times/Output/RDS/Times.RDS")


d <- c()
for(i in 1:length(NCOL)){ # replications
    for(j in 1:length(NOBS)){
      q = cbind(J = NCOL[i], N = NOBS[j]*2, data = rep(1:4,each=50) , time =  unlist(RESULTS[[i,j]]))
      d = rbind(d,q)
    }
  }


d %>%
  as_tibble() %>%
  mutate(data = paste("T =",data), J = paste("J =",J)) %>%
  ggplot(aes(x = factor(N), y = time/10)) +
  geom_boxplot(
    aes(
#      col  = factor(V3),
#      group = interaction(V2, V3)
    )
  ) +
  facet_grid(data ~ J)
#  scale_y_log10()


Q <- d %>%
  as_tibble() %>%   mutate(data = paste("T =",data), J = paste("J =",J),
                           N = N) %>%
  ggplot(aes(x = factor(N), y = time/10)) +
  geom_boxplot(
    aes(
      #      col  = factor(V3),
      #      group = interaction(V2, V3)
    )
  ) +
  facet_grid(data ~ J)+
 scale_y_log10()+
  theme(text = element_text(size=15))+
  ylab("Average time per iteration in seconds (log10 scale)")+xlab("N")

Q+ggview::canvas(h=7,w=7)
ggsave("Simulation_Study_Times/Output/PLOT/boxplots_seconds.png",h=7,w=7)
ggsave("Simulation_Study_Times/Output/PLOT/boxplots_seconds.pdf",h=7,w=7)
ggsave("Simulation_Study_Times/Output/PLOT/boxplots_seconds.eps",h=7,w=7)
