library(tidyverse)



tmp = load("Simulation_VaryingSpatial_competitors/results/POSE_results_summary.RData")
results_summary_POSE <- results_summary

tmp = load("Simulation_VaryingSpatial_competitors/results/BayesSpace_results_summary.RData")
results_summary_BS <- results_summary

tmp = load("Simulation_VaryingSpatial_competitors/results/UTAG_results_summary_str.RData")
results_summary_UTAG <- results_summary

combined <- rbind(results_summary_POSE %>%
                    mutate(method = "Pose") %>%
                    pivot_longer(cols = c(ARI, VI, ARI_bicl, VI_bicl, RMSE_bicl),
                                 names_to = "quantity",
                                 values_to = "value"),

                  results_summary_BS %>%
                    mutate(method = "BayesSpace") %>%
                    pivot_longer(cols = c(ARI, VI),
                                 names_to = "quantity",
                                 values_to = "value"),

                  results_summary_UTAG %>%
                    mutate(method = "UTAG") %>%
                    select(i_sim, i_str, str, method, best_ARI, best_VI) %>%
                    rename(ARI = best_ARI, VI = best_VI) %>%
                    pivot_longer(cols = c(ARI, VI),
                                 names_to = "quantity",
                                 values_to = "value"))


theme_set(theme_bw())

P1 <- combined %>%
  mutate(method = factor(method, levels = c("Pose", "BayesSpace", "UTAG"))) %>%
  filter(quantity == "ARI") %>%
  ggplot(aes(x = as.factor(str), y = value)) +
  geom_boxplot() +
  xlab("Degree of spatial dependence")+
  ylab("ARI")+
  theme(text = element_text(size=15))+
  facet_wrap(~method)
P1+ggview::canvas(h=4,w=12)

ggsave(plot = P1,
       filename = "Simulation_VaryingSpatial_competitors/results/ARI_plot.png",
       h=4,w=12)



P2 <- combined %>%
  mutate(method = factor(method, levels = c("Pose", "BayesSpace", "UTAG"))) %>%
  filter(quantity == "ARI") %>%
  mutate(str = paste("Spatial dependence level:",str),
         str = factor(str, levels =
                        c("Spatial dependence level: 0.5",
                          "Spatial dependence level: 1",
                          "Spatial dependence level: 2",
                          "Spatial dependence level: 10"))) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  xlab("Degree of spatial dependence")+
  ylab("ARI")+
  theme(text = element_text(size=15))+
  facet_wrap(~as.factor(str),nrow=1)
P2+ggview::canvas(h=4,w=15)
ggsave(plot = P2,
       filename = "Simulation_VaryingSpatial_competitors/results/ARI_plot_v2.png",
       h=4,w=15)


P3 <- combined %>%
  filter(method == "Pose") %>%
  filter(quantity %in% c("ARI_bicl","RMSE_bicl")) %>%
  mutate(quantity = case_when(quantity == "ARI_bicl" ~ "Pose: biclustering ARI",
                              quantity == "RMSE_bicl" ~ "Pose: biclustering RMSE",
  )) %>%
  ggplot(aes(x = as.factor(str), y = value)) +
  geom_boxplot() +
  xlab("Degree of spatial dependence")+
  ylab("Value")+
  theme(text = element_text(size=15))+
  facet_wrap(~as.factor(quantity), scales = "free_y")
P3+ggview::canvas(h=4,w=12)
ggsave(plot = P3,
       filename = "Simulation_VaryingSpatial_competitors/results/Bicl_plot.png",
       h=4,w=12)
