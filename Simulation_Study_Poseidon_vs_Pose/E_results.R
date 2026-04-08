library(Poser)
library(purrr)
library(tidyverse)
library(stringr)
library(latex2exp)
library(ggview)

sers <-  readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/posers_dec25.RDS")
idon <-  readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/poseid_dec25.RDS")
datas <- readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/data_and_true_clustering.RDS")

# -------------------------------------------------------------------------
length(sers)

D0_cols = lapply(sers[[1]], function(x) estimate_clust_pose(res =  x)$col_cl)
D1_cols = lapply(sers[[2]], function(x) estimate_clust_pose(res =  x)$col_cl)
D2_cols = lapply(sers[[3]], function(x) estimate_clust_pose(res =  x)$col_cl)
D3_cols = lapply(sers[[4]], function(x) estimate_clust_pose(res =  x)$col_cl)
D4_cols = lapply(sers[[5]], function(x) estimate_clust_pose(res =  x)$col_cl)

D0_bicl = lapply(sers[[1]], function(x) estimate_clust_pose(res =  x)$bicl)
D1_bicl = lapply(sers[[2]], function(x) estimate_clust_pose(res =  x)$bicl)
D2_bicl = lapply(sers[[3]], function(x) estimate_clust_pose(res =  x)$bicl)
D3_bicl = lapply(sers[[4]], function(x) estimate_clust_pose(res =  x)$bicl)
D4_bicl = lapply(sers[[5]], function(x) estimate_clust_pose(res =  x)$bicl)

mse_0 = map2_vec(datas$D0list,D0_bicl, ~mean((.x-.y)^2))
mse_1 = map2_vec(datas$D1list,D1_bicl, ~mean((.x-.y)^2))
mse_2 = map2_vec(datas$D2list,D2_bicl, ~mean((.x-.y)^2))
mse_3 = map2_vec(datas$D3list,D3_bicl, ~mean((.x-.y)^2))
mse_4 = map2_vec(datas$D4list,D4_bicl, ~mean((.x-.y)^2))

boxplot(cbind(mse_1,mse_2,mse_3,mse_4,mse_0))

ari_0 = map2_vec(datas$CClist,D0_cols, ~mcclust::arandi(.x,.y))
ari_1 = map2_vec(datas$CClist,D1_cols, ~mcclust::arandi(.x,.y))
ari_2 = map2_vec(datas$CClist,D2_cols, ~mcclust::arandi(.x,.y))
ari_3 = map2_vec(datas$CClist,D3_cols, ~mcclust::arandi(.x,.y))
ari_4 = map2_vec(datas$CClist,D4_cols, ~mcclust::arandi(.x,.y))
boxplot(cbind(ari_1,ari_2,ari_3,ari_4,ari_0))


D0_bicl_IND = lapply(sers[[1]], function(x) as.numeric(as.factor(estimate_clust_pose(res =  x)$bicl)))
D1_bicl_IND = lapply(sers[[2]], function(x) as.numeric(as.factor(estimate_clust_pose(res =  x)$bicl)))
D2_bicl_IND = lapply(sers[[3]], function(x) as.numeric(as.factor(estimate_clust_pose(res =  x)$bicl)))
D3_bicl_IND = lapply(sers[[4]], function(x) as.numeric(as.factor(estimate_clust_pose(res =  x)$bicl)))
D4_bicl_IND = lapply(sers[[5]], function(x) as.numeric(as.factor(estimate_clust_pose(res =  x)$bicl)))


ari0_observ_cl = map2_vec(D0_bicl_IND,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari1_observ_cl = map2_vec(D1_bicl_IND,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari2_observ_cl = map2_vec(D2_bicl_IND,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari3_observ_cl = map2_vec(D3_bicl_IND,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari4_observ_cl = map2_vec(D4_bicl_IND,datas$RClist, ~mcclust::arandi(.x,c(.y)))

ncl0_observ_cl = map_vec(D0_bicl_IND, ~length(unique(c(.x))))
ncl1_observ_cl = map_vec(D1_bicl_IND, ~length(unique(c(.x))))
ncl2_observ_cl = map_vec(D2_bicl_IND, ~length(unique(c(.x))))
ncl3_observ_cl = map_vec(D3_bicl_IND, ~length(unique(c(.x))))
ncl4_observ_cl = map_vec(D4_bicl_IND, ~length(unique(c(.x))))

boxplot(cbind(ari1_observ_cl,
              ari2_observ_cl,
              ari3_observ_cl,
              ari4_observ_cl,
              ari0_observ_cl))

boxplot(cbind(ncl1_observ_cl,
              ncl2_observ_cl,
              ncl3_observ_cl,
              ncl4_observ_cl,
              ncl0_observ_cl))

# Poseidon ----------------------------------------------------------------

D1_idon_cols = lapply(idon[[1]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)
D2_idon_cols = lapply(idon[[2]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)
D3_idon_cols = lapply(idon[[3]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)
D4_idon_cols = lapply(idon[[4]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)

D1_idon_bicl = lapply(idon[[1]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)
D2_idon_bicl = lapply(idon[[2]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)
D3_idon_bicl = lapply(idon[[3]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)
D4_idon_bicl = lapply(idon[[4]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)


image(D1_idon_bicl[[1]][[1]])
image(D1_idon_bicl[[1]][[2]])
image(datas$D4list[[1]])

mse_1_a = map2_vec(datas$D0list,D1_idon_bicl, ~mean((.x-.y[[1]])^2))
mse_1_b = map2_vec(datas$D1list,D1_idon_bicl, ~mean((.x-.y[[2]])^2))
mse_2_a = map2_vec(datas$D0list,D2_idon_bicl, ~mean((.x-.y[[1]])^2))
mse_2_b = map2_vec(datas$D2list,D2_idon_bicl, ~mean((.x-.y[[2]])^2))
mse_3_a = map2_vec(datas$D0list,D3_idon_bicl, ~mean((.x-.y[[1]])^2))
mse_3_b = map2_vec(datas$D3list,D3_idon_bicl, ~mean((.x-.y[[2]])^2))
mse_4_a = map2_vec(datas$D0list,D4_idon_bicl, ~mean((.x-.y[[1]])^2))
mse_4_b = map2_vec(datas$D4list,D4_idon_bicl, ~mean((.x-.y[[2]])^2))

boxplot(cbind(mse_1_a,
              mse_2_a,
              mse_3_a,
              mse_4_a,
              mse_1_b,
              mse_2_b,
              mse_3_b,
              mse_4_b))


ari_idon_1 = map2_vec(datas$CClist,D1_idon_cols, ~mcclust::arandi(.x,.y))
ari_idon_2 = map2_vec(datas$CClist,D2_idon_cols, ~mcclust::arandi(.x,.y))
ari_idon_3 = map2_vec(datas$CClist,D3_idon_cols, ~mcclust::arandi(.x,.y))
ari_idon_4 = map2_vec(datas$CClist,D4_idon_cols, ~mcclust::arandi(.x,.y))

boxplot(cbind(ari_idon_1,
              ari_idon_2,
              ari_idon_3,
              ari_idon_4))


D1_idon_bicl_IND_a = lapply(idon[[1]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[1]])))
D2_idon_bicl_IND_a = lapply(idon[[2]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[1]])))
D3_idon_bicl_IND_a = lapply(idon[[3]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[1]])))
D4_idon_bicl_IND_a = lapply(idon[[4]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[1]])))
D1_idon_bicl_IND_b = lapply(idon[[1]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[2]])))
D2_idon_bicl_IND_b = lapply(idon[[2]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[2]])))
D3_idon_bicl_IND_b = lapply(idon[[3]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[2]])))
D4_idon_bicl_IND_b = lapply(idon[[4]], function(x) as.numeric(as.factor(POSEIDON::estimate_clust_poseidon(q = x)$bicl[[2]])))

plot(D4_idon_bicl_IND_b[[1]])
plot(c(datas$RClist[[1]]))
image(matrix(D4_idon_bicl_IND_b[[1]],200,50))
image(matrix(datas$RClist[[1]],200,50))

mcclust::arandi(D4_idon_bicl_IND_b[[1]],datas$RClist[[1]])

image(D4_idon_bicl[[1]][[1]])
image(D4_idon_bicl[[1]][[2]])


ari1_idon_observ_cl_a = map2_vec(D1_idon_bicl_IND_a,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari2_idon_observ_cl_a = map2_vec(D2_idon_bicl_IND_a,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari3_idon_observ_cl_a = map2_vec(D3_idon_bicl_IND_a,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari4_idon_observ_cl_a = map2_vec(D4_idon_bicl_IND_a,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari1_idon_observ_cl_b = map2_vec(D1_idon_bicl_IND_b,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari2_idon_observ_cl_b = map2_vec(D2_idon_bicl_IND_b,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari3_idon_observ_cl_b = map2_vec(D3_idon_bicl_IND_b,datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari4_idon_observ_cl_b = map2_vec(D4_idon_bicl_IND_b,datas$RClist, ~mcclust::arandi(.x,c(.y)))


ncl1_idon_observ_cl_a = map_vec(D1_idon_bicl_IND_a, ~length(unique(c(.x))))
ncl2_idon_observ_cl_a = map_vec(D2_idon_bicl_IND_a, ~length(unique(c(.x))))
ncl3_idon_observ_cl_a = map_vec(D3_idon_bicl_IND_a, ~length(unique(c(.x))))
ncl4_idon_observ_cl_a = map_vec(D4_idon_bicl_IND_a, ~length(unique(c(.x))))
ncl1_idon_observ_cl_b = map_vec(D1_idon_bicl_IND_b, ~length(unique(c(.x))))
ncl2_idon_observ_cl_b = map_vec(D2_idon_bicl_IND_b, ~length(unique(c(.x))))
ncl3_idon_observ_cl_b = map_vec(D3_idon_bicl_IND_b, ~length(unique(c(.x))))
ncl4_idon_observ_cl_b = map_vec(D4_idon_bicl_IND_b, ~length(unique(c(.x))))



boxplot(cbind(ari1_idon_observ_cl_a,
              ari2_idon_observ_cl_a,
              ari3_idon_observ_cl_a,
              ari4_idon_observ_cl_a,
              ari1_idon_observ_cl_b,
              ari2_idon_observ_cl_b,
              ari3_idon_observ_cl_b,
              ari4_idon_observ_cl_b))
boxplot(cbind(ncl1_idon_observ_cl_a,
              ncl2_idon_observ_cl_a,
              ncl3_idon_observ_cl_a,
              ncl4_idon_observ_cl_a,
              ncl1_idon_observ_cl_b,
              ncl2_idon_observ_cl_b,
              ncl3_idon_observ_cl_b,
              ncl4_idon_observ_cl_b))


boxplot(cbind(ari1_observ_cl,
              ari2_observ_cl,
              ari3_observ_cl,
              ari4_observ_cl))
boxplot(cbind(ari1_idon_observ_cl_b,
              ari2_idon_observ_cl_b,
              ari3_idon_observ_cl_b,
              ari4_idon_observ_cl_b),col=2,add=T)



data.frame(ari = ari_0, data = "0" )
data.frame(ari = ari_1, data = "1" )
data.frame(ari = ari_2, data = "2" )
data.frame(ari = ari_3, data = "3" )
data.frame(ari = ari_4, data = "4" )


# ari_col and row_clu -------------------------------------------------------------
X_POSE <- data.frame(rbind(c(
  #c(mean(ari_0), sd(ari_0),
    mean(ari_4), sd(ari_4),
    mean(ari_3), sd(ari_3),
    mean(ari_2), sd(ari_2),
    mean(ari_1), sd(ari_1)),
  #mean(ari0_observ_cl), sd(ari0_observ_cl),
  c( mean(ari4_observ_cl), sd(ari4_observ_cl),
    mean(ari3_observ_cl), sd(ari3_observ_cl),
    mean(ari2_observ_cl), sd(ari2_observ_cl),
    mean(ari1_observ_cl), sd(ari1_observ_cl))))

# recall that you switched the order 4-3-2-1 in the dataset
X_POSE <- X_POSE %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))

knitr::kable(X_POSE,format = "latex")


####################################################################
# POSE ON REFERENCE DATASET
mean(ari_0); sd(ari_0)
mean(ari0_observ_cl); sd(ari0_observ_cl)
####################################################################
X_POSE
####################################################################


X_IDON <- data.frame(
  rbind(c(mean(ari_idon_4), sd(ari_idon_4),
          mean(ari_idon_3), sd(ari_idon_3),
          mean(ari_idon_2), sd(ari_idon_2),
          mean(ari_idon_1), sd(ari_idon_1)),
  c( mean(ari4_idon_observ_cl_a), sd(ari4_idon_observ_cl_a),
     mean(ari3_idon_observ_cl_a), sd(ari3_idon_observ_cl_a),
     mean(ari2_idon_observ_cl_a), sd(ari2_idon_observ_cl_a),
     mean(ari1_idon_observ_cl_a), sd(ari1_idon_observ_cl_a)),
  c( mean(ari4_idon_observ_cl_b), sd(ari4_idon_observ_cl_b),
     mean(ari3_idon_observ_cl_b), sd(ari3_idon_observ_cl_b),
     mean(ari2_idon_observ_cl_b), sd(ari2_idon_observ_cl_b),
     mean(ari1_idon_observ_cl_b), sd(ari1_idon_observ_cl_b))))


# recall that you switched the order 4-3-2-1 in the dataset
X_IDON <- X_IDON %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))

knitr::kable(X_IDON,format = "latex",digits = 3)


rbind(X_POSE,X_IDON) %>%  knitr::kable(format = "latex",digits = 3)



############### PLOTTING ROW CLUSTERING
# obs_col_clu -------------------------------------------------------------
YY <- data.frame(ari1_observ_cl,
                 ari2_observ_cl,
                 ari3_observ_cl,
                 ari4_observ_cl,
                 ari0_observ_cl)

YY <- (data.frame(
  ari1_idon_observ_cl_a,
  ari2_idon_observ_cl_a,
  ari3_idon_observ_cl_a,
  ari4_idon_observ_cl_a,
  ari1_idon_observ_cl_b,
  ari2_idon_observ_cl_b,
  ari3_idon_observ_cl_b,
  ari4_idon_observ_cl_b))


YY <- reshape::melt(YY) %>% mutate(var = str_sub(variable,-1,-1))
YY <- YY %>% mutate(D = case_when(substr(variable,1,4)=="ari1"~("4"),
                                  substr(variable,1,4)=="ari2"~("3"),
                                  substr(variable,1,4)=="ari3"~("2"),
                                  substr(variable,1,4)=="ari4"~("1")),
                    single_ari = case_when(substr(variable,1,4)=="ari1"~mean(ari1_observ_cl),
                                           substr(variable,1,4)=="ari2"~mean(ari2_observ_cl),
                                           substr(variable,1,4)=="ari3"~mean(ari3_observ_cl),
                                           substr(variable,1,4)=="ari4"~mean(ari4_observ_cl)),
                    lab_x = case_when(var == "a" ~ "Reference",
                                      var == "b" ~ "Noisy"))

appender0 <- function(labs) {
  lapply(labs, function(s) {
    bquote( "{" ~ bold(Y)^(0) * "," ~ bold(Y)^(.(s)) ~ "}" )
  })
}
YY$lab_x= as.factor(YY$lab_x)
YY$lab_x = factor(YY$lab_x, levels = levels(YY$lab_x)[2:1])

ggplot(YY)+theme_bw()+
  geom_hline(aes(yintercept = single_ari),col="darkgray",lty=2)+
  geom_boxplot(aes(y=value,x=lab_x,group=variable))+
  theme(text = element_text(size=18))+
  facet_wrap(~D,nrow = 1,labeller =
               as_labeller(appender0,
                           default = label_parsed))+
  xlab("")+ylab("ARI - Row Clustering")+
  canvas(w=16,h=4)


ggsave("Simulation_Study_Poseidon_vs_Pose/Output/PLOT/SIM_ari_RC_dec25.png",width = 16,height = 4)
ggsave("Simulation_Study_Poseidon_vs_Pose/Output/PLOT/SIM_ari_RC_dec25.pdf",width = 16,height = 4)
ggsave("Simulation_Study_Poseidon_vs_Pose/Output/PLOT/SIM_ari_RC_dec25.eps",width = 16,height = 4)



D = rbind( cbind(reshape2::melt(datas$D0list[[1]]),d= 0),
           cbind(reshape2::melt(datas$D1list[[1]]),d= 1),
           cbind(reshape2::melt(datas$D2list[[1]]),d= 2),
           cbind(reshape2::melt(datas$D3list[[1]]),d= 3),
           cbind(reshape2::melt(datas$D4list[[1]]),d= 4))
plot(D$value)

# appender <- function(string)
#   TeX(paste("${D}^{(",string,")}$"))

appender <- function(labs) {
  lapply(labs, function(s) {
    bquote( ~ bold(Y)^(.(s)) )
  })
}


appender(1)

ggplot(D)+theme_bw()+
  geom_tile(aes(Var2,rev(Var1),fill=value))+
  facet_wrap(~d,nrow = 1,labeller =
               as_labeller(appender,
                           default = label_parsed))+
  scale_fill_gradient2("", low = "red",high = 4,midpoint = .0)+
  geom_segment(aes(x=10.5,xend=10.5,y=0.5,yend=50.5))+
  theme(legend.position = "bottom",text = element_text(size=18))+
  ggtitle("Single realization of the data used in the simulation study")+
  geom_segment(aes(y = 25.5,yend=25.5,x=0,xend=10.5))+
  xlab("")+ylab("")+
  canvas(w=20,h=7)


ggsave("Simulation_Study_Poseidon_vs_Pose/Output/PLOT/dataset_dec25.pdf",width = 20,height = 7)
