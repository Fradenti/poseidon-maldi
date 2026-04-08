library(Poser)
library(purrr)
library(tidyverse)
library(stringr)
library(latex2exp)
library(ggview)

sers  <- readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/posers_dec25.RDS")
stack <- readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/posers_stacked_dec25.RDS")
idon  <- readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/poseid_dec25.RDS")
datas <- readRDS("Simulation_Study_Poseidon_vs_Pose/Output/RDS/data_and_true_clustering.RDS")


a = estimate_clust_pose(stack[[1]][[1]])
image(a$bicl)

# Posers -------------------------------------------------------------------------
# Estimate column clusters
D0_cols <- lapply(sers[[1]], function(x) estimate_clust_pose(res = x)$col_cl)
D1_cols <- lapply(sers[[2]], function(x) estimate_clust_pose(res = x)$col_cl)
D2_cols <- lapply(sers[[3]], function(x) estimate_clust_pose(res = x)$col_cl)
D3_cols <- lapply(sers[[4]], function(x) estimate_clust_pose(res = x)$col_cl)
D4_cols <- lapply(sers[[5]], function(x) estimate_clust_pose(res = x)$col_cl)

# Estimate biclustering
D0_bicl <- lapply(sers[[1]], function(x) estimate_clust_pose(res = x)$bicl)
D1_bicl <- lapply(sers[[2]], function(x) estimate_clust_pose(res = x)$bicl)
D2_bicl <- lapply(sers[[3]], function(x) estimate_clust_pose(res = x)$bicl)
D3_bicl <- lapply(sers[[4]], function(x) estimate_clust_pose(res = x)$bicl)
D4_bicl <- lapply(sers[[5]], function(x) estimate_clust_pose(res = x)$bicl)

mse_0 <- map2_vec(datas$D0list, D0_bicl, ~mean((.x-.y)^2))
mse_1 <- map2_vec(datas$D1list, D1_bicl, ~mean((.x-.y)^2))
mse_2 <- map2_vec(datas$D2list, D2_bicl, ~mean((.x-.y)^2))
mse_3 <- map2_vec(datas$D3list, D3_bicl, ~mean((.x-.y)^2))
mse_4 <- map2_vec(datas$D4list, D4_bicl, ~mean((.x-.y)^2))

ari_0 <- map2_vec(datas$CClist, D0_cols, ~mcclust::arandi(.x,.y))
ari_1 <- map2_vec(datas$CClist, D1_cols, ~mcclust::arandi(.x,.y))
ari_2 <- map2_vec(datas$CClist, D2_cols, ~mcclust::arandi(.x,.y))
ari_3 <- map2_vec(datas$CClist, D3_cols, ~mcclust::arandi(.x,.y))
ari_4 <- map2_vec(datas$CClist, D4_cols, ~mcclust::arandi(.x,.y))
boxplot(cbind(ari_1,ari_2,ari_3,ari_4,ari_0))


D0_bicl_IND = lapply(sers[[1]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))
D1_bicl_IND = lapply(sers[[2]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))
D2_bicl_IND = lapply(sers[[3]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))
D3_bicl_IND = lapply(sers[[4]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))
D4_bicl_IND = lapply(sers[[5]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))


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


# Stacked -------------------------------------------------------------------------

D01s_cols = lapply(stack[[1]], function(x) estimate_clust_pose(res = x)$col_cl)
D02s_cols = lapply(stack[[2]], function(x) estimate_clust_pose(res = x)$col_cl)
D03s_cols = lapply(stack[[3]], function(x) estimate_clust_pose(res = x)$col_cl)
D04s_cols = lapply(stack[[4]], function(x) estimate_clust_pose(res = x)$col_cl)

D01s_bicl = lapply(stack[[1]], function(x) estimate_clust_pose(res = x)$bicl)
D02s_bicl = lapply(stack[[2]], function(x) estimate_clust_pose(res = x)$bicl)
D03s_bicl = lapply(stack[[3]], function(x) estimate_clust_pose(res = x)$bicl)
D04s_bicl = lapply(stack[[4]], function(x) estimate_clust_pose(res = x)$bicl)

D01s = map2(datas$D0list,datas$D1list, ~ rbind(.x,.y))
D02s = map2(datas$D0list,datas$D2list, ~ rbind(.x,.y))
D03s = map2(datas$D0list,datas$D3list, ~ rbind(.x,.y))
D04s = map2(datas$D0list,datas$D4list, ~ rbind(.x,.y))

mses_01 = map2_vec(D01s,D01s_bicl, ~mean((.x-.y)^2))
mses_02 = map2_vec(D02s,D02s_bicl, ~mean((.x-.y)^2))
mses_03 = map2_vec(D03s,D03s_bicl, ~mean((.x-.y)^2))
mses_04 = map2_vec(D04s,D04s_bicl, ~mean((.x-.y)^2))
boxplot(cbind(mses_01,mses_02,mses_03,mses_04))

ari_01s = map2_vec(datas$CClist,D01s_cols, ~mcclust::arandi(.x,.y))
ari_02s = map2_vec(datas$CClist,D02s_cols, ~mcclust::arandi(.x,.y))
ari_03s = map2_vec(datas$CClist,D03s_cols, ~mcclust::arandi(.x,.y))
ari_04s = map2_vec(datas$CClist,D04s_cols, ~mcclust::arandi(.x,.y))
boxplot(cbind(ari_01s,ari_02s,ari_03s,ari_04s))


D01s_bicl_IND = lapply(stack[[1]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))
D02s_bicl_IND = lapply(stack[[2]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))
D03s_bicl_IND = lapply(stack[[3]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))
D04s_bicl_IND = lapply(stack[[4]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl)))

# -------------------------------------------------------------------------
RCstack = map2(datas$RClist, datas$RClist, ~rbind(.x,.y))

ari01s_observ_cl = map2_vec(D01s_bicl_IND, RCstack, ~mcclust::arandi(.x,c(.y)))
ari02s_observ_cl = map2_vec(D02s_bicl_IND, RCstack, ~mcclust::arandi(.x,c(.y)))
ari03s_observ_cl = map2_vec(D03s_bicl_IND, RCstack, ~mcclust::arandi(.x,c(.y)))
ari04s_observ_cl = map2_vec(D04s_bicl_IND, RCstack, ~mcclust::arandi(.x,c(.y)))
ncl01s_observ_cl = map_vec(D01s_bicl_IND, ~length(unique(c(.x))))
ncl02s_observ_cl = map_vec(D02s_bicl_IND, ~length(unique(c(.x))))
ncl03s_observ_cl = map_vec(D03s_bicl_IND, ~length(unique(c(.x))))
ncl04s_observ_cl = map_vec(D04s_bicl_IND, ~length(unique(c(.x))))

boxplot(cbind(ari01s_observ_cl,
              ari02s_observ_cl,
              ari03s_observ_cl,
              ari04s_observ_cl),ylim=c(0,1))

boxplot(cbind(ncl01s_observ_cl,
              ncl02s_observ_cl,
              ncl03s_observ_cl,
              ncl04s_observ_cl))

## Separate results ------------------------------------------------------------
# obtain the matrices
D01s_bicl_IND_ref <- lapply(stack[[1]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[1:50,])))
D02s_bicl_IND_ref <- lapply(stack[[2]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[1:50,])))
D03s_bicl_IND_ref <- lapply(stack[[3]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[1:50,])))
D04s_bicl_IND_ref <- lapply(stack[[4]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[1:50,])))

D01s_bicl_IND_noi <- lapply(stack[[1]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[51:100,])))
D02s_bicl_IND_noi <- lapply(stack[[2]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[51:100,])))
D03s_bicl_IND_noi <- lapply(stack[[3]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[51:100,])))
D04s_bicl_IND_noi <- lapply(stack[[4]], function(x) as.numeric(as.factor(estimate_clust_pose(res = x)$bicl[51:100,])))

ari01s_observ_cl_ref <- map2_vec(D01s_bicl_IND_ref, datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari02s_observ_cl_ref <- map2_vec(D02s_bicl_IND_ref, datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari03s_observ_cl_ref <- map2_vec(D03s_bicl_IND_ref, datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari04s_observ_cl_ref <- map2_vec(D04s_bicl_IND_ref, datas$RClist, ~mcclust::arandi(.x,c(.y)))

ari01s_observ_cl_noi <- map2_vec(D01s_bicl_IND_noi, datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari02s_observ_cl_noi <- map2_vec(D02s_bicl_IND_noi, datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari03s_observ_cl_noi <- map2_vec(D03s_bicl_IND_noi, datas$RClist, ~mcclust::arandi(.x,c(.y)))
ari04s_observ_cl_noi <- map2_vec(D04s_bicl_IND_noi, datas$RClist, ~mcclust::arandi(.x,c(.y)))


boxplot(cbind(ari01s_observ_cl_ref,
              ari02s_observ_cl_ref,
              ari03s_observ_cl_ref,
              ari04s_observ_cl_ref))
boxplot(cbind(ari01s_observ_cl_noi,
              ari02s_observ_cl_noi,
              ari03s_observ_cl_noi,
              ari04s_observ_cl_noi))

# Poseidon ----------------------------------------------------------------

D1_idon_cols = lapply(idon[[1]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)
D2_idon_cols = lapply(idon[[2]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)
D3_idon_cols = lapply(idon[[3]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)
D4_idon_cols = lapply(idon[[4]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$col_cl)

D1_idon_bicl = lapply(idon[[1]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)
D2_idon_bicl = lapply(idon[[2]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)
D3_idon_bicl = lapply(idon[[3]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)
D4_idon_bicl = lapply(idon[[4]], function(x) POSEIDON::estimate_clust_poseidon(q = x)$bicl)


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


ari_idon_1 = map2_vec(datas$CClist, D1_idon_cols, ~mcclust::arandi(.x,.y))
ari_idon_2 = map2_vec(datas$CClist, D2_idon_cols, ~mcclust::arandi(.x,.y))
ari_idon_3 = map2_vec(datas$CClist, D3_idon_cols, ~mcclust::arandi(.x,.y))
ari_idon_4 = map2_vec(datas$CClist, D4_idon_cols, ~mcclust::arandi(.x,.y))

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


ari1_idon_observ_cl_a = map2_vec(D1_idon_bicl_IND_a, datas$RClist, ~mcclust::arandi(.x,c(.y)))
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
boxplot(cbind(ari01s_observ_cl,
              ari02s_observ_cl,
              ari03s_observ_cl,
              ari04s_observ_cl),col=3,add=T)
boxplot(cbind(ari01s_observ_cl_noi,
              ari02s_observ_cl_noi,
              ari03s_observ_cl_noi,
              ari04s_observ_cl_noi),col=4,add=T)



data.frame(ari = ari_0, data = "0" )
data.frame(ari = ari_1, data = "1" )
data.frame(ari = ari_2, data = "2" )
data.frame(ari = ari_3, data = "3" )
data.frame(ari = ari_4, data = "4" )


# COLLECTING RESULTS IN TABLE --------------------------------------------------

X_POSE <- data.frame(rbind(
  #mean(ari_0), sd(ari_0),
  c(mean(ari_4), sd(ari_4),
    mean(ari_3), sd(ari_3),
    mean(ari_2), sd(ari_2),
    mean(ari_1), sd(ari_1)),
  #mean(ari0_observ_cl), sd(ari0_observ_cl),
    c( mean(ari4_observ_cl), sd(ari4_observ_cl),
    mean(ari3_observ_cl), sd(ari3_observ_cl),
    mean(ari2_observ_cl), sd(ari2_observ_cl),
    mean(ari1_observ_cl), sd(ari1_observ_cl))))
rownames(X_POSE) <- c("ARI CC", "ARI RC")
# recall that you switched the order 4-3-2-1 in the dataset
X_POSE <- X_POSE %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))
knitr::kable(X_POSE,format = "latex")

X_IDON <- data.frame(
  rbind(c(mean(ari_idon_4), sd(ari_idon_4),
          mean(ari_idon_3), sd(ari_idon_3),
          mean(ari_idon_2), sd(ari_idon_2),
          mean(ari_idon_1), sd(ari_idon_1)),
        c(mean(ari4_idon_observ_cl_a), sd(ari4_idon_observ_cl_a),
          mean(ari3_idon_observ_cl_a), sd(ari3_idon_observ_cl_a),
          mean(ari2_idon_observ_cl_a), sd(ari2_idon_observ_cl_a),
          mean(ari1_idon_observ_cl_a), sd(ari1_idon_observ_cl_a)),
        c(mean(ari4_idon_observ_cl_b), sd(ari4_idon_observ_cl_b),
          mean(ari3_idon_observ_cl_b), sd(ari3_idon_observ_cl_b),
          mean(ari2_idon_observ_cl_b), sd(ari2_idon_observ_cl_b),
          mean(ari1_idon_observ_cl_b), sd(ari1_idon_observ_cl_b))))



# recall that you switched the order 4-3-2-1 in the dataset
X_IDON <- X_IDON %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))
rownames(X_IDON) <- c("ARI CC", "ARI RC - REF", "ARI RC - OTHER")

knitr::kable(X_IDON,format = "latex",digits = 3)




X_STACK <- data.frame(rbind(
  c(mean(ari_04s), sd(ari_04s),
    mean(ari_03s), sd(ari_03s),
    mean(ari_02s), sd(ari_02s),
    mean(ari_01s), sd(ari_01s)),
  c(mean(ari04s_observ_cl), sd(ari04s_observ_cl),
    mean(ari03s_observ_cl), sd(ari03s_observ_cl),
    mean(ari02s_observ_cl), sd(ari02s_observ_cl),
    mean(ari01s_observ_cl), sd(ari01s_observ_cl))))

# recall that you switched the order 4-3-2-1 in the dataset
X_STACK <- X_STACK %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))

rownames(X_STACK) <- c("ARI CC", "ARI RC")
knitr::kable(X_STACK,format = "latex",digits = 3)


## REF VS NOISY di dataset stacked ma valutati separati

X_STACKED_SEP <- data.frame(rbind(
  c(mean(ari04s_observ_cl_ref), sd(ari04s_observ_cl_ref),
    mean(ari03s_observ_cl_ref), sd(ari03s_observ_cl_ref),
    mean(ari02s_observ_cl_ref), sd(ari02s_observ_cl_ref),
    mean(ari01s_observ_cl_ref), sd(ari01s_observ_cl_ref)),
  c(mean(ari04s_observ_cl_noi), sd(ari04s_observ_cl_noi),
    mean(ari03s_observ_cl_noi), sd(ari03s_observ_cl_noi),
    mean(ari02s_observ_cl_noi), sd(ari02s_observ_cl_noi),
    mean(ari01s_observ_cl_noi), sd(ari01s_observ_cl_noi))))

# recall that you switched the order 4-3-2-1 in the dataset
X_STACKED_SEP <- X_STACKED_SEP %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))

rownames(X_STACKED_SEP) <- c("ARI RC - REF", "ARI RC - OTHER")
knitr::kable(X_STACKED_SEP,format = "latex",digits = 3)


rbind(
X_POSE,
X_IDON,
X_STACK,
X_STACKED_SEP) %>% knitr::kable(format = "latex",digits = 3)


data_list <- list(
  "POSE"        = X_POSE,
  "IDON"        = X_IDON,
  "STACK"       = X_STACK,
  "STACKED_SEP" = X_STACKED_SEP
)

bind_rows(data_list, .id = "Method")  %>%
  knitr::kable(format = "latex",
        digits = 3,
        booktabs = TRUE)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Number of clusters -----------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

X_colclu <- data.frame(rbind(
  c(mean(ncl4_observ_cl), sd(ncl4_observ_cl),
    mean(ncl3_observ_cl), sd(ncl3_observ_cl),
    mean(ncl2_observ_cl), sd(ncl2_observ_cl),
    mean(ncl1_observ_cl), sd(ncl1_observ_cl),
    mean(ncl0_observ_cl), sd(ncl0_observ_cl))
))


# recall that you switched the order 4-3-2-1 in the dataset
X_colclu <- X_colclu %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"),
                                               X10 = paste0("(",X10,")"))

knitr::kable(X_colclu,format = "latex")


X_colclu <- data.frame(
  rbind(c(mean(ncl4_idon_observ_cl_a), sd(ncl4_idon_observ_cl_a),
          mean(ncl3_idon_observ_cl_a), sd(ncl3_idon_observ_cl_a),
          mean(ncl2_idon_observ_cl_a), sd(ncl2_idon_observ_cl_a),
          mean(ncl1_idon_observ_cl_a), sd(ncl1_idon_observ_cl_a)),
        c(mean(ncl4_idon_observ_cl_b), sd(ncl4_idon_observ_cl_b),
          mean(ncl3_idon_observ_cl_b), sd(ncl3_idon_observ_cl_b),
          mean(ncl2_idon_observ_cl_b), sd(ncl2_idon_observ_cl_b),
          mean(ncl1_idon_observ_cl_b), sd(ncl1_idon_observ_cl_b))))


# recall that you switched the order 4-3-2-1 in the dataset
X_colclu <- X_colclu %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))

knitr::kable(X_colclu,format = "latex",digits = 3)

X_colclu <- data.frame(
  rbind(c(mean(ncl04s_observ_cl), sd(ncl04s_observ_cl),
          mean(ncl03s_observ_cl), sd(ncl03s_observ_cl),
          mean(ncl02s_observ_cl), sd(ncl02s_observ_cl),
          mean(ncl01s_observ_cl), sd(ncl01s_observ_cl))))


# recall that you switched the order 4-3-2-1 in the dataset
X_colclu <- X_colclu %>% round(.,3) %>% mutate(X2 = paste0("(",X2,")"),
                                               X4 = paste0("(",X4,")"),
                                               X6 = paste0("(",X6,")"),
                                               X8 = paste0("(",X8,")"))

knitr::kable(X_colclu,format = "latex",digits = 3)
