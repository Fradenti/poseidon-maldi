library(tidyverse)
library(raster)
library(Poser)
library(patchwork)

# Lipids ------------------------------------------------------------------

K <- 30
L <- 40
pixels     <- as.matrix(readRDS("CRCC_Analysis/Data/FinalMatrix_Coor.RDS"))
val        <- as.matrix(readRDS("CRCC_Analysis/Data/FinalMatrix_Lipids.RDS"))
rownames(val) <- colnames(val) <- NULL
nnList_cpp <- readRDS("CRCC_Analysis/Data/list_ind_resolution_4123_cpp.RDS")

# check NN indexes are ok ---------------------------------------------------------

plot(pixels)
points(pixels[1,2]~pixels[1,1],col=2,pch=21,bg=4)
points(pixels[nnList_cpp[[1]]+1,],col=2,pch=21,bg="1")
points(pixels[10,2]~pixels[10,1],col=2,pch=21,bg=2)
points(pixels[nnList_cpp[[10]]+1,],col=2,pch=21,bg="1")


# hyperpar ---------------------------------------------------------

hyperpar = list(
  a  = 1e-4,
  b  = 1e-4,
  m0 = 0,
  k0 = .01,
  c0 = 3,
  d0 = 2,
  s1 = 1,
  s2 = 1,
  s3 = 1,
  s4 = 1)


# transform data -----------------------------------------------------------
mapper2 <- function(Y) {
  m = min(Y)-1e-3
  M = max(Y)+1e-3
  z = (Y-m)/(M-m)
  qnorm(z) - mean(qnorm(z))
}

r = mapper2(t(val))


# Multistart --------------------------------------------------------------
nsim <- 250
for(i in 1:nsim){
  set.seed(1234554*i)
  par <- maotai::kmeanspp(t(r),k=K)
  RHOi <- from_partition_to_probs(par,
                                  K = K,
                                  perturb = .5)
  N <- nrow(r)
  XI_ikl <- array(NA,c(N,K,L))
  for(k in 1:K){
    par_row <- maotai::kmeanspp(as.matrix((r[,par==k])),k=5)
    XI_ikl[,k,] = from_partition_to_probs(par_row, K=L, perturb = .2)
  }

  MOD = Poser::POSEr_fiSAN(Y = r,
                           NN_inds = nnList_cpp,
                           RHO_init = RHOi,
                           XI_init = XI_ikl,
                           L = L,
                           K = K,
                           verbose = FALSE,
                           itemp = 1,
                           epsilon = 1e-8,
                           nsim = 500,
                           estimate_itemp = TRUE,
                           hyperpar = hyperpar,
                           itemp_grid =  seq(0,2,by=.1))
  cat(paste("|--------------- Done with iteration:", i ,"-------------|\n"))
  saveRDS(MOD,paste0("CRCC_Analysis/Output/RDS/All_runs_Lipids_Pose/Pose_4kLipids_fisan_K30L40_epsilon1e-8_run#",i,".RDS"))
  rm(MOD)
  gc()
}

## Extract best run ------------------------------------------------------------

ELBOS <- list()

for(i in 1:nsim){
  x <- readRDS(paste0("CRCC_Analysis/Output/RDS/All_runs_Lipids_Pose/Pose_4kLipids_fisan_K30L40_epsilon1e-8_run#",i,".RDS"))
  ELBOS[[i]] <- x$Elbo_val
  cat(i)
}

plot(ELBOS[[1]],ylim=range(unlist(ELBOS)),type="l",xlim=c(0,max(unlist(lapply(ELBOS, function(x) length(x))))))
for(i in 1:nsim){
  points(ELBOS[[i]],type="l")
}
ind.max = which.max(unlist(lapply(ELBOS, max)))
points(ELBOS[[ind.max]],col=2)



mod1 <- readRDS(paste0("CRCC_Analysis/Output/RDS/All_runs_Lipids_Pose/Pose_4kLipids_fisan_K30L40_epsilon1e-8_run#",ind.max,".RDS"))
saveRDS(mod1,"CRCC_Analysis/Output/RDS/All_runs_Lipids_Pose/Best_Pose_4kLipids.RDS")

colcl = apply(mod1$RHO,1,which.max)
plot(pixels,pch=22,bg=colcl,cex=2)

D_hpm <- data.frame(pixels,v = colcl) %>% mutate(v2 = (v))

ggplot(D_hpm)+
  geom_tile(aes(x=pixels[,1],
                y=pixels[,2],fill = as.factor(v)),col=1) +
  scale_fill_viridis_d() +
  theme_void() + theme(legend.position = "bottom") + labs(fill = "cluster")

