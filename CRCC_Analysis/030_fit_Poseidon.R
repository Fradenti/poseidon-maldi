# -------------------------------------------------------------------------

library(tidyverse)
library(raster)
library(Poser)
library(POSEIDON)
library(patchwork)

K <- 30
L <- 40
pixels     <- readRDS("CRCC_Analysis/Data/FinalMatrix_Coor.RDS")
nnList_cpp <- readRDS("CRCC_Analysis/Data/list_ind_resolution_4123_cpp.RDS")

val1        <- readRDS("CRCC_Analysis/Data/FinalMatrix_Glycans.RDS")
val2        <- readRDS("CRCC_Analysis/Data/FinalMatrix_Lipids.RDS")
val3        <- readRDS("CRCC_Analysis/Data/FinalMatrix_Peptides.RDS")

colnames(val1) = rownames(val1) = NULL
colnames(val2) = rownames(val2) = NULL
colnames(val3) = rownames(val3) = NULL

# check NN are ok ---------------------------------------------------------

plot(pixels)
points(pixels[1,],col=2,pch=21,bg=4)
points(pixels[nnList_cpp[[1]]+1,],col=2,pch=21,bg="1")

points(pixels[134,],col=2,pch=21,bg=2)
points(pixels[nnList_cpp[[134]]+1,],col=2,pch=21,bg="1")

points(pixels[1134,],col=2,pch=21,bg=2)
points(pixels[nnList_cpp[[1134]]+1,],col=2,pch=21,bg="1")

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

r1 = mapper2(t(val1))
r2 = mapper2(t(val2))
r3 = mapper2(t(val3))

r = list(r1,r2,r3)
sum(is.na(r))




# Multistart --------------------------------------------------------------
nsim <- 250
for(i in 1:nsim){
# intial values -----------------------------------------------------------
  set.seed(123143*i)
  par <- maotai::kmeanspp(t(rbind(r1,r2,r2)),
                        k = K)

  RHOi <- Poser::from_partition_to_probs(par,#sample(1:K,size = J,replace = T),
                                       K = K,
                                       perturb = .5)
  XXX <- list()
  for(tt in 1:3){
  N <- nrow(r[[tt]])
  XI_ikl <- array(NA,c(N,K,L))


  for(k in 1:K){
    par_row <- maotai::kmeanspp(as.matrix(r[[tt]][,par==k]),k=5)
    XI_ikl[,k,] = Poser::from_partition_to_probs(par_row, K=L, perturb = .2)
  }
  XXX[[tt]] = XI_ikl
}

#  spatial ---------------------------------------------------------------------

MOD = POSEIDON::POSEIDON_fiSAN(Y_list = r,
                                   NN_inds = nnList_cpp,
                                   RHOi = RHOi,
                                   XIi = XXX,
                                   L = L,
                                   K = K,
                                   verbose = FALSE,
                                   itemp = 1,
                                   epsilon = 1e-8,
                                   nsim = 500,
                                   estimate_itemp = TRUE,
                                   hyperpar = hyperpar,
                                   itemp_grid = seq(0,2,by=.1))

  cat(paste("|--------------- Done with iteration:", i ,"-------------|\n"))
  saveRDS(MOD,paste0("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Poseidon_4k_fisan_K30L40_epsilon1e-8_run#",i,".RDS"))
  rm(MOD)
  gc()
}
## ------------------------------------
ELBOS <- list()
for(i in 1:nsim){
  x <- readRDS(paste0("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Poseidon_4k_fisan_K30L40_epsilon1e-8_run#",i,".RDS"))
  ELBOS[[i]] <- x$Elbo_val
  cat(i)
}

plot(ELBOS[[1]],ylim=range(unlist(ELBOS)),type="l",xlim=c(0,max(unlist(lapply(ELBOS, function(x) length(x))))))
for(i in 1:nsim){
  points(ELBOS[[i]],type="l")
}
ind.max = which.max(unlist(lapply(ELBOS, max)))
points(ELBOS[[ind.max]],col=2)



mod1 <- readRDS(paste0("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Poseidon_4k_fisan_K30L40_epsilon1e-8_run#",ind.max,".RDS"))
saveRDS(mod1,"CRCC_Analysis/Output/RDS/All_runs_Poseidon/Best_Poseidon.RDS")


colcl = apply(mod1$RHO,1,which.max)
plot(pixels,pch=22,bg=colcl,cex=1)
plot(pixels,pch=21,bg=clustering$col_cl/2,cex=1)

D_hpm <- data.frame(pixels,v = colcl) %>% mutate(v2 = (v))

ggplot(D_hpm)+
  geom_tile(aes(x=pixels[,1],
                y=pixels[,2],fill = as.factor(v)),col=1) +
  scale_fill_discrete() +
  theme_void() + theme(legend.position = "bottom") + labs(fill = "cluster")


