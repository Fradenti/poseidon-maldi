source("CRCC_Analysis/000_Auxiliary.R")

library(tidyverse)
library(desplot)
library(ggsci)
theme_set(theme_bw())
library(raster)
library(patchwork)
library(ggsci)


# Read Data ---------------------------------------------------------------
K <- 30
L <- 40
pixels     <- (readRDS("CRCC_Analysis/Data/FinalMatrix_Coor.RDS"))
nnList_cpp <- readRDS("CRCC_Analysis/Data/list_ind_resolution_4123_cpp.RDS")
Gli <- val1        <- readRDS("CRCC_Analysis/Data/FinalMatrix_Glycans.RDS")
Lip <- val2        <- readRDS("CRCC_Analysis/Data/FinalMatrix_Lipids.RDS")
Pep <- val3        <- readRDS("CRCC_Analysis/Data/FinalMatrix_Peptides.RDS")

####
ind <- pixels %>% mutate(n = 1:nrow(pixels)) %>%
  filter(X==50,Y==25) %>% dplyr::select(n)
ind <- pull(ind)
p <- t(Pep)[,ind]
l <- t(Lip)[,ind]
g <- t(Gli)[,ind]

p <- t(Pep)[,ind]
l <- t(Lip)[,ind]
g <- t(Gli)[,ind]

ind__g <- round(seq(from=5,to=length(names(g))-5,length.out=5))
ind__p <- round(seq(from=5,to=length(names(p))-5,length.out=5))
ind__l <- round(seq(from=5,to=length(names(l))-5,length.out=5))

rou__g <- round(as.numeric(names(g)[ind__g]),2)
rou__p <- round(as.numeric(names(p)[ind__p]),2)
rou__l <- round(as.numeric(names(l)[ind__l]),2)

allx__g <- factor(round(as.numeric(names(g)),2))
allx__p <- factor(round(as.numeric(names(p)),2))
allx__l <- factor(round(as.numeric(names(l)),2))

limx__g <- factor(rou__g)
limx__p <- factor(rou__p)
limx__l <- factor(rou__l)

mz1 <- colnames(val1)
mz2 <- colnames(val2)
mz3 <- colnames(val3)

colnames(val1) <- rownames(val1) <- NULL
colnames(val2) <- rownames(val2) <- NULL
colnames(val3) <- rownames(val3) <- NULL

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

r1 <- mapper2(t(val1))
r2 <- mapper2(t(val2))
r3 <- mapper2(t(val3))

all <- list(r1,r2,r3)

