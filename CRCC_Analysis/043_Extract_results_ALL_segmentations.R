library(tidyverse)
library(Poser)
library(desplot)
library(patchwork)
library(POSEIDON)



# REAL IMAGE-------------------------------------------------------------------------

library("spdep")
library("patchwork")
library("tidyverse")
library(imager)
normalizer <- function(x) (x-min(x))/(max(x)-min(x))
orig = load.image("CRCC_Analysis/S1_ccrcc_legend.png")

sc <- 1
dim(orig)
down <- resize(orig,round( width(orig)/sc),
               round(height(orig)/sc))

dim(down)
prod(dim(down)[1:2])


Dr <- as.data.frame(down) %>% filter(cc==1) %>% dplyr::select(-cc)
Dg <- as.data.frame(down) %>% filter(cc==2) %>% dplyr::select(-cc)
Db <- as.data.frame(down) %>% filter(cc==3) %>% dplyr::select(-cc)
Dt <- as.data.frame(down) %>% filter(cc==4) %>% dplyr::select(-cc)

D = Dr %>% rename(r = value) %>% mutate(g = Dg$value, b = Db$value) %>% filter(Dt$value>0)
#kme = kmeans(D[,3:5],10)
#plot(down)

#kme <- kmeans(D[,3:5],10)
#kme$centers
#ind <- which.min(kme$centers[,1])

D_only <- D %>% #mutate(cl = kme$cluster) %>% filter(cl != ind) %>%
  mutate(col = rgb(r,g,b))
D_only


Fig = ggplot(D_only)+
  #  ylim(10,210)+
  theme_bw()+
  geom_tile(aes(x=x,y=y),fill=D_only$col)+
  scale_y_reverse()+
  facet_wrap(~"H&E stain image")+
  theme(#axis.title.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+
  theme(legend.position = "bottom", text=element_text(size = 18))+
  theme(legend.position = "none")+# +
  xlab("x coord.")+  ylab("y coord.")

# -------------------------------------------------------------------------















pixels <- readRDS("CRCC_Analysis/Data/pixels.RDS")
out3   <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Best_Poseidon.RDS")
outA   <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Stacked_Pose/Best_Pose_4kStacked.RDS")
outL   <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Lipids_Pose/Best_Pose_4kLipids.RDS")
outG   <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Glycans_Pose/Best_Pose_4kGlycans.RDS")
outP   <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Peptides_Pose/Best_Pose_4kPeptides.RDS")


# -------------------------------------------------------------------------
times_poseidon <- numeric(250)
times_poserall <- numeric(250)
for(j in 1:250){
  a <- readRDS(paste0("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Poseidon_4k_fisan_K30L40_epsilon1e-8_run#",j,".RDS"))
  b <- readRDS(paste0("CRCC_Analysis/Output/RDS/All_runs_Stacked_Pose/Pose_4kStacked_fisan_K30L40_epsilon1e-8_run#",j,".RDS"))
  times_poseidon[j] <- as.numeric(a$time,units="secs")
  times_poserall[j] <- as.numeric(b$time,units="secs")
}
boxplot(cbind(times_poseidon,times_poserall))
mean(times_poseidon)/60
sd(times_poseidon)/60
# -------------------------------------------------------------------------
clA <- estimate_clust_pose(outA)
clL <- estimate_clust_pose(outL)
clG <- estimate_clust_pose(outG)
clP <- estimate_clust_pose(outP)
cl3 <- estimate_clust_poseidon(out3)

length(table(clA$col_cl))
summary(as.numeric(table(clA$col_cl)))
length(table(cl3$col_cl))
summary(as.numeric(table(cl3$col_cl)))

table(clL$col_cl)
table(clG$col_cl)
table(clP$col_cl)
table(cl3$col_cl)

mean(((outA$Y) - clA$bicl)^2)

mean(c(
mean((out3$Y_list[[1]] - cl3$bicl[[1]])^2),
mean((out3$Y_list[[2]] - cl3$bicl[[2]])^2),
mean((out3$Y_list[[3]] - cl3$bicl[[3]])^2)
))
mean((outL$Y - clL$bicl)^2)
mean((outG$Y - clG$bicl)^2)
mean((outP$Y - clP$bicl)^2)

image(clA$bicl)

inds <- apply(clA$row_cl,2,function(x) sort(x,index=TRUE))
SORTED <- c()
sa <- sort(unique(clA$col_cl))
for(i in seq_along(sa)){

  SORTED <- cbind(SORTED,outA$Y[inds[[i]]$ix,clA$col_cl == sa[i]])

}
image(SORTED)

# -------------------------------------------------------------------------



fun <- function(out, title){

  if(title == "(E) Poseidon"){

    clustering <- estimate_clust_poseidon(out)
    table(clustering$col_cl)
    D <- colMeans(do.call(rbind,out$Y)) %>% as_tibble() %>% mutate(cl = clustering$col_cl) %>%
      group_by(cl) %>% mutate(m = mean(value)) %>% ungroup() %>%
      mutate(colf = as.numeric(factor(m)))
    Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = factor(D$colf))

  }else{

    clustering <- estimate_clust_pose(out)
    inds <- sort(tapply(colMeans(out$Y), clustering$col_cl, mean),index = T)
    ord_cl <- clustering$col_cl[inds$ix]

    D <- colMeans(out$Y) %>% as_tibble() %>% mutate(cl = clustering$col_cl) %>%
      group_by(cl) %>% mutate(m = mean(value)) %>% ungroup() %>%
      mutate(colf = as.numeric(factor(m)))
    Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = factor(D$colf))

  }


  Grd <- expand_grid(x=0:(max(pixels[,1])+1),
                     y=0:(max(pixels[,2])+1)) %>%
    mutate(p = paste0(x,"-",y)) %>% left_join(Px,by = "p")

Grd <- Grd %>% mutate(ccl2 = as.factor(ifelse(is.na(ccl),0,ccl))) %>% mutate(title = title)
palette <- colorRampPalette(c("darkblue", "white", "tomato3"))(length(table(Grd$ccl2)))

POS = ggplot(Grd)+
  geom_tile(aes(x=-X,
                y=Y,
                fill = ccl2),
            col="gray",lwd=.1,alpha=.5) +
  geom_tileborder(
    aes(x=-x,
        y=y,
        group=1,
        grp=ccl2), lwd=.5)+
  scale_fill_viridis_d(option = "B",direction = 1)+
  facet_wrap(~title) +
  xlab("x coord.") + ylab("y coord.") +
  theme_void()+
  theme(legend.position = "none",text = element_text(size=15))
POS
}

LI <- fun(out = outL,title = "(A) Lipids")+coord_fixed()
GL <- fun(out = outG,title = "(B) N-Glycans")+coord_fixed()
PE <- fun(out = outP,title = "(C) Peptides")+coord_fixed()
AA <- fun(out = outA,title = "(D) Pose - Stacked datasets")+coord_fixed()
PP <- fun(out = out3,title = "(E) Poseidon")+coord_fixed()

FF <-  Fig +  theme_void()+
  theme(legend.position = "none",text = element_text(size=15))+
  facet_wrap(~"(F) Reference") + coord_fixed(ratio=1.558282)

FF



(LI|GL|PE)/(AA|PP) + ggview::canvas(h=7,w=10)

ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v2.pdf", h=7, w=10)
ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v2.png", h=7, w=10)
ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v2.eps", h=7, w=10, device = cairo_ps)









# -------------------------------------------------------------------------






LI <- fun(out = outL,title = "(A) Lipids")+coord_fixed()
GL <- fun(out = outG,title = "(B) N-Glycans")+coord_fixed()
PE <- fun(out = outP,title = "(C) Peptides")+coord_fixed()
AA <- fun(out = outA,title = "(D) Pose - Stacked datasets")+coord_fixed()
PP <- fun(out = out3,title = "(E) Poseidon")+coord_fixed()

FF <-  Fig +  theme_void()+
  theme(legend.position = "none",text = element_text(size=15))+
  facet_wrap(~"(F) Reference") + coord_fixed(ratio=1.558282)

FF



(LI|GL|PE)/(AA|PP|FF) + ggview::canvas(h=7,w=10)

ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v31.pdf", h=7, w=10)
ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v31.png", h=7, w=10)
ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v31.eps", h=7, w=10, device = cairo_ps)



# -------------------------------------------------------------------------




fun2 <- function(out, title){

  if(title == "(E) Poseidon"){

    clustering <- estimate_clust_poseidon(out)
    table(clustering$col_cl)
    D <- colMeans(do.call(rbind,out$Y)) %>% as_tibble() %>% mutate(cl = clustering$col_cl) %>%
      group_by(cl) %>% mutate(m = mean(value)) %>% ungroup() %>%
      mutate(colf = as.numeric(factor(m)))
    Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = factor(D$colf))

  }else{

    clustering <- estimate_clust_pose(out)
    inds <- sort(tapply(colMeans(out$Y), clustering$col_cl, mean),index = T)
    ord_cl <- clustering$col_cl[inds$ix]

    D <- colMeans(out$Y) %>% as_tibble() %>% mutate(cl = clustering$col_cl) %>%
      group_by(cl) %>% mutate(m = mean(value)) %>% ungroup() %>%
      mutate(colf = as.numeric(factor(m)))
    Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = factor(D$colf))

  }


  Grd <- expand_grid(x=0:(max(pixels[,1])+1),
                     y=0:(max(pixels[,2])+1)) %>%
    mutate(p = paste0(x,"-",y)) %>% left_join(Px,by = "p")

  Grd <- Grd %>% mutate(ccl2 = as.factor(ifelse(is.na(ccl),0,ccl))) %>% mutate(title = title)
  palette <- colorRampPalette(c("darkblue", "white", "tomato3"))(length(table(Grd$ccl2)))

  POS = ggplot(Grd)+
    geom_tile(aes(x=-X,
                  y=Y,
                  fill = ccl2),
              col="gray",lwd=.1,alpha=.5) +
    geom_tileborder(
      aes(x=-x,
          y=y,
          group=1,
          grp=ccl2), lwd=.5)+
    scale_fill_viridis_d(option = "B",direction = 1)+
    facet_wrap(~title) +
    xlab("x coord.") + ylab("y coord.") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 0.8),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",text = element_text(size=15),
          aspect.ratio = 1)
  POS
}

LI <- fun2(out = outL,title = "(A) Lipids")+coord_fixed(ratio = 1)
GL <- fun2(out = outG,title = "(B) N-Glycans")+coord_fixed(ratio = 1)
PE <- fun2(out = outP,title = "(C) Peptides")+coord_fixed(ratio = 1)
AA <- fun2(out = outA,title = "(D) Pose - Stacked datasets")+coord_fixed(ratio = 1)
PP <- fun2(out = out3,title = "(E) Poseidon")+coord_fixed(ratio = 1)

FF <-  Fig +    theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )+
  theme(legend.position = "none",text = element_text(size=15),
        aspect.ratio = 1)+
  facet_wrap(~"(F) Reference")


(LI|GL|PE)/(AA|PP|FF) + ggview::canvas(h=8,w=12)

ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v31.pdf", h=8, w=10)
ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v31.png", h=8, w=10)
ggsave("CRCC_Analysis/Output/PLOT/Model_results/ALL_segmentations_v31.eps", h=8, w=10, device = cairo_ps)
