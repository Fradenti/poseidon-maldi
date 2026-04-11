library("spdep")
library("patchwork")
library("tidyverse")
normalizer <- function(x) (x-min(x))/(max(x)-min(x))

# Orginal image -----------------------------------------------------------
library(imager)
library(tidyverse)
library(patchwork)
annA <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(A)"),
  hjustvar = c(1) ,
  vjustvar = c(1))
annB <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(B)"),
  hjustvar = c(1) ,
  vjustvar = c(1))
annC <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(C)"),
  hjustvar = c(1) ,
  vjustvar = c(1))
annD <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(D)"),
  hjustvar = c(1) ,
  vjustvar = c(1))

orig = load.image("CRCC_Analysis/kidney3.png")
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
  geom_tile(aes(x=-x,y=y),fill=D_only$col)+
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
  xlab("x coord.")+  ylab("y coord.")+
  geom_text(data = annA, aes(x=Inf,y=0,hjust=1.5,
                             vjust=0.5,label=annotateText),
            size=7)


Fig

# 4k --------------------------------------------------------------------------
Pep <- readRDS("CRCC_Analysis/Data/4k/FinalMatrix_Peptides.RDS")
Gli <- readRDS("CRCC_Analysis/Data/4k/FinalMatrix_Glycans.RDS")
Lip <- readRDS("CRCC_Analysis/Data/4k/FinalMatrix_Lipids.RDS")
pix <- readRDS("CRCC_Analysis/Data/4k/FinalMatrix_Coor.RDS")
dim(Gli)
dim(Lip)
dim(Pep)

dim(pix)
plot(pix)
# pix <- Rotation(pix, angle=270*pi/180)
pix <- as.data.frame(pix) #as.data.frame(round(pix,1))
colnames(pix) = c("X","Y")


theme_set(theme_bw())
a1 = ggplot()+
  geom_tile(aes(x=-pix$X,
                y=pix$Y,
                fill= normalizer(apply(Pep,1,median))),col=1,lwd=.02)+
  scale_fill_distiller(palette="Spectral")+
  facet_wrap(~"Peptides")+
  theme(#axis.title.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+
  xlab("x coord.")+  ylab("y coord.")+
  theme(legend.position = "none",
        text=element_text(size=18))+
  geom_text(data = annD, aes(x=Inf,y=Inf,hjust=1.5,
                             vjust=1.5,label=annotateText),
            size=7)


a2 = ggplot()+
  geom_tile(aes(x=-pix$X,
                y=pix$Y,
                fill= normalizer(apply(Lip,1,median))),col=1,lwd=.02)+
  scale_fill_distiller(palette="Spectral")+
  facet_wrap(~"Lipids")+
  theme(#axis.title.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+
  xlab("x coord.")+  ylab("y coord.")+
  theme(legend.position = "none",
        text=element_text(size=18))+
  geom_text(data = annB, aes(x=Inf,y=Inf,hjust=1.5,
                             vjust=1.5,label=annotateText),
            size=7)


a3 = ggplot()+
  geom_tile(aes(x=-pix$X,
                y=pix$Y,
                fill= normalizer(apply(Gli,1,median))),col=1,lwd=.02)+
  scale_fill_distiller(palette="Spectral")+
  facet_wrap(~"N-Glycans")+
  theme(#axis.title.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+
  xlab("x coord.")+  ylab("y coord.")+
  theme(legend.position = "none",
        text=element_text(size=18))+
  geom_text(data = annC, aes(x=Inf,y=Inf,hjust=1.5,
                             vjust=1.5,label=annotateText),
            size=7)


Fig + a2 + a3 + a1 +
  plot_layout(nrow=1) +
  ggview::canvas(w=20,h=5)


ggsave("CRCC_Analysis/Output/PLOT/stain_median_v3.pdf",width = 20,height = 5)
ggsave("CRCC_Analysis/Output/PLOT/stain_median_v3.png",width = 20,height = 5)
ggsave("CRCC_Analysis/Output/PLOT/stain_median_v3.eps",width = 20,height = 5)

# COL1 <- Fig + a2 + a3 + a1 +
#   plot_layout(ncol=1)



# -------------------------------------------------------------------------

inds = 125000

annE <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(E)"),
  hjustvar = c(1) ,
  vjustvar = c(1))
annF <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(F)"),
  hjustvar = c(1) ,
  vjustvar = c(1))
annG <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(G)"),
  hjustvar = c(1) ,
  vjustvar = c(1))
annH <- data.frame(
  xpos =  c(Inf),
  ypos =  c(Inf),
  annotateText = c("(H)"),
  hjustvar = c(1) ,
  vjustvar = c(1))




Fig = ggplot(D_only)+
  #  ylim(10,210)+
  theme_bw()+
  geom_tile(aes(x=-x,y=y),fill=D_only$col,alpha=.5)+
  geom_vline(xintercept = -D_only[inds,1], lty=3, col="black")+
  geom_hline(yintercept = D_only[inds,2], lty=3, col="black")+
  geom_point(y = -D_only[inds,2],
             x = -D_only[inds,1],
             col="black",cex=2,pch=0)+
  # geom_tile(data= D_only[inds,],
  #           aes(x=-x,y=y),fill=1,width = 2,height = 2)+
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
  xlab("x coord.")+  ylab("y coord.")+
  geom_text(data = annE, aes(x=Inf,y=0,hjust=1.5,
                             vjust=0.5,label=annotateText),
            size=7)

Fig

plot(pix)
points(20,55,bg=2,pch=21)

which_x = which(pix$X==20)
which_y = which(pix$Y==55)


ind <- pix %>% mutate(n = 1:nrow(pix)) %>% filter(X==20,Y==55) %>% select(n)
ind <- pull(ind)
p = t(Pep)[,ind]
l = t(Lip)[,ind]
g = t(Gli)[,ind]

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



a1 = ggplot()+
  geom_point(aes(y = p,
              x = allx__p))+
  geom_segment(aes(y = p,
                   yend = 0,
                   x = allx__p,
                   xend =  allx__p))+
  scale_x_discrete("m/z",breaks = limx__p )+
  geom_hline(yintercept = 0, col="gray")+
  facet_wrap(~"Peptides")+
  coord_cartesian(ylim=c(0,max(c(p,l,g))))+
  xlab("m/z")+  ylab("Abundance")+
  theme(legend.position = "none",
        text=element_text(size=18),
        axis.text.x = element_text(size=10))+
  geom_text(data = annH, aes(x=Inf,y=Inf,hjust=1.5,
                             vjust=1.5,label=annotateText),
            size=7)
a1
a2 = ggplot()+
  geom_point(aes(y = l,
                 x = allx__l))+
  geom_segment(aes(y = l, yend=0,
                   x = allx__l))+
  geom_hline(yintercept = 0, col="gray")+
  facet_wrap(~"Lipids")+
  coord_cartesian(ylim=c(0,max(c(p,l,g))))+
  scale_x_discrete("m/z",breaks = limx__l )+
  xlab("m/z")+  ylab("Abundance")+
  theme(legend.position = "none",
        text=element_text(size=18),
        axis.text.x = element_text(size=10))+
  geom_text(data = annF, aes(x=Inf,y=Inf,hjust=1.5,
                             vjust=1.5,label=annotateText),
            size=7)

a3 = ggplot()+
  geom_point(aes(y = g,
                 x = allx__g))+
  geom_segment(aes(y = g, yend=0,
                   x = allx__g))+
  geom_hline(yintercept = 0, col="gray")+
  facet_wrap(~"N-Glycans")+
  coord_cartesian(ylim=c(0,max(c(p,l,g))))+
  scale_x_discrete("m/z",breaks = limx__g )+
  ylab("Abundance")+
  theme(legend.position = "none",
        text=element_text(size=18),
        axis.text.x = element_text(size=10))+
  geom_text(data = annG, aes(x=Inf,y=Inf,hjust=1.5,
                             vjust=1.5,label=annotateText),
            size=7)




Fig + a2 + a3 + a1 +
  plot_layout(nrow=1) +
  ggview::canvas(w=20,h=5)


ggsave("CRCC_Analysis/Output/PLOT/stain_spectr_v3.pdf",width = 20,height = 5)
ggsave("CRCC_Analysis/Output/PLOT/stain_spectr_v3.png",width = 20,height = 5)
ggsave("CRCC_Analysis/Output/PLOT/stain_spectr_v3.eps",width = 20,height = 5,
       device = cairo_ps)


COL2 <- Fig + a2 + a3 + a1 +
  plot_layout(ncol=1)

(COL1|COL2)+
  ggview::canvas(w=15,h=20)


ggsave("CRCC_Analysis/Output/PLOT/Columns_v3.pdf",width = 10,height = 20)
ggsave("CRCC_Analysis/Output/PLOT/Columns_v3.png",width = 10,height = 20)
ggsave("CRCC_Analysis/Output/PLOT/Columns_v3.eps",width = 10,height = 20,
       device = cairo_ps)
