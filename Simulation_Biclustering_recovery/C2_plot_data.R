data_bicl_1 <- readRDS("Simulation_Biclustering_recovery/OUTPUT/RDS/DATASETS/data_bicl_1.RDS")

image(data_bicl_1$true_mean[1:120,])
inds <- sort(data_bicl_1$cl_pix,index = TRUE)
image(data_bicl_1$true_mean[,inds$ix])

dim(data_bicl_1$data1)

col <- reshape2::melt(matrix(data_bicl_1$cl_pix,20,20))
xin = c(0.5,cumsum(table(data_bicl_1$cl_pix))+.5)


A <- ggplot(col)+theme_bw()+
  scale_fill_viridis_d("",option="A",direction = -1)+
  geom_tile(aes(x=Var2, y =Var1, fill = factor(value)),col=1)+
#  coord_fixed()+
  facet_wrap(~"CC structure")+
  theme(legend.position = "bottom",
        text = element_text(size=16))+
  xlab("")+ylab("")+
  guides(fill = guide_legend(nrow = 1))

A

data <- data_bicl_1$true_mean[,inds$ix]
md <-  reshape2::melt(data)

B <- ggplot(md)+theme_bw()+
  scale_fill_viridis_c("",option="A",direction = -1)+
  geom_tile(aes(x=Var2, y =Var1, fill =value))+
  xlab("")+ylab("")+
  facet_wrap(~"True values")+
  theme(legend.position = "bottom",
        text = element_text(size=16))+theme(
          # Adjust width of the key (make it a longer rectangle)
          legend.key.width = unit(1.5, "cm"),
          # Optional: adjust height
          legend.key.height = unit(.5, "cm")
        )+
  geom_vline(xintercept = xin,col="gray")

B

data2 <- data_bicl_1$X1[,inds$ix]
md2   <-  reshape2::melt(data2)

C <- ggplot(md2)+theme_bw()+
  scale_fill_viridis_c("",option="A",direction = -1)+
  geom_tile(aes(x=Var2, y =Var1, fill =value))+
  facet_wrap(~"Noisy measurements")+
  xlab("")+ylab("")+
  theme(legend.position = "bottom",
        text = element_text(size=16))+theme(
          # Adjust width of the key (make it a longer rectangle)
          legend.key.width = unit(1.5, "cm"),
          # Optional: adjust height
          legend.key.height = unit(.5, "cm")
        )+
  geom_vline(xintercept = xin,col="gray")


library(patchwork)
Q <- A+B+C

Q + ggview::canvas(h=4,w=12)
ggsave("Simulation_Biclustering_recovery/OUTPUT/PLOT/plot_data.pdf",h=4,w=12)
ggsave("Simulation_Biclustering_recovery/OUTPUT/PLOT/plot_data.png",h=4,w=12)
ggsave("Simulation_Biclustering_recovery/OUTPUT/PLOT/plot_data.eps",h=4,w=12)
