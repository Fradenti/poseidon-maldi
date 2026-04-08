library(tidyverse)
library(raster)
library(gridExtra)
library(RColorBrewer)

source("Simulation_VaryingSpatial_competitors/A_Aux_Function_1.R")
source("Simulation_VaryingSpatial_competitors/A_Aux_Function_2.R")
source("Simulation_VaryingSpatial_competitors/B_Hyperpar.R")

# -------------------------------------------------------------------------

n_side_x <- 20
n_side_y <- 20
n_pix = n_pixel = n_side_x*n_side_y

# coords = expand.grid(x = 1:n_side, y = 1:n_side)

rast <- raster::raster(matrix(1,nrow = n_side_y,ncol = n_side_x))
coords <- as.data.frame(raster::rowColFromCell(rast, 1:raster::ncell(rast)))
coords$cell <- cellFromRowCol(
  rast,
  coords$row,
  coords$col
)
colnames(coords) <- c("y", "x", "cell")
coords <- coords[,c("cell", "x", "y")]


dist_dahl <- as.matrix(stats::dist(coords[,c("x", "y")], method = ifelse(ROOK_vs_QUEEN,"manhattan","maximum")))

## this code generates multiple (Nsample) partitions! one with their own permutation
Nsample <- 10
i_sim <- 8 # 1, 7, 8

plot_list <- list()

for(i_str in seq_along(strength)){
  set.seed(1234)
  perm_matrix <- sapply(1:Nsample, function (x) sample(n_pix)) # one permutation for each column (n x Nsample)
  str <- strength[i_str]
  alpha_DP <- 2; lambda_dahl <- exp(-dist_dahl*str)
  tmp <- apply(perm_matrix, 2, function(perm) sample_dahl(perm,lambda_dahl,alpha_DP) )
  parts <- tmp[1:n_pix,] # if you need the log prior probabilities, they're stored in tmp[n_pix+1,]

  cl_pix <- parts[,i_sim]
  K_pix <- length(unique(cl_pix))

  max_cluster <- 13 # this you need to figure out manually by checking all the partitions visualized
  my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(max_cluster)

  p1 = coords %>% mutate(cluster = cl_pix) %>%
    ggplot() + geom_tile(aes(x=x,y=y,fill = as.factor(cluster)), color = "gray") +
    # scale_fill_brewer(palette = "Set3", name = "Cluster") +
    scale_fill_manual(values = my_colors, name = "Cluster")+
    xlab("")+ ylab("") +
    theme_minimal() + theme(axis.text = element_blank(),
                            axis.ticks = element_blank()) +
    ggtitle(paste0("Strength s = ", str))
  plot_list[[i_str]] <- p1

}

# p <- do.call("grid.arrange", c(plot_list, ncol=2))
# ggsave("Simulation_VaryingSpatial_competitors/plots/EPA_clusters.png", plot = p,
#        width = 16, height = 16, units = "cm", dpi = 300, scale = 1.08)
p <- do.call("grid.arrange", c(plot_list, ncol=4))
ggsave("Simulation_VaryingSpatial_competitors/plots/EPA_clusters.png", plot = p,
       width = 32, height = 8, units = "cm", dpi = 300, scale = 1.08)
