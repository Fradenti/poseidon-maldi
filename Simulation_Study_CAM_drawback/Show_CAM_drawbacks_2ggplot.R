library(tidyverse )

dmixnorm = function(x,pi,ucl,Res,k){
  sum(pi[,ucl[k]] * dnorm(x, Res[,1], sqrt(Res[,2])))
}

draw_dens_pose_gg <- function(res,
                              k = 1,
                              xlim = NULL,
                              n_points_dens = 5000,
                              hist_breaks = 25) {

  # 1. Data Preparation (Logic remains from your original function)
  if (is.null(xlim)) {
    xlim <- c(min(res$Y), max(res$Y))
  }

  clu <- estimate_clust_pose(res, bicl_mat = FALSE)
  ucl <- unique(clu$col_cl)
  Res <- post_exp_val(res)
  L   <- nrow(Res)

  if (res$model %in% c("fiSAN", "fSAN")) {
    pi <- apply(res$B_star, 2, function(t) t / sum(t))
  } else if (res$model == "CAM") {
    S  <- res$a_bar_lk + res$b_bar_lk
    pi <- res$a_bar_lk / S * apply(rbind(1, (res$b_bar_lk / S)[-L, ]), 2, cumprod)
  }

  # 2. Extract specific data for the k-th cluster
  y_vals <- res$Y[, clu$col_cl == ucl[k]]
  df_hist <- data.frame(y = y_vals)

  # 3. Create Density Curve Data
  SEQ <- seq(xlim[1], xlim[2], length.out = n_points_dens)
  yy  <- sapply(SEQ, function(x) dmixnorm(x = x, pi = pi, ucl = ucl, Res = Res, k = k))
  df_dens <- data.frame(x = SEQ, y = yy)

  # 4. Create Components Data (The points and lines at the bottom)
  # Normalizing Res[,3] for point size as you did with 'norm'
  df_comp <- data.frame(
    mu = Res[,1],
    pi_k = pi[, ucl[k]] ) %>% filter(pi_k > 1e-2)
  df_hist <- as_tibble(as_tibble(unlist(df_hist)))
  # 5. Build the ggplot
  ggplot() +
    # Histogram layer
    geom_histogram(data = df_hist, aes(x = value, y = after_stat(density)),
                   bins = hist_breaks, fill = "lightgray", color = "darkblue",alpha=.5,lwd=.1) +
    # Density line layer
    geom_line(data = df_dens, aes(x = x, y = y),
              color = "royalblue3", linewidth = 0.5) +
    # Lollipop plot for mixture components
    geom_segment(data = df_comp, aes(x = mu, xend = mu, y = 0, yend = df_comp$pi_k),
                 color = "black") +
    geom_point(data = df_comp, aes(x = mu, y = pi_k), size = df_comp$pi_k,
               color = "black") +
    # Formatting
    coord_cartesian(xlim = xlim) +
    facet_wrap(~paste0("CC #", k))+
    labs(x = "Y", y = "Posterior density") +
    theme_bw() +
    theme(text = element_text(size=18))+
    guides(size = "none") # Hide the legend for the point sizing
}

# -------------------------------------------------------------------------
res_cam <- readRDS("Simulation_Study_CAM_drawback/Output/RDS/run_CAM_30obs.RDS")
res_fisan <- readRDS("Simulation_Study_CAM_drawback/Output/RDS/run_fiSAN_30obs.RDS")
best_cam <- extract_best_run_from_list(res_list = res_cam)
best_fisan <- extract_best_run_from_list(res_list = res_fisan)


a1 <- draw_dens_pose_gg(res = best_cam,k = 1,n_points_dens = 400,hist_breaks = 25)
a2 <- draw_dens_pose_gg(res = best_cam,k = 2,n_points_dens = 400,hist_breaks = 25)
a3 <- draw_dens_pose_gg(res = best_cam,k = 3,n_points_dens = 400,hist_breaks = 25)
library(patchwork)
Camplot <- a1+a2+a3+
  plot_annotation(
    title = 'CC posterior density estimate with common atoms',
#    subtitle = 'Comparison between Cluster 1 and Cluster 2',
#    caption = 'Source: Posterior Predictive Results',
    theme = theme(plot.title = element_text(size = 18, hjust = 0))
  )
Camplot+ ggview::canvas(h=5,w=15)

ggsave("Simulation_Study_CAM_drawback/Output/PLOT/dens_cam.png",h=5,w=15)
ggsave("Simulation_Study_CAM_drawback/Output/PLOT/dens_cam.pdf",h=5,w=15)

a1 <- draw_dens_pose_gg(res = best_fisan,k = 1,n_points_dens = 400,hist_breaks = 25)
a2 <- draw_dens_pose_gg(res = best_fisan,k = 2,n_points_dens = 400,hist_breaks = 25)
a3 <- draw_dens_pose_gg(res = best_fisan,k = 3,n_points_dens = 400,hist_breaks = 25)
library(patchwork)
fisplot <- a1+a2+a3+
  plot_annotation(
    title = 'CC posterior density estimate with shared atoms',
    #    subtitle = 'Comparison between Cluster 1 and Cluster 2',
    #    caption = 'Source: Posterior Predictive Results',
    theme = theme(plot.title = element_text(size = 18, hjust = 0))
  )
fisplot+ ggview::canvas(h=5,w=15)
ggsave("Simulation_Study_CAM_drawback/Output/PLOT/dens_fisan.png",h=5,w=15)
ggsave("Simulation_Study_CAM_drawback/Output/PLOT/dens_fisan.pdf",h=5,w=15)

