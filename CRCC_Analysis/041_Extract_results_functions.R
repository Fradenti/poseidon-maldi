# POSEIDON functions for results -------------------------------------------------------------------------
IS <- function(out, ind, pixels){

  clustering <- relabeling_poseidon(pixels = pixels,res = out)

  Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = factor(clustering$col_cl))
  Grd <- expand_grid(x=0:(max(pixels[,1])+1),
                     y=0:(max(pixels[,2])+1)) %>%
    mutate(p = paste0(x,"-",y)) %>% left_join(Px,by = "p")


  Q <- ggplot(Grd)+
    geom_tile(data=Grd,aes(x = -X,
                           y = Y),alpha=.4,
              col="gray",fill="lightgray")+
    geom_tile(data = Grd %>% filter(ccl == ind),
              aes(x = -X,
                  y = Y), fill= "darkgray",
              alpha=.9, col=1)+
    scale_fill_brewer(palette = "Spectral") +
    facet_wrap(~paste0("Poseidon - Image Segmentation - CC-",ind)) +
    theme(legend.position = "none",text = element_text(size=16)) +
    xlab("x coord.") + ylab("y coord.")

  return(Q)
}


gg_draw_dens_poseidon <- function(res,
                                  pixels,
                                  k = 1,
                                  xlim = NULL,
                                  n_points_dens = 5000,
                                  hist_breaks = 25,
                                  molec = c("N-Glycans", "Lipids", "Peptides")) {

  # 1. Data Preparation (Logic remains from your original function)
  if (is.null(xlim)) {
    xlim <- lapply(res$Y_list, function(Y) c(min(Y), max(Y)))
  }

  clu <- relabeling_poseidon(res = res, pixels = pixels)
  ucl <- as.numeric(clu$Data[,1])
  L   <- dim(res$XI[[1]])[3]

  PLOTS <- list()
  Q <- IS(out = out,ind = k,pixels = pixels)

  for(j in 1:dim(res$MKCD)[3]){

    Res <-  cbind(mu = res$MKCD[,1,j], sigma2 = res$MKCD[,4,j]/(res$MKCD[,3,j]-1))

    if (res$model %in% c("fiSAN", "fSAN")) {
      pi <- apply(res$B_star[,,j], 2, function(t) t / sum(t))
    } else if (res$model == "CAM") {
      # S  <- res$a_bar_lk + res$b_bar_lk
      # pi <- res$a_bar_lk / S * apply(rbind(1, (res$b_bar_lk / S)[-L, ]), 2, cumprod)
      stop("no ggplot for cam version yet")
    }

    # 2. Extract specific data for the k-th cluster
    y_vals <- res$Y_list[[j]][, clu$col_cl == k]
    df_hist <- data.frame(y = y_vals)

    # 3. Create Density Curve Data
    SEQ <- seq(xlim[[j]][1], xlim[[j]][2], length.out = n_points_dens)
    yy  <- sapply(SEQ, function(x) dmixnorm(x = x, pi = pi, ucl = ucl, Res = Res, k = k))
    df_dens <- data.frame(x = SEQ, y = yy)

    # 4. Create Components Data (The points and lines at the bottom)
    # Normalizing Res[,3] for point size as you did with 'norm'
    df_comp <- data.frame(
      mu = Res[,1],
      pi_k = pi[, ucl[k]] ) %>%
      filter(pi_k > 1e-2)
    df_hist <- as_tibble(as_tibble(unlist(df_hist)))
    # 5. Build the ggplot
    P <- local({
      # --- CAPTURE CURRENT DATA VALUES ---
      # These assignments force evaluation at the current loop index 'j'
      cur_df_hist <- df_hist
      cur_df_dens <- df_dens
      cur_df_comp <- df_comp
      # You also need to capture 'j' if you use it in labels, text, etc.
      # cur_j <- j # Not strictly needed for your current plot code, but good practice.
      cur_j <- j

      ggplot() +
        # Histogram layer - USE THE CAPTURED VARIABLE
        geom_histogram(data = cur_df_hist, aes(x = value, y = after_stat(density)),
                       bins = hist_breaks, fill = "lightgray",
                       color = "darkblue",alpha=.5,lwd=.1) +
        # Density line layer - USE THE CAPTURED VARIABLE
        geom_line(data = cur_df_dens, aes(x = x, y = y),
                  color = "royalblue3", linewidth = 0.5) +
        # Lollipop plot layers - USE THE CAPTURED VARIABLE
        geom_segment(data = cur_df_comp, aes(x = mu, xend = mu, y = 0, yend = pi_k), # Removed cur_df_comp$ prefix for yend
                     color = "black") +
        geom_point(data = cur_df_comp, aes(x = mu, y = pi_k), size = cur_df_comp$pi_k,
                   color = "black") +
        # Formatting
        coord_cartesian(xlim = xlim[[cur_j]]) +
        # Use the captured k here
        labs(x = "Transformed abundances", y = "Posterior density") +
        theme_bw() +
        theme(text = element_text(size=18))+
        guides(size = "none")+ facet_wrap(~paste0("CC-", k, " - ", molec[cur_j] ))
    })
    PLOTS[[j]] <- P
  }

  G0 <- PLOTS[[2]]|PLOTS[[1]]|PLOTS[[3]]  # ORDER in PLOT IS LIP-NGLYC-PEP
  G  <- Q+G0
  G <-  G + plot_layout(widths = c(1, 3))
  return(G)
}





gg_draw_dots <- function(pixels, out, ind,  median = T, allx, limx){

  clustering <- relabeling_poseidon(pixels = pixels,res = out)

  n <- length(ind)
  cl_fact <- factor(clustering$col_cl)

  L1 = c()
  L2 = c()
  L3 = c()

  if(median){

    L1 = apply(out$Y_list[[1]][,cl_fact==ind],1,median)
    L2 = apply(out$Y_list[[2]][,cl_fact==ind],1,median)
    L3 = apply(out$Y_list[[3]][,cl_fact==ind],1,median)

    ss <-   ylab("Medians of transformed abundances")

  }else{
    L1 = rowMeans(out$Y_list[[1]][,cl_fact == ind])
    L2 = rowMeans(out$Y_list[[2]][,cl_fact == ind])
    L3 = rowMeans(out$Y_list[[3]][,cl_fact == ind])

    ss <-   ylab("Means of transformed abundances")

  }

  ind_col <- which(colnames(clustering$row_cl[[1]]) %in% paste0("COL_CL#",ind))
  glis <- c(clustering$row_cl[[1]][,ind_col])
  ind_col <- which(colnames(clustering$row_cl[[2]]) %in% paste0("COL_CL#",ind))
  lips <- c(clustering$row_cl[[2]][,ind_col])
  ind_col <- which(colnames(clustering$row_cl[[3]]) %in% paste0("COL_CL#",ind))
  peps <- c(clustering$row_cl[[3]][,ind_col])


  LL1 <- as.data.frame(L1) %>% mutate(xs = allx[[1]]) %>% reshape2::melt() %>% mutate(type = "N-Glycans",  rc = (glis))
  LL2 <- as.data.frame(L2) %>% mutate(xs = allx[[2]]) %>% reshape2::melt() %>% mutate(type = "Lipids",   rc = (lips))
  LL3 <- as.data.frame(L3) %>% mutate(xs = allx[[3]]) %>% reshape2::melt() %>% mutate(type = "Peptides", rc = (peps))


  LLD <- bind_rows(LL1,LL2,LL3) %>% group_by(type,rc) %>% mutate(mv = mean(value))


  library(RColorBrewer)

  # Base palette (max 11 colors from RColorBrewer)
  base_spectral <- brewer.pal(11, "Spectral")
  # Extend to 20 colors using interpolation
  spectral_extended <- colorRampPalette(base_spectral)(length(table(LLD$mv)))

  Q <- IS(out = out,ind = ind,pixels = pixels)
  S0 <- ggplot(LLD)+
    geom_segment(aes(x = xs,xend = xs,
                     y=0,yend = value,
                     group = variable))+
    geom_point(aes(x=xs,y=value,group = variable,
                   fill = cut(mv,length(table(LLD$mv)))),
               size=2,pch=21,col=1)+
    facet_grid(~type, scales= "free_x")+
    ylab(ss)+
    scale_x_discrete("m/z",breaks = c(limx[[1]],limx[[2]],limx[[3]]) )+
    theme_bw()+
    scale_fill_manual(values = spectral_extended)+
    theme(legend.position = "none",text = element_text(size=16))

  G <- Q + S0
  G <- G + plot_layout(widths = c(1, 3))
  return(G)
}




gg_draw_dots_detailed_axis <- function(pixels, out, ind,  median = T, allx, limx){

  clustering <- relabeling_poseidon(pixels = pixels,res = out)

  n <- length(ind)
  cl_fact <- factor(clustering$col_cl)

  L1 = c()
  L2 = c()
  L3 = c()

  if(median){

    L1 = apply(out$Y_list[[1]][,cl_fact==ind],1,median)
    L2 = apply(out$Y_list[[2]][,cl_fact==ind],1,median)
    L3 = apply(out$Y_list[[3]][,cl_fact==ind],1,median)

    ss <-   ylab("Medians of transformed abundances")

  }else{
    L1 = rowMeans(out$Y_list[[1]][,cl_fact == ind])
    L2 = rowMeans(out$Y_list[[2]][,cl_fact == ind])
    L3 = rowMeans(out$Y_list[[3]][,cl_fact == ind])

    ss <-   ylab("Means of transformed abundances")

  }

  ind_col <- which(colnames(clustering$row_cl[[1]]) %in% paste0("COL_CL#",ind))
  glis <- c(clustering$row_cl[[1]][,ind_col])
  ind_col <- which(colnames(clustering$row_cl[[2]]) %in% paste0("COL_CL#",ind))
  lips <- c(clustering$row_cl[[2]][,ind_col])
  ind_col <- which(colnames(clustering$row_cl[[3]]) %in% paste0("COL_CL#",ind))
  peps <- c(clustering$row_cl[[3]][,ind_col])


  LL1 <- as.data.frame(L1) %>% mutate(xs = allx[[1]]) %>% reshape2::melt() %>% mutate(type = "N-glycans",  rc = (glis))
  LL2 <- as.data.frame(L2) %>% mutate(xs = allx[[2]]) %>% reshape2::melt() %>% mutate(type = "Lipids",   rc = (lips))
  LL3 <- as.data.frame(L3) %>% mutate(xs = allx[[3]]) %>% reshape2::melt() %>% mutate(type = "Peptides", rc = (peps))


  LLD <- bind_rows(LL1,LL2,LL3) %>% group_by(type,rc) %>% mutate(mv = mean(value))


  library(RColorBrewer)

  # Base palette (max 11 colors from RColorBrewer)
  base_spectral <- brewer.pal(11, "Spectral")
  # Extend to 20 colors using interpolation
  spectral_extended <- colorRampPalette(base_spectral)(length(table(LLD$mv)))

  Q <- IS(out = out,ind = ind,pixels = pixels)
  S0 <- ggplot(LLD)+
    geom_segment(aes(x = xs,xend = xs,
                     y=0,yend = value,
                     group = variable))+
    geom_point(aes(x=xs,
                   y=value,
                   group = variable,
                   fill = cut(mv,length(table(LLD$mv)))),
               size=2,pch=21,col=1)+
    facet_grid(~type, scales= "free_x")+
    ylab(ss)+
    scale_x_discrete("m/z" )+
    theme_bw()+
    scale_fill_manual(values = spectral_extended)+
    theme(legend.position = "none",text = element_text(size=10),
          axis.text.x = element_text(angle = 90, hjust = 1))
  #S0
  G <- Q + S0
  G <- G + plot_layout(widths = c(1, 3))
  return(G)
}
