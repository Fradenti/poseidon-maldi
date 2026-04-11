

estimate_clust_pose = function(VV, bicl_mat = T){

  col_cl = apply(VV$RHO,1,which.max)
  ucl = unique(col_cl)
  row_cl = matrix(NA,nrow(VV$Y), length(ucl))

  for(k in seq_along(ucl)){
    row_cl[,k] = apply(VV$XI[,ucl[k],],1,which.max)
  }
  colnames(row_cl) <-
    ifelse(ucl<=9, paste0("COL_CL#0",ucl), paste0("COL_CL#",ucl))
  ix = sort(colnames(row_cl),index=T)$ix

  ord_row_cl <- row_cl[,ix]

  if(bicl_mat){
    biclmat_mean = VV$Y

    for(i in 1:ncol(VV$Y)){
      id = which(ucl==col_cl[i])
      biclmat_mean[,i] = VV$MKCD[as.matrix(row_cl)[,id],1]
    }
  }else{
    biclmat_mean = NULL
  }
  list(col_cl = col_cl,row_cl = ord_row_cl,bicl = biclmat_mean)
}


estimate_clust_poseidon <- function(q, bicl_mat = T){

  col_cl = (apply(q$RHO,1,which.max))
  ucl = (unique(col_cl))

  row_cl = list()

  for(x in 1:nrow(q$Y_list)){
    subrow = matrix(NA,
                    nrow(q$Y_list[[x]]),
                    length(seq_along(ucl)))
    for(k in seq_along(ucl)){
      subrow[,k] = apply(q$XI[[x]][,ucl[k],],
                         1,
                         which.max)
      colnames(subrow) = paste0("COL_CL#",ucl)
    }
    row_cl[[x]] <- subrow
  }


  if(bicl_mat){
    biclmat_all <- list()
    for(x in 1:nrow(q$Y_list)){
      biclmat_mean <- q$Y[[x]]

      for(i in 1:ncol(q$Y[[x]])){
        id = which(ucl==col_cl[i])
        biclmat_mean[,i] <- q$MKCD[row_cl[[x]][,id],1,x]
      }
      biclmat_all[[x]] <-  biclmat_mean
    }
  }else{
    biclmat_all = NULL
  }
  list(col_cl = col_cl, row_cl = row_cl, bicl = biclmat_all)
}


mapper2 <- function(Y) {
  m = min(Y)-1e-3
  M = max(Y)+1e-3
  z = (Y-m)/(M-m)
  qnorm(z) - mean(qnorm(z))
}

extract_median_cl <- function(colcl, data){
  ucl <- unique(colcl)
  M   <- cbind(ucl,ucl)
  for(j in seq_along(ucl)){
    M[j,2] <- mean(apply(as.matrix(data[,colcl == ucl[j]]),2,median))
  }
  MM <- as_tibble(M) %>% arrange(-ucl)
  colcl2 <- colcl
  for(j in seq_along(ucl)){
    colcl2[colcl==ucl[j]] <- pull(MM[j,2])
  }
  return(colcl2)
}


mapping_cl_pixels = function(ind,clord){
  num_clord <- as.numeric(factor(clord))
  Q <- ggplot(D_hpm)+
    geom_tile(aes(x=pixels[,1],
                  y=pixels[,2]),fill="lightgray",col=1,alpha=.6) +
    geom_tile(data = D_hpm[num_clord %in% ind,],
              aes(x=pixels[num_clord %in% ind,1],
                  y=pixels[num_clord %in% ind,2],
                  fill=factor(v)),col=1) +
    theme_void() + scale_fill_brewer(palette = "Set1")+
    theme(legend.position = "bottom") + labs(fill = "cluster")
  return(Q)
}

mapping_cl_mz = function(ind,clord, ylim = NULL, name  = NULL){

  num_clord <- as.numeric(factor(clord))
  D  <- as_tibble(t(r[,num_clord %in% ind])) %>%
    mutate(cl = num_clord[num_clord %in% ind])

  mD <- reshape2::melt(D,"cl") %>% mutate(var = as.numeric(factor(variable)))
  if(!is.null(name)){
    mD <- mD %>% mutate(var = factor(name[var],levels = name))%>% ungroup()
  }

  msd <- mD %>% group_by(var,cl) %>% summarise(me = median(value)) %>% ungroup()
  s <- ggplot()+
    geom_boxplot(data=mD,
                 aes(x=factor(var),
                     y=value, col =factor(cl) ),alpha=.1) +
    geom_line(data = msd,
              mapping =
                aes(x=as.numeric(var),
                    y=me,
                    col =factor(cl) ),lwd=1) +
    theme_bw()+scale_color_brewer(palette = "Set1")+
    theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  if(is.null(ylim)){
    return(s)
  }else{
    s <- s+coord_cartesian(ylim = ylim)
    return(s)
  }
}

# poseidon ----------------------------------------------------------------


# extract_median_cl <- function(colcl){
#   ucl <- unique(colcl)
#   M   <- cbind(ucl,ucl)
#   for(j in seq_along(ucl)){
#     M[j,2] <- median(c(r[,colcl == ucl[j]]))
#   }
#   MM <- as_tibble(M) #%>% arrange(-V2)
#   colcl2 <- colcl
#   for(j in seq_along(ucl)){
#     colcl2[colcl==ucl[j]] <- pull(MM[j,2])
#   }
#   return(colcl2)
# }

estimate_clust_posidon <- function(q, bicl_mat = T){

  col_cl = (apply(q$RHO,1,which.max))
  ucl = (unique(col_cl))



  row_cl = list()

  for(x in 1:nrow(q$Y_list)){
    subrow = matrix(NA,nrow(q$Y_list[[x]]),length(seq_along(ucl)))
    for(k in seq_along(ucl)){
      subrow[,k] = apply(q$XI[[x]][,ucl[k],],1,which.max)
      colnames(subrow) = paste0("COL_CL#",ucl)
    }
    row_cl[[x]] <- subrow
  }


  if(bicl_mat){
    biclmat_all <- list()
    for(x in 1:nrow(q$Y_list)){
      biclmat_mean <- q$Y[[x]]

      for(i in 1:ncol(q$Y[[x]])){
        id = which(ucl==col_cl[i])
        biclmat_mean[,i] <- q$MKCD[row_cl[[x]][,id],1,x]
      }
      biclmat_all[[x]] <-  biclmat_mean
    }
  }else{
    biclmat_all = NULL
  }
  list(col_cl = col_cl,row_cl = row_cl,bicl = biclmat_all)
}


# bicl2 <- function(clustering, ind, all,
#                   palette_CC="davos",
#                   palette_RC="imola"){
#
#   if( !all(ind %in% ucl)){
#     stop(paste("wrong clustering index: use one of",
#                paste(paste("-", sort(ucl)), collapse = '') ) )
#   }
#   i = sort(clustering$col_cl,ind=T)
#
#   S1   <- all[[1]][,i$x==ucl[ind]]
#   S2   <- all[[2]][,i$x==ucl[ind]]
#   S3   <- all[[3]][,i$x==ucl[ind]]
#
#
#   mS1  <- reshape2::melt(S1) %>% mutate(data = "Lipids")
#   mS2  <- reshape2::melt(S2) %>% mutate(data = "Glicans")
#   mS3  <- reshape2::melt(S3) %>% mutate(data = "Peptides")
#
#   pp <- as_tibble(pixels) %>% mutate(colcl = clustering$col_cl)
#
#   aa <- ggplot()+
#     geom_tile(data = pp,aes(x=pixels[,1],
#                             y=pixels[,2]),alpha=.4,col=1,fill="lightgray")+
#     geom_tile(aes(x=pixels[pp$colcl==ucl[ind],1],
#                   y=pixels[pp$colcl==ucl[ind],2]),alpha=.9,col=1)+
#     scico::scale_fill_scico(palette = palette_CC)+
#     facet_wrap(~"segmentation")+
#     theme_bw()
#
#   bb1 <- ggplot(mS1)+
#     geom_tile(aes(x=Var2,y=-Var1,fill=value))+
#     scico::scale_fill_scico(palette = palette_RC,
#                             limits= c(min(unlist(all[[1]])),max(unlist(all[[1]]))))+
#     theme_bw()+facet_wrap(~data)+
#     theme(legend.position = "none")
#   bb2 <- ggplot(mS2)+
#     geom_tile(aes(x=Var2,y=-Var1,fill=value))+
#     scico::scale_fill_scico(palette = palette_RC,
#                             limits= c(min(unlist(all[[2]])),max(unlist(all[[2]]))))+
#     theme_bw()+facet_wrap(~data)+
#     theme(legend.position = "none")
#   bb3 <- ggplot(mS3)+
#     geom_tile(aes(x=Var2,y=-Var1,fill=value))+
#     scico::scale_fill_scico(palette = palette_RC,
#                             limits= c(min(unlist(all[[3]])),max(unlist(all[[3]]))))+
#     theme_bw()+facet_wrap(~data)+
#     theme(legend.position = "none")
#
#
#   aa|bb1|bb2|bb3
# }
#

# bicl1 <- function(data, clustering, ind, all){
#   i = sort(clustering$col_cl,ind=T)
#
#   o   <- apply(clustering$row_cl[[data]],2,function(x) sort(x,index=T))
#   ucl <- unique(clustering$col_cl)
#   S   <- all[[data]][o[[ind]]$ix,i$x==ucl[ind]]
#   mS  <- reshape2::melt(S)
#
#   ggplot(mS)+
#     geom_tile(aes(x=Var2,y=-Var1,fill=value))+
#     scale_fill_viridis_c(option = "H")+
#     geom_hline(yintercept = -c(cumsum(table(o[[ind]]$x)))-.5)+
#     theme_bw()
# }

# ex bicl3 pose
Pose_cl_boxpl <- function(clustering, ind, molec = c("gly","lip","pep"),
                          pal_CC = "davos",
                          pal_RC = "oslo"){

  ucl <- unique(clustering$col_cl)
  if( !all(ind %in% ucl)){
    stop(paste("wrong clustering index: use one of",
               paste(paste("-", sort(ucl)), collapse = '') ) )
  }
  i = sort(clustering$col_cl,ind=T)

  S1   <- all[[1]][,i$x==ind]
  S2   <- all[[2]][,i$x==ind]
  S3   <- all[[3]][,i$x==ind]

  mS1  <- reshape2::melt(S1) %>% mutate(data = "Glicans")
  mS2  <- reshape2::melt(S2) %>% mutate(data = "Lipids")
  mS3  <- reshape2::melt(S3) %>% mutate(data = "Peptides")

  ucl <- unique(clustering$col_cl)
  if( !all(ind %in% ucl)){
    stop(paste("wrong clustering index: use one of",
               paste(paste("-", sort(ucl)), collapse = '') ) )
  }
  i = sort(clustering$col_cl,ind=T)

  S1   <- all[[1]][,i$x==ind]
  mS1  <- reshape2::melt(S1) %>% mutate(data = "Glycans")
  S2   <- all[[2]][,i$x==ind]
  mS2  <- reshape2::melt(S2) %>% mutate(data = "Lipids")
  S3   <- all[[3]][,i$x==ind]
  mS3  <- reshape2::melt(S3) %>% mutate(data = "Peptides")



  Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = (clustering$col_cl))
  Grd <- expand_grid(x=0:(max(pixels[,1])+1),
                     y=0:(max(pixels[,2])+1)) %>%
    mutate(p = paste0(x,"-",y)) %>% left_join(Px,by = "p")

  Grd <- Grd %>% mutate(ccl2 = ifelse(is.na(ccl),0,ccl))

  aa <- ggplot()+
    geom_tile(data=Grd,aes(x = X,
                           y = Y),alpha=.4,
              col=1,fill="lightgray")+
    geom_tile(data = Grd %>% filter(ccl2%in%ind),
              aes(x = X,
                  y = Y, fill=factor(ccl)),
              alpha=.9, col=1) +
    geom_tileborder(data = Grd,
                    aes(x = x,
                        y = y,
                        group = 1,
                        grp = ifelse( ccl2 == ind, 1, 0)), lwd=.5)+
    scale_fill_simpsons()+
    facet_wrap(~ paste("Poseidon - Column Cluster",ind))+
    theme_bw()+
    theme(legend.position = "none",text = element_text(size=16)) +
    xlab("x coord.") + ylab("y coord.")


  q <-   ifelse(ind<=9, paste0("COL_CL#0",ind), paste0("COL_CL#",ind))

  if(molec == "gly"){
    ind1 <- which(paste0("COL_CL#",(sort(ucl))) %in% paste0("COL_CL#",ind))
    rcl1 <- clustering$row_cl[,ind1]
    mS1 <- mS1 %>% mutate(Var3 = rcl1[Var1], x= allx__g[Var2])
    bb <- ggplot(mS1)+
      geom_boxplot(aes(x=x,y=value,group=factor(Var1),
                       fill=factor(Var3)),lwd=.4,
                   outlier.size = .25)+
      theme_bw()+facet_wrap(~data) +ylim(c(-2,6))+
      scico::scale_fill_scico_d(palette = pal_RC,begin=.1)+
      scale_x_discrete("m/z",breaks = limx__g )+
      theme(legend.position = "none",text = element_text(size=16)) +
      ylab("Abundance")


  }else if(molec == "lip"){
    ind2 <- which(paste0("COL_CL#",(sort(ucl))) %in% paste0("COL_CL#",ind))
    rcl2 <- clustering$row_cl[,ind2]
    mS2 <- mS2 %>% mutate(Var3 = rcl2[Var1], x= allx__g[Var2])
    bb <- ggplot(mS2)+
      geom_boxplot(aes(x=(x),y=value,group=Var1,
                       fill=factor(Var3)),lwd=.4,
                   outlier.size = .25)+
      theme_bw()+facet_wrap(~data) +ylim(c(-2,6))+
      scico::scale_fill_scico_d(palette = pal_RC,end = .8)+
      scale_x_discrete("m/z",breaks = limx__g )+
      theme(legend.position = "none",text = element_text(size=16)) +
      ylab("Abundance")
  }else if(molec == "pep"){
    ind3 <- which(paste0("COL_CL#",(sort(ucl))) %in% paste0("COL_CL#",ind))

    rcl3 <- clustering$row_cl[,ind3]
    mS3 <- mS3 %>% mutate(Var3 = rcl3[Var1], x= allx__p[Var2])
    bb <- ggplot(mS3)+
      geom_boxplot(aes(x=x,y=value,group=Var1,
                       fill=factor(Var3)),lwd=.4,
                   outlier.size = .25)+
      theme_bw()+facet_wrap(~data) +ylim(c(-2,6))+
      scico::scale_fill_scico_d(palette = pal_RC,end = .8)+
      scale_x_discrete("m/z",
                       breaks = limx__p )+
      theme(legend.position = "none",
            text = element_text(size=16)) +
      ylab("Abundance")


  }

  (aa|bb)
}


# ex bicl3
Poseidon_cl_boxpl <- function(clustering, ind,
                              pal_CC = "davos",
                              pal_RC = "oslo"){

  ucl <- unique(clustering$col_cl)
  if( !all(ind %in% ucl)){
    stop(paste("wrong clustering index: use one of",
               paste(paste("-", sort(ucl)), collapse = '') ) )
  }
  i = sort(clustering$col_cl,ind=T)

  S1   <- all[[1]][,i$x==ind]
  S2   <- all[[2]][,i$x==ind]
  S3   <- all[[3]][,i$x==ind]

  mS1  <- reshape2::melt(S1) %>% mutate(data = "Glycans")
  mS2  <- reshape2::melt(S2) %>% mutate(data = "Lipids")
  mS3  <- reshape2::melt(S3) %>% mutate(data = "Peptides")


  Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = (clustering$col_cl))
  Grd <- expand_grid(x=0:(max(pixels[,1])+1),
                     y=0:(max(pixels[,2])+1)) %>%
    mutate(p = paste0(x,"-",y)) %>% left_join(Px,by = "p")

  Grd <- Grd %>% mutate(ccl2 = ifelse(is.na(ccl),0,ccl))

  aa <- ggplot()+
    geom_tile(data=Grd,aes(x = X,
                           y = Y),alpha=.4,
              col=1,fill="lightgray")+
    geom_tile(data = Grd %>% filter(ccl2%in%ind),
              aes(x = X,
                  y = Y, fill=factor(ccl)),
              alpha=.9, col=1) +
    geom_tileborder(data = Grd,
                    aes(x = x,
                        y = y,
                        group = 1,
                        grp = ifelse( ccl2 == ind, 1, 0)), lwd=.5)+
    scale_fill_simpsons()+
    facet_wrap(~ paste("Poseidon - Column Cluster",ind))+
    theme_bw()+
    theme(legend.position = "none",text = element_text(size=16)) +
    xlab("x coord.") + ylab("y coord.")

  ind1 <- which(colnames(clustering$row_cl[[1]])==paste0("COL_CL#",ind[1]))
  rcl1 <- clustering$row_cl[[1]][,ind1]
  mS1 <- mS1 %>% mutate(Var3 = rcl1[Var1], x= allx__g[Var2])
  bb1 <- ggplot(mS1)+
    geom_boxplot(aes(x=x,y=value,group=factor(Var1),
                     fill=factor(Var3)),lwd=.4,
                 outlier.size = .25)+
    theme_bw()+facet_wrap(~data) +ylim(c(-2,6))+
    scico::scale_fill_scico_d(palette = pal_RC,begin=.1)+
    scale_x_discrete("m/z",breaks = limx__g )+
    theme(legend.position = "none",text = element_text(size=16)) +
    ylab("Abundance")

  ind2 <- which(colnames(clustering$row_cl[[2]])==paste0("COL_CL#",ind[1]))
  rcl2 <- clustering$row_cl[[2]][,ind1]
  mS2 <- mS2 %>% mutate(Var3 = rcl2[Var1], x= allx__g[Var2])
  bb2 <- ggplot(mS2)+
    geom_boxplot(aes(x=(x),y=value,group=Var1,
                     fill=factor(Var3)),lwd=.4,
                 outlier.size = .25)+
    theme_bw()+facet_wrap(~data) +ylim(c(-2,6))+
    scico::scale_fill_scico_d(palette = pal_RC,end = .8)+
    scale_x_discrete("m/z",breaks = limx__g )+
    theme(legend.position = "none",text = element_text(size=16)) +
    ylab("Abundance")

  ind3 <- which(colnames(clustering$row_cl[[3]])==paste0("COL_CL#",ind[1]))
  rcl3 <- clustering$row_cl[[3]][,ind1]
  mS3 <- mS3 %>% mutate(Var3 = rcl3[Var1], x= allx__p[Var2])
  bb3 <- ggplot(mS3)+
    geom_boxplot(aes(x=x,y=value,group=Var1,
                     fill=factor(Var3)),lwd=.4,
                 outlier.size = .25)+
    theme_bw()+facet_wrap(~data) +ylim(c(-2,6))+
    scico::scale_fill_scico_d(palette = pal_RC,end = .8)+
    scale_x_discrete("m/z",
                     breaks = limx__p )+
    theme(legend.position = "none",
          text = element_text(size=16)) +
    ylab("Abundance")


  (aa|bb1)/(bb2|bb3)
}

# ex bicl0
Multiple_CC <- function(clustering, ind, model = "Poseidon"){

  ucl <- unique(clustering$col_cl)
  if( !all(ind %in% ucl)){
    stop(paste("wrong clustering index: use one of",
               paste(paste("-", sort(ucl)), collapse = '') ) )
  }
  i = sort(clustering$col_cl,ind=T)

  Px <- pixels %>% mutate(p = paste0(X,"-",Y), ccl = (clustering$col_cl))
  Grd <- expand_grid(x=0:(max(pixels[,1])+1),
                     y=0:(max(pixels[,2])+1)) %>%
    mutate(p = paste0(x,"-",y)) %>% left_join(Px,by = "p")

  Grd <- Grd %>% mutate(ccl2 = ifelse(is.na(ccl),0,ccl))
  #View(Grd)
  #pp <- as_tibble(pixels) %>% mutate(colcl = clustering$col_cl)

  aa <- ggplot()+
    geom_tile(data=Grd,aes(x = X,
                           y = Y),alpha=.4,
              col=1,fill="lightgray")+
    geom_tile(data = Grd %>% filter(ccl2%in%ind),
              aes(x = X,
                  y = Y, fill=factor(ccl2)),
              alpha=.9, col=1) +
    geom_tileborder(data = Grd,
                    aes(x = x,
                        y = y,
                        group = 1,
                        grp = ifelse( ccl2%in%ind, 1, 0)), lwd=1)+
    scale_fill_brewer(palette = "Spectral") +
    facet_wrap(~ paste0(model," - Multiple Column Clusters"))+
    #    facet_wrap(~"Poseidon - Image Segmentation") +
    theme(legend.position = "none",text = element_text(size=16)) +
    xlab("x coord.") + ylab("y coord.")

  aa
}

# ex bicl4
Medians_means_RC <- function(clustering, inds,
                             all, median = TRUE){
  n <- length(inds)
  ucl <- unique(clustering$col_cl)
  if( !all(inds %in% ucl)){
    stop(paste("wrong clustering index: use one of",
               paste(paste("-", sort(ucl)), collapse = '') ) )
  }

  L1 = c()
  L2 = c()
  L3 = c()

  if(median){
    for(i in seq_along(inds)){

      L1 = cbind(L1,apply(all[[1]][,clustering$col_cl==inds[i]],1,median))
      L2 = cbind(L2,apply(all[[2]][,clustering$col_cl==inds[i]],1,median))
      L3 = cbind(L3,apply(all[[3]][,clustering$col_cl==inds[i]],1,median))

      ss <-   ylab("Median Abundance")

    }
  }else{
    for(i in seq_along(inds)){

      L1 = cbind(L1,rowMeans(all[[1]][,clustering$col_cl==inds[i]]))
      L2 = cbind(L2,rowMeans(all[[2]][,clustering$col_cl==inds[i]]))
      L3 = cbind(L3,rowMeans(all[[3]][,clustering$col_cl==inds[i]]))

      ss <-   ylab("Mean Abundance")
    }

  }

  colnames(L1) = inds
  colnames(L2) = inds
  colnames(L3) = inds

  ind_col <- which(colnames(clustering$row_cl[[1]]) %in% paste0("COL_CL#",inds))
  glis <- c(clustering$row_cl[[1]][,ind_col])
  ind_col <- which(colnames(clustering$row_cl[[2]]) %in% paste0("COL_CL#",inds))
  lips <- c(clustering$row_cl[[2]][,ind_col])
  ind_col <- which(colnames(clustering$row_cl[[3]]) %in% paste0("COL_CL#",inds))
  peps <- c(clustering$row_cl[[3]][,ind_col])


  LL1 <- as.data.frame(L1) %>% mutate(xs = allx__g) %>% reshape2::melt() %>% mutate(type = "Glycans",  rc = (glis))
  LL2 <- as.data.frame(L2) %>% mutate(xs = allx__l) %>% reshape2::melt() %>% mutate(type = "Lipids",   rc = (lips))
  LL3 <- as.data.frame(L3) %>% mutate(xs = allx__p) %>% reshape2::melt() %>% mutate(type = "Peptides", rc = (peps))


  LLD <- bind_rows(LL1,LL2,LL3)


  p1 <- ggplot(LL1)+
    geom_line(aes(x=xs,y=value,group = variable,col=variable))+
    geom_point(aes(x=xs,y=value,group = variable,col=variable,
                   shape = factor(rc)),
               size=2)+
    facet_wrap(~type, scales = "free_x")+
    scale_shape_manual(values = 0:length(unique(factor(glis))),
                       guide= "none")+
    ggsci::scale_color_simpsons(name="Column Cluster")+
    scale_x_discrete("m/z",breaks = limx__g )+
    theme_bw()+
    ss
  p1

  p2 <- ggplot(LL2)+
    geom_line(aes(x=xs,y=value,group = variable,col=variable))+
    geom_point(aes(x=xs,y=value,group = variable,col=variable,
                   shape = factor(rc)),size=2)+
    facet_wrap(~type, scales = "free_x")+
    scale_shape_manual(values = 0:length(unique(factor(lips))),
                       guide= "none")+
    ggsci::scale_color_simpsons(name="Column Cluster")+
    scale_x_discrete("m/z",breaks = limx__l )+
    theme_bw()+
    ss
  p3 <- ggplot(LL3)+
    geom_line(aes(x=xs,y=value,group = variable,col=variable))+
    geom_point(aes(x=xs,y=value,group = variable,col=variable,
                   shape = factor(rc)),size=2)+
    facet_wrap(~type, scales = "free_x")+
    #scale_color_brewer("Column Cluster", palette = "Set1")+
    scale_shape_manual(values = 0:length(unique(factor(peps))),
                       guide= "none")+
    ggsci::scale_color_simpsons(name="Column Cluster")+
    scale_x_discrete("m/z",breaks = limx__p )+
    theme_bw()+
    ss


  (p1/p2/p3)+
    plot_layout(guides = "collect") & theme(legend.position = 'bottom',
                                            text=element_text(size=16))
  # par(mfrow=c(3,1))
  # matplot(L1,type="b",lwd=2,lty=1,ylab = "Glicani")
  # matplot(L2,type="b",lwd=2,lty=1,ylab = "Lipidi")
  # matplot(L3,type="b",lwd=2,lty=1,ylab = "Peptidi")
  # par(mfrow=c(1,1))

}

#
# bicl4_pose <- function(clustering, inds, data, median = TRUE){
#   n <- length(inds)
#
#   L1 = c()
#
#   if(median){
#     for(i in seq_along(inds)){
#
#       L1 = cbind(L1,apply(data[,clustering$col_cl==inds[i]],1,median))
#
#     }
#   }else{
#     for(i in seq_along(inds)){
#
#       L1 = cbind(L1,rowMeans(data[,clustering$col_cl==inds[i]]))
#
#     }
#
#   }
#
#   matplot(L1,type="b",lwd=2,lty=1,ylab = "Lip")
#
# }
#



Pose_Medians_means_RC <- function(clustering, inds,
                                  all, median = TRUE, molec = c("gly","lip","pep")){
  n <- length(inds)
  ucl <- unique(clustering$col_cl)
  if( !all(inds %in% ucl)){
    stop(paste("wrong clustering index: use one of",
               paste(paste("-", sort(ucl)), collapse = '') ) )
  }

  L1 = c()
  L2 = c()
  L3 = c()

  if(median){
    for(i in seq_along(inds)){

      L1 = cbind(L1,apply(all[[1]][,clustering$col_cl==inds[i]],1,median))
      L2 = cbind(L2,apply(all[[2]][,clustering$col_cl==inds[i]],1,median))
      L3 = cbind(L3,apply(all[[3]][,clustering$col_cl==inds[i]],1,median))

      ss <-   ylab("Median Abundance")

    }
  }else{
    for(i in seq_along(inds)){

      L1 = cbind(L1,rowMeans(all[[1]][,clustering$col_cl==inds[i]]))
      L2 = cbind(L2,rowMeans(all[[2]][,clustering$col_cl==inds[i]]))
      L3 = cbind(L3,rowMeans(all[[3]][,clustering$col_cl==inds[i]]))

      ss <-   ylab("Mean Abundance")
    }

  }

  colnames(L1) = inds
  colnames(L2) = inds
  colnames(L3) = inds






  if(molec == "gly"){


    ind_col <- which(paste0("COL_CL#",(sort(ucl))) %in% paste0("COL_CL#",inds))
    glis <- c(clustering$row_cl[,ind_col])
    LL1 <- as.data.frame(L1) %>% mutate(xs = allx__g) %>% reshape2::melt() %>% mutate(type = "Glycans",  rc = (glis))

    pp2 <- ggplot(LL1)+
      geom_line(aes(x=xs,y=value,group = variable,col=variable))+
      geom_point(aes(x=xs,y=value,group = variable,col=variable,
                     shape = factor(rc)),
                 size=2)+
      facet_wrap(~type, scales = "free_x")+
      scale_shape_manual(values = 0:length(unique(factor(glis))),
                         guide= "none")+
      ggsci::scale_color_simpsons(name="Column Cluster")+
      scale_x_discrete("m/z",breaks = limx__g )+
      theme_bw()+
      ss

  }else if(molec== 'lip'){
    ind_col <- which(paste0("COL_CL#",(sort(ucl))) %in% paste0("COL_CL#",inds))
    lips <- c(clustering$row_cl[,ind_col])
    LL2 <- as.data.frame(L2) %>% mutate(xs = allx__l) %>% reshape2::melt() %>% mutate(type = "Lipids",   rc = (lips))

    pp2 <- ggplot(LL2)+
      geom_line(aes(x=xs,y=value,group = variable,col=variable))+
      geom_point(aes(x=xs,y=value,group = variable,col=variable,
                     shape = factor(rc)),size=2)+
      facet_wrap(~type, scales = "free_x")+
      scale_shape_manual(values = 0:length(unique(factor(lips))),
                         guide= "none")+
      ggsci::scale_color_simpsons(name="Column Cluster")+
      scale_x_discrete("m/z",breaks = limx__l )+
      theme_bw()+
      ss
  }else if(molec == "pep"){
    ind_col <- which(paste0("COL_CL#",(sort(ucl))) %in% paste0("COL_CL#",inds))
    peps <- c(clustering$row_cl[,ind_col])
    LL3 <- as.data.frame(L3) %>% mutate(xs = allx__p) %>% reshape2::melt() %>% mutate(type = "Peptides", rc = (peps))

    pp2 <- ggplot(LL3)+
      geom_line(aes(x=xs,y=value,group = variable,col=variable))+
      geom_point(aes(x=xs,y=value,group = variable,col=variable,
                     shape = factor(rc)),size=2)+
      facet_wrap(~type, scales = "free_x")+
      #scale_color_brewer("Column Cluster", palette = "Set1")+
      scale_shape_manual(values = 0:length(unique(factor(peps))),
                         guide= "none")+
      ggsci::scale_color_simpsons(name="Column Cluster")+
      scale_x_discrete("m/z",breaks = limx__p )+
      theme_bw()+
      ss
  }

  (pp2)+
    plot_layout(guides = "collect") & theme(legend.position = 'bottom',
                                            text=element_text(size=16))
  # par(mfrow=c(3,1))
  # matplot(L1,type="b",lwd=2,lty=1,ylab = "Glicani")
  # matplot(L2,type="b",lwd=2,lty=1,ylab = "Lipidi")
  # matplot(L3,type="b",lwd=2,lty=1,ylab = "Peptidi")
  # par(mfrow=c(1,1))

}
