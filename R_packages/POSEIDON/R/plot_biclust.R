#' Title
#'
#' @param q ...
#' @param data ...
#'
#' @return ...
#' @export
#'
plot_bicluster_poseidon <- function(poseidon_fit, data=1){

  Y   <- poseidon_fit$Y_list[[data]]
  # pheatmap::pheatmap(Y,cluster_rows = F,cluster_cols = F)

  cls <- estimate_clust_poseidon(poseidon_fit)
  a1  <- sort(cls$col_cl,index=T)
  ucl <- unique(cls$col_cl)

  #indcol <- which(diff(as.numeric(as.factor(a1$x)))>0)
  b1     <- apply(cls$row_cl[[data]],2,function(r) sort(r,index=T))



  indrows <- lapply(b1,function(x) which(diff(as.numeric(as.factor(x$x)))>0))

  tempY <- c()
  incol2 <- c()
  for(i in 1:length(b1)){
    tempY <- cbind( tempY, Y[b1[[i]]$ix,a1$ix[a1$x==ucl[i]]] )
    check <- length(ncol(Y[b1[[i]]$ix,a1$ix[a1$x==ucl[i]]]))
    incol2[i] = ifelse(check==0,0,ncol(Y[b1[[i]]$ix,a1$ix[a1$x==ucl[i]]]))
  }
  #pheatmap::pheatmap(Y[b1[[i]]$ix,a1$ix[a1$x==ucl[i]]])

  mY <- reshape2::melt(tempY)
  table(cls$col_cl)

  ggplot(mY)+theme_bw()+
    geom_tile(aes(x=Var2, y=rev(Var1), fill=value))+
    geom_vline(xintercept = cumsum(incol2) + .5, lwd=1)

}
