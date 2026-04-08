
#' Title
#'
#' @param MKCD
#'
#' @return
#' @export
#'
#' @examples
post_exp_val = function(res){
  MKCD <- res$MKCD
  A <-  cbind(post_mean = MKCD[,1],
        post_var = MKCD[,4]/(MKCD[,3]-1),
        avg_obs_assigned = colSums(res$RHO%*%apply(res$XI,c(2,3),sum)))
  print(round(A,4))
  return(A)
}

dmixnorm = function(x,pi,ucl,Res,k){
    sum(pi[,ucl[k]] * dnorm(x, Res[,1], sqrt(Res[,2])))
  }
#' Title
#'
#' @param res
#' @param k
#' @param xlim
#' @param n
#' @param brbr
#'
#' @return
#' @export
#'
#' @examples
draw_dens_pose <-function(res,
                          k=1,
                          xlim=NULL,
                          n_points_dens = 5000,
                          hist_breaks=25){

  if(is.null(xlim)) {xlim=c(min(res$Y),max(res$Y))}
  clu = estimate_clust_pose(res,bicl_mat = F)
  ucl = unique(clu$col_cl)
  Res = post_exp_val(res)
  L   = nrow(Res)
  if(res$model %in% c("fiSAN","fSAN")){
    pi = apply(res$B_star,2,function(t) t/sum(t))
  }else if(res$model=="CAM"){
    S  = res$a_bar_lk+res$b_bar_lk
    pi = res$a_bar_lk/S * apply(rbind(1,(res$b_bar_lk/S)[-L,]),2,cumprod)
  }

  hist(res$Y[,clu$col_cl==ucl[k]],freq = F,breaks=hist_breaks,xlim=xlim,
       fill="lightgray", border='darkblue',main = paste("Est. density for CC",ucl[k]))

  SEQ = seq(xlim[1],xlim[2],length.out=n_points_dens)
  yy = sapply(SEQ, function(x) dmixnorm(x = x,pi =
                                          pi,ucl = ucl,Res = Res,k=k))
  lines(yy~SEQ,lwd=1,col=4)
  points(pi[,ucl[k]]~Res[,1],type="h",lwd=2)
  norm <- (Res[,3]-min(Res[,3]))/(max(Res[,3])-min(Res[,3]))
  points(pi[,ucl[k]]~Res[,1],lwd=2, cex = norm*2)
}



#' Title
#'
#' @param res
#'
#' @return
#' @export
#'
#' @examples
estimate_clust_pose = function(res, bicl_mat = T){

  col_cl = apply(res$RHO,1,which.max)
  ucl = unique(col_cl)
  row_cl = matrix(NA,nrow(res$Y), length(ucl))

  for(k in seq_along(ucl)){
    row_cl[,k] = apply(res$XI[,ucl[k],],1,which.max)
  }
  colnames(row_cl) <-
    ifelse(ucl<=9, paste0("COL_CL#0",ucl), paste0("COL_CL#",ucl))
  ix = sort(colnames(row_cl),index=T)$ix

  ord_row_cl <- row_cl[,ix]

  if(bicl_mat){
    biclmat_mean = res$Y

    for(i in 1:ncol(res$Y)){
      id = which(ucl==col_cl[i])
      biclmat_mean[,i] = res$MKCD[as.matrix(row_cl)[,id],1]
    }
  }else{
    biclmat_mean = NULL
  }
  list(col_cl = col_cl,row_cl = ord_row_cl,bicl = biclmat_mean)
}



#' Title
#'
#' @param pixels
#' @param res
#'
#' @return
#' @export
#'
#' @examples
plot_col_clust_pose = function(pixels,res,ratio=1){
  cl <- estimate_clust_pose(res)
  coul <- colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(length(unique(cl$col_cl)))
  colnames(pixels)[1:2] = c("Var1","Var2")
  b <- as_tibble(pixels) %>% mutate(m = cl$col_cl)
  ggplot(b)+
    geom_tile(aes(x=Var1,
                  y=Var2,
                  fill=factor(as.numeric(factor(m)))),lwd=.05)+
    coord_fixed(ratio = ratio)+
    scale_fill_manual("cl",values = coul)+theme_dark()

}



#' Title
#'
#' @param res
#'
#' @return
#' @export
#'
#' @examples
filter_elbo <- function(res){

  n <- length(res$Elbo_val)
  inds_to_save <- which(1:n %% res$CEE == 0)

  cbind(ind = inds_to_save,res$Elbo_comp[inds_to_save,],
        elbo = res$Elbo_val[inds_to_save,],
        diff = c(0,diff(res$Elbo_val[inds_to_save,])))


  }



#' Title
#'
#' @param res_list
#'
#' @returns
#' @export
#'
#' @examples
extract_best_run_from_list <- function(res_list){

  el <- lapply(res_list,function(x) x$Elbo_val)
  len <- lapply(el,function(x) length(x))
  MAX <- lapply(el,function(x) max(x))

  len <- range(unlist(len))
  ran <- range(unlist(el))

  plot(el[[1]],
       ylim = c(ran[1],ran[2]),
       xlim = c(0,len[2]),
       type="l", col="lightgray")
  for(i in 2:length(el)){
    lines(el[[i]], col="lightgray")
    }

  ind <- which.max(MAX)
  lines(el[[ind]], col=4)
  return(res_list[[ind]])
}
