#' Title
#'
#' @param x  ...
#'
#' @return ...
#' @export
#'
#' @examples
#' 2+2
norm <- function(x){
  (x-min(x))/((max(x)-min(x)))
  }


#' Title
#'
#' @param MKCD
#'
#' @return
#' @export
#'
#' @examples
post_exp_val = function(MKCD){
  cbind(MKCD[,1],MKCD[,4]/(MKCD[,3]-1))
}

#' Title
#'
#' @param VV
#' @param k
#' @param xlim
#' @param n
#' @param brbr
#'
#' @return
#' @export
#'
#' @examples
draw_dens_poseidon = function(VV,
                              k=1,
                              xlim=NULL,
                              n_points_dens = 5000,
                              hist_breaks=25){

  if(is.null(xlim)) {
    xlim = list()
    for(i in 1:length(VV$Y_list)){
      xlim[[i]]=c(min(VV$Y_list[[i]]),max(VV$Y_list[[i]]))
    }
  }

  clu = estimate_clust_poseidon(VV,bicl_mat = F)
  ucl = unique(clu$col_cl)

  par(mfrow=c(2,1))
  L   = nrow(VV$MKCD[,,1])

  for(i in 1:length(VV$Y_list)){


    Res = post_exp_val(MKCD = VV$MKCD[,,i])
    pi = apply(VV$B_star[,,i], 2, function(t) t/sum(t))
    hist(VV$Y_list[[i]][,clu$col_cl==ucl[k]],freq = F,breaks=hist_breaks,xlim=xlim[[i]])
    ff = function(x){
      sum(pi[,ucl[k]]*dnorm(x,Res[,1],sqrt(Res[,2])))
    }

    SEQ = seq(xlim[[i]][1],xlim[[i]][2],length.out=n_points_dens)
    yy = sapply(SEQ, function(x) ff(x))
    lines(yy~SEQ,lwd=1,col=4)
    points(norm(pi[,ucl[k]])*max(yy)/2~Res[,1],type="h",lwd=2)
    points(norm(pi[,ucl[k]])*max(yy)/2~Res[,1],lwd=2)
  }
}



#' Title
#'
#' @param VV
#'
#' @return
#' @export
#'
#' @examples
estimate_clust_poseidon <- function(q, bicl_mat = T){

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
  list(col_cl = col_cl, row_cl = row_cl, bicl = biclmat_all)
}





#' Title
#'
#' @param VV
#'
#' @return
#' @export
#'
#' @examples
filter_elbo <- function(VV){

  n <- length(VV$Elbo_val)
  inds_to_save <- which(1:n %% VV$CEE == 0)

  cbind(ind = inds_to_save,VV$Elbo_comp[inds_to_save,],
        elbo = VV$Elbo_val[inds_to_save,],
        diff = c(0,diff(VV$Elbo_val[inds_to_save,])))


  }

