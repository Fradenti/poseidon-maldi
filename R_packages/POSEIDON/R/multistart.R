

#' Extract best run from a list of Poseidons
#'
#' @param VVm
#' @param min
#'
#' @return
#' @export
#'
extract_best_poseidon = function(VVm,min = NULL, plot = FALSE){

    MIN = min(unlist(lapply(VVm, function(x) min(x$Elbo_val))))
    MAX = max(unlist(lapply(VVm, function(x) max(x$Elbo_val))))
    MAXLEN = max(unlist(lapply(VVm, function(x) length(x$Elbo_val))))
    ind = which.max(unlist(lapply(VVm, function(x) max(x$Elbo_val))))
    if(!is.null(min)){
      MIN = min
    }
    if(plot){
      plot(VVm[[1]]$Elbo_val,ylim=c(MIN,MAX),xlim=c(1,MAXLEN),type="l")
      for(g in 1:length(VVm)){
        lines(VVm[[g]]$Elbo_val)
      }
      VV = VVm[[ind]]
      points(VV$Elbo_val,col=4,type="b",cex=.5)
    }else{
      VV = VVm[[ind]]
    }
  return(VV)
}


#' Title
#'
#' @param VVm
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
elbo_diff_poseidon = function(VVm, ylim=c(-.4,.1)){
  MAXLEN = max(unlist(lapply(VVm, function(x) length(x$Elbo_val))))
  plot(diff(VVm[[1]]$Elbo_val),
       ylim=ylim,
       xlim=c(1,MAXLEN),
       type="l")
  for(g in 1:length(VVm)){

    lines(diff(VVm[[g]]$Elbo_val))

  }
  abline(h=0,col=2)
}


#' Title
#'
#' @param VVm
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
elbo_ms_poseidon = function(VVm, ylim=NULL){
  MAXLEN = max(unlist(lapply(VVm, function(x) length(x$Elbo_val))))
  plot((VVm[[1]]$Elbo_val),
       ylim=ylim,
       xlim=c(1,MAXLEN),
       type="l")
  for(g in 1:length(VVm)){

    lines(VVm[[g]]$Elbo_val)

  }
  abline(h=0,col=2)
}
