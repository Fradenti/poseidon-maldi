#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.msPoser <- function(x,...){


  all_val = unlist(lapply(x,function(g) g$Elbo_val))
  all_len = unlist(lapply(x,function(g) length(g$Elbo_val)))

  e1 = x[[1]]$Elbo_val[!(x[[1]]$Elbo_val==0)]
  plot(e1,type="l",col="gray",
       ylim = c(min(all_val),max(all_val)),
       xlim=c(1,max(all_len)))
  for(j in 2:length(x)){
    e1 = x[[j]]$Elbo_val[!(x[[j]]$Elbo_val==0)]
    lines(e1,col="gray")
  }

}


