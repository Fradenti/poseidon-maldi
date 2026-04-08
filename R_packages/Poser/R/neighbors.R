#' Title
#'
#' @param NeighMat
#'
#' @return
#' @export
#'
#' @examples
create_cpp_neigh_list <- function(NeighMat){
  rast = raster::raster(NeighMat)
  ADJ = raster::adjacent(rast,
                 cells = 1:raster::ncell(rast),
                 directions=8,
                 sorted=T)

  nnList = list()
  for(i in 1:ncell(rast)){
    nnList[[i]] = ADJ[which(ADJ[,1]==i),2]
  }
  nnList_cpp = lapply(nnList,function(x)x-1)
  return(nnList_cpp)
}
