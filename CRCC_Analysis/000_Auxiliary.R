mapper2 <- function(Y) {
  m = min(Y)-1e-3
  M = max(Y)+1e-3
  z = (Y-m)/(M-m)
  qnorm(z) - mean(qnorm(z))
}

dmixnorm = function(x,pi,ucl,Res,k){
  sum(pi[,ucl[k]] * dnorm(x, Res[,1], sqrt(Res[,2])))
}

# SINGLE DATASET ==============================================================

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


relabeling_poseidon <- function(res,pixels){

  clustering <- POSEIDON::estimate_clust_poseidon(res,bicl_mat = FALSE)

  med_pix <- tapply(pixels$Y,    clustering$col_cl, median)
  sorted <- sort(med_pix, ind=TRUE)

  Data <- data.frame(old = names(sorted$x),new = 1:length(names(sorted$x)))

  map <- setNames(Data$new, Data$old)
  Y <- unname(map[as.character(clustering$col_cl)])

  uNX <- colnames(clustering$row_cl[[1]])
  uNY <- paste0("COL_CL#",Data$new)

  map <- setNames(uNY, uNX)
  NY <- unname(map[as.character(uNX)])

  clustering$col_cl <- Y
  clustering$Data <- Data
  colnames(clustering$row_cl[[1]]) <- NY
  colnames(clustering$row_cl[[2]]) <- NY
  colnames(clustering$row_cl[[3]]) <- NY

  return(clustering)
}







relabeling_pose <- function(res,pixels){

  clustering <- Poser::estimate_clust_pose(res = res,bicl_mat = FALSE)

  med_pix <- tapply(pixels$Y,    clustering$col_cl, median)
  sorted <- sort(med_pix, ind=TRUE)

  Data <- data.frame(old = names(sorted$x),new = 1:length(names(sorted$x)))

  map <- setNames(Data$new, Data$old)
  Y <- unname(map[as.character(clustering$col_cl)])

  uNX <- colnames(clustering$row_cl[[1]])
  uNY <- ifelse(Data$new<=9,
                paste0("COL_CL#0",Data$new), paste0("COL_CL#",Data$new))


  map <- setNames(uNY, uNX)
  NY <- unname(map[as.character(uNX)])

  clustering$col_cl <- Y
  clustering$Data <- Data
  colnames(clustering$row_cl) <- NY

  return(clustering)
}

