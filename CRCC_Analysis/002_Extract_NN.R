library(raster)
library(tidyverse)

# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------
# 4000 pixels -------------------------------------------------------------


pix <- readRDS("CRCC_Analysis/Data/FinalMatrix_Coor.RDS") + 1 
colnames(pix) <- c("X","Y")
number_of_pix <- nrow(pix)

summary(pix)
plot(pix)
plot(pix,pch=".")
dim(pix)
# get NN
rast <- raster::raster(matrix(1,max(pix$X),max(pix$Y)))
ADJ  <- adjacent(rast,
                 cells = 1:ncell(rast),
                 directions=8,sorted=T)
pixels <- expand_grid(1:max(pix$X),1:max(pix$Y))

nnList = list()
for(i in 1:ncell(rast)){
  nnList[[i]] = ADJ[which(ADJ[,1]==i),2]
  if( i %%500 ==0)cat(i)
}

nnList[[1]]
Dat = as_tibble(cbind(pixels , raster::values(rast)))
nnList_cpp = lapply(nnList,function(x)x-1)


plot(pixels, pch=".")
points(pix, col=2)
points(pixels[1,], col=5,pch=21,bg=2)
points(pixels[nnList[[1]],], pch="x") #
# these up here are the nn on the square - we need to map it to the new reference values

# -------------------------------------------------------------------------
# remove problematic pixels (no neighbors)
nnList <- lapply(nnList_cpp,function(x)x+1)
IND_all <- pixels

POTTS_list  = list()
orig_px_num = numeric(nrow(pix))
problematic_pix = numeric(2000) ; q=0

for(i in 1:nrow(pix)){
  
  found = which(pix[i,2]==IND_all[,2] & pix[i,1]==IND_all[,1])
  if(length(found) != 1){
    POTTS_list[[i]] =NA
    orig_px_num[i] = NA
    q=q+1
    problematic_pix[q] = i
  } else{
    POTTS_list[[i]] = nnList[[found]]
    orig_px_num[i] = found
  }
  if(i %%100 == 0 )cat(paste(i,"\n"))
}

problematic_pix
sum(is.na(orig_px_num))



##


POTTS_list_sub = POTTS_list
for(i in 1:nrow(pix)){
  POTTS_list_sub[[i]] =
    POTTS_list[[i]][POTTS_list[[i]] %in% orig_px_num]
  if(i %%100 == 0 )cat(paste(i,"\n"))
}


LIST_INDEX_PICS = list()
for(i in 1:nrow(pix)){
  
  pixe = IND_all[POTTS_list_sub[[i]],]
  if(is.null(nrow(pix))){
    LIST_INDEX_PICS[[i]] = which(pixe[1]== pix[,1] & pixe[2]==pix[,2])
  }else{
    LIST_INDEX_PICS[[i]] = apply(pixe,1,function(x) which(c(x)[1]==pix[,1] &
                                                            c(x)[2]==pix[,2]))
  }
  if(i %% 500 ==0) cat(i)
}

plot(pix,pch=".")
q = 10
points(pix[q,2]~pix[q,1],pch="x")
points(pix[LIST_INDEX_PICS[[q]],],pch="x",col=4)
name1 <- paste0("CRCC_Analysis/Data/list_ind_resolution_",number_of_pix,".RDS")
saveRDS(LIST_INDEX_PICS,name1)

name2 <- paste0("CRCC_Analysis/Data/list_ind_resolution_",number_of_pix,"_cpp.RDS")

LIST_INDEX_PICS_cpp = lapply(LIST_INDEX_PICS,function(x)x-1)
saveRDS(LIST_INDEX_PICS_cpp,name2)



plot(pix,pch=".")
q = 11
points(pix[q,2]~pix[q,1],pch="x")
points(pix[LIST_INDEX_PICS[[q]],],pch="x",col=4)
plot(pix,pch=".")
q = 1000
points(pix[q,2]~pix[q,1],pch="x")
points(pix[LIST_INDEX_PICS_cpp[[q]]+1,],pch="x",col=4)

