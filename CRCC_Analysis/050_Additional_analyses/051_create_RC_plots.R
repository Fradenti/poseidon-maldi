source("CRCC_Analysis/000_Auxiliary.R")
library(POSEIDON)

out <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Best_Poseidon.RDS")
K <- 30
L <- 40

nnList_cpp  <- readRDS("CRCC_Analysis/Data/list_ind_resolution_4123_cpp.RDS")
Gli <- val1 <- readRDS("CRCC_Analysis/Data/FinalMatrix_Glycans.RDS")
Lip <- val2 <- readRDS("CRCC_Analysis/Data/FinalMatrix_Lipids.RDS")
Pep <- val3 <- readRDS("CRCC_Analysis/Data/FinalMatrix_Peptides.RDS")
pixels      <- readRDS("CRCC_Analysis/Data/FinalMatrix_Coor.RDS")

cl         <- relabeling_poseidon(out,pixels)

median_Lip <- (apply(Lip,2,function(x) tapply(x,cl$col_cl,median)))
median_Gly <- (apply(Gli,2,function(x) tapply(x,cl$col_cl,median)))
median_Pep <- (apply(Pep,2,function(x) tapply(x,cl$col_cl,median)))

rownames(median_Lip) <- paste("Cluster",1:15)
rownames(median_Gly) <- paste("Cluster",1:15)
rownames(median_Pep) <- paste("Cluster",1:15)


RC_heats <- function(data, CC, ind_data){
  data_sub <- data[,cl$col_cl==CC]



ind <- sort(cl$row_cl[[ind_data]][,CC],ind=TRUE)

annotation_df <- data.frame(
  ROW_CL = factor(cl$row_cl[[ind_data]][,CC][ind$ix])
)
rownames(annotation_df) <- rownames(data[ind$ix,])


mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(15)

# Optional: Specify custom colors for your annotations
ann_colors <- list(
  COL_CL = mycolors
)
cluster_counts <- table(cl$row_cl[[ind_data]][,CC][ind$ix])
# cumsum gives the end index of each block: 10, 20, 30
gap_locations <- cumsum(cluster_counts)
# We remove the last one because we don't need a gap after the final row
gap_locations <- gap_locations[-length(gap_locations)]

ph_obj <- pheatmap(log(data_sub[ind$ix,]),
                   cluster_rows = FALSE,     # Keep rows in our defined order to match annotations
                   cluster_cols = FALSE,
                   gaps_row = gap_locations,
                   annotation_row = annotation_df # Add the annotation side bar
)
return(ph_obj)
}



# Glycans ------------------------------------------------------------------
data <- t(as.matrix(Gli))
dim(data)
dim(cl$row_cl[[1]])
for(CC in 1:15){

  pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Glycans/Pixels_of_CC#",CC,".pdf"),
      width = 8, height = 8)
  plot(pixels,pch=".")
  points(pixels[cl$col_cl == CC,])
  dev.off()
  pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Glycans/Row_clusters_of_CC#",CC,".pdf"), width = 10, height = 20)
  a1 = RC_heats(data,CC = CC,ind_data = 1)
  print(a1)
  dev.off()
  cat(CC)
}

# Lipids ------------------------------------------------------------------
data <- t(as.matrix(Lip))
dim(data)
dim(cl$row_cl[[2]])
for(CC in 1:15){

  pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Lipids/Pixels_of_CC#",CC,".pdf"),
      width = 8, height = 8)
  plot(pixels,pch=".")
  points(pixels[cl$col_cl == CC,])
  dev.off()
  pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Lipids/Row_clusters_of_CC#",CC,".pdf"), width = 10, height = 20)
  a1 = RC_heats(data,CC = CC,ind_data = 2)
  print(a1)
  dev.off()
  cat(CC)
}

# Pepties ------------------------------------------------------------------
data <- t(as.matrix(Pep))
dim(data)
dim(cl$row_cl[[3]])
for(CC in 1:15){

  pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Peptides/Pixels_of_CC#",CC,".pdf"),
      width = 8, height = 8)
  plot(pixels,pch=".")
  points(pixels[cl$col_cl == CC,])
  dev.off()
  pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Peptides/Row_clusters_of_CC#",CC,".pdf"), width = 10, height = 25)
  a1 = RC_heats(data,CC = CC,ind_data = 3)
  print(a1)
  dev.off()
  cat(CC)
}






lapply(out$Y_list,dim)
# -------------------------------------------------------------------------
cl <- relabeling_poseidon(out,pixels)
# -------------------------------------------------------------------------
data <- t(Gli)
MAT  <- matrix(NA,nrow(data),15)
for(i in 1:15){
  MAT[,i] <-  out$MKCD[cl$row_cl[[1]][,i],1,1]
}
rownames(MAT) <- rownames(data)
colnames(MAT) <- paste("CC#",1:15)
a = pheatmap::pheatmap(MAT[sort(MAT[,1],index=T)$ix,],cluster_rows = F,cluster_cols = FALSE,
                   display_numbers = TRUE,
                   number_color = "black",
                   fontsize_number = 4,
                   color = hcl.colors(300, "Spectral"), main = "N-Glycans")

pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Glycans_summaryCapitoli.pdf"), width = 10, height = 15)
print(a)
dev.off()
# -------------------------------------------------------------------------
data <- t(Lip)
MAT  <- matrix(NA,nrow(data),15)
for(i in 1:15){
  MAT[,i] <-  out$MKCD[cl$row_cl[[2]][,i],1,1]
}
rownames(MAT) <- rownames(data)
colnames(MAT) <- paste("CC#",1:15)
a = pheatmap::pheatmap(MAT[sort(MAT[,1],index=T)$ix,],cluster_rows = F,cluster_cols = FALSE,
                   display_numbers = TRUE,
                   number_color = "black",
                   fontsize_number = 4,
                   color = hcl.colors(300, "Spectral"), main = "Lipids")

pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Lipids_summaryCapitoli.pdf"), width = 10, height = 15)
print(a)
dev.off()
dev.off()
# -------------------------------------------------------------------------
data <- t(Pep)
MAT  <- matrix(NA,nrow(data),15)
for(i in 1:15){
  MAT[,i] <-  out$MKCD[cl$row_cl[[3]][,i],1,1]
}
rownames(MAT) <- rownames(data)
colnames(MAT) <- paste("CC#",1:15)
a = pheatmap::pheatmap(MAT[sort(MAT[,1],index=T)$ix,],cluster_rows = F,cluster_cols = FALSE,
                   display_numbers = TRUE,
                   number_color = "black",
                   fontsize_number = 4,
                   color = hcl.colors(300, "Spectral"), main = "Peptides")

pdf(file = paste0("CRCC_Analysis/050_Additional_analyses/RCs/Peptides_summaryCapitoli.pdf"), width = 10, height = 20)
print(a)
dev.off()
