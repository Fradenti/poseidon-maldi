out <- readRDS("CRCC_Analysis/Output/RDS/All_runs_Poseidon/Best_Poseidon.RDS")

cl <- relabeling_poseidon(out,pixels)

median_Lip <- (apply(Lip,2,function(x) tapply(x,cl$col_cl,median)))
median_Gly <- (apply(Gli,2,function(x) tapply(x,cl$col_cl,median)))
median_Pep <- (apply(Pep,2,function(x) tapply(x,cl$col_cl,median)))

rownames(median_Lip) <- paste("Cluster",1:15)
rownames(median_Gly) <- paste("Cluster",1:15)
rownames(median_Pep) <- paste("Cluster",1:15)

write_csv2(as.data.frame(median_Lip),file = "CRCC_Analysis/050_Additional_analyses/Median_lip.csv")
write_csv2(as.data.frame(median_Gly),file = "CRCC_Analysis/050_Additional_analyses/Median_Gly.csv")
write_csv2(as.data.frame(median_Pep),file = "CRCC_Analysis/050_Additional_analyses/Median_Pep.csv")


# -------------------------------------------------------------------------


library(pheatmap)
# --- 1. Create sample data ---
heat <- function(data, title){

mat <- data
rownames(mat) <- paste0("Pixel", 1:nrow(mat))
colnames(mat)
ind <- sort(cl$col_cl,ind=TRUE)
# --- 2. Create row annotations ---
# This data frame must have rownames that match the rownames of your matrix
annotation_df <- data.frame(
  COL_CL = factor(cl$col_cl[ind$ix])
)
rownames(annotation_df) <- rownames(mat[ind$ix,])


mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(15)

names(mycolors) = as.character(1:15)
# Optional: Specify custom colors for your annotations
ann_colors <- list(
  COL_CL = mycolors
)

cluster_counts <- table(cl$col_cl[ind$ix])
# cumsum gives the end index of each block: 10, 20, 30
gap_locations <- cumsum(cluster_counts)
# We remove the last one because we don't need a gap after the final row
gap_locations <- gap_locations[-length(gap_locations)]

# Check indices:
print(gap_locations)
# --- 3. Plot with annotation_row ---
ph_obj <- pheatmap(log(mat[ind$ix,]),
         cluster_rows = FALSE,     # Keep rows in our defined order to match annotations
         cluster_cols = FALSE,
         #filename = paste0("CRCC_Analysis/050_Additional_analyses/",title,"_heatmap.pdf"),
         annotation_row = annotation_df, # Add the annotation side bar
         annotation_colors = ann_colors, # Add custom colors
         gaps_row = gap_locations,
         show_rownames = FALSE,     # Show the actual text labels
         main = title )
return(ph_obj)
}


pdf(file = "CRCC_Analysis/050_Additional_analyses/Lipids_heatmap.pdf", width = 20, height = 10)
a1 = heat(Lip,"Lipids")
print(a1)
dev.off()
pdf(file = "CRCC_Analysis/050_Additional_analyses/Glycans_heatmap.pdf", width = 20, height = 10)
a2 = heat(Gli,"N-Glycans")
print(a2)
dev.off()
pdf(file = "CRCC_Analysis/050_Additional_analyses/Peptides_heatmap.pdf", width = 20, height = 10)
a3 = heat(Pep,"Peptides")
print(a3)
dev.off()

