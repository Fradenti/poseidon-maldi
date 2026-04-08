source("A_Aux_Function_1.R")
source("A_Aux_Function_2.R")
source("B_Hyperpar.R")

#BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
#BiocManager::install("scater")
library(scater)
#BiocManager::install("BayesSpace",force=TRUE)
library(BayesSpace)

library(mcclust)

# -------------------------------------------------------------------------
# Simulation Script with i_sim input and result saving
# -------------------------------------------------------------------------
# i_sim <- 5

# Accept i_sim from environment (or default to 1)
i_sim <- as.integer(Sys.getenv("i_sim", unset = "1"))
i_str <- as.integer(Sys.getenv("i_str", unset = "1"))

cat("Running simulation number:", i_sim, " with i_str ", i_str,"\n")

str <- strength[i_str]

full_data = read.csv(file = paste0("data/simulated_dataC/simulated_data_str",str,"_",i_sim,".csv"))
metadata_cols = c('cell_id', 'x', 'y', 'true_cluster', 'sample_id')
expression_cols = setdiff(colnames(full_data),metadata_cols)
X = apply(full_data[,expression_cols], 2, as.numeric)
obs = full_data[,metadata_cols]
obs$true_cluster <- as.numeric(obs$true_cluster)

kClusters <- length(unique(obs$true_cluster))

tm <- Sys.time()

sce <- SingleCellExperiment(assays = list(logcounts = Matrix::Matrix(t(X), sparse = TRUE)))
colData(sce)$array_row <- as.integer(obs$y)
colData(sce)$array_col <- as.integer(obs$x)

sce <- BayesSpace::spatialPreprocess(sce, platform="ST",
                                     skip.PCA=TRUE) # we need to do PCA manually

kPCA = 15 # default recommended by BayesSpace
sce <- scater::runPCA(sce, ntop = nrow(sce),
                      ncomponents = kPCA, exprs_values="logcounts")

sce <- spatialCluster(sce, q=kClusters, platform="ST", d=kPCA,
                      init.method="mclust", model="t", gamma=2,
                      nrep=1000, burn.in=100,
                      save.chain=FALSE)

tot_time <- Sys.time() - tm

est_colcl <- colData(sce)$spatial.cluster
ARI_colcl <- mcclust::arandi(obs$true_cluster,est_colcl)
VI_colcl <- mcclust::vi.dist(obs$true_cluster,est_colcl)
est_bicl <- NA
RMSE_bicl <- NA
ARI_bicl <- NA

save_list <- c("i_sim",
               "sce","kClusters", "ARI_colcl", "VI_colcl", "RMSE_bicl", "ARI_bicl",
               "est_colcl", "est_bicl", "tot_time")
save_file <- paste0("results/RData/BayesSpace/BayesSpace_result_str",str,"_", i_sim, ".RData")
save(list = save_list, file = save_file)
