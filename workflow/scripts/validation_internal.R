#### load all clustering results and metadata and determine internal cluster indices ####

#### load libraries
library("clusterCrit")
library("stats")

### configurations

# input
data_path <- snakemake@input[["data"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/data/digits_data.csv"
metadata_path <- snakemake@input[["metadata"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/data/digits_labels.csv"
clustering_path <- snakemake@input[["clustering"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/metadata_clusterings.csv"
pca_path <- snakemake@input[["pca"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_default_data_small.csv"
pca_var_path <- snakemake@input[["pca_var"]] #"/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_default_2_var.csv"

# output
result_path <- snakemake@output[["internal_indices"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/internal_indices.csv"

# parameters
samples_by_features <- as.integer(snakemake@params['samples_by_features']) #1

### load data
data <- read.csv(file=file.path(data_path), row.names=1, header=TRUE)
metadata <- read.csv(file=file.path(metadata_path), row.names=1, header=TRUE)
clusterings <- read.csv(file=file.path(clustering_path), row.names=1, header=TRUE)
pca <- read.csv(file=file.path(pca_path), row.names=1, header=TRUE)
pca_var <- read.csv(file=file.path(pca_var_path), row.names=1, header=TRUE)

# check and fix orientation
if(samples_by_features==0){
    data <- t(data)
}

# transform metadata
na_cols <- c()
for (col in colnames(metadata)){
    # if NA -> remove column and move on
    if (any(is.na(metadata[[col]]))){
        na_cols <- c(na_cols, col)
        next
    }
    # if a metadata class is empty ("") fill with "unknown"
    if (!any(is.na(metadata[[col]]))){
        if (any(metadata[[col]]=="")){
            metadata[metadata[[col]]=="", col] <- "unknown"
        }
    }
    # convert metadata to categorical if less than 25 unique integer values
    if (is.numeric(metadata[[col]]) & length(unique(metadata[[col]]))<=25){
        if(all(metadata[[col]] == round(metadata[[col]]))){
            metadata[col] <- as.factor(metadata[[col]])
        }
    }
}
# remove columns with NA
metadata <- metadata[, !(colnames(metadata) %in% na_cols),drop=FALSE]
# add categorical metadata to clustering results with prefix "metadata_"
metadata_cat <- metadata[,sapply(metadata, function(x) !is.numeric(x)), drop=FALSE]
# Convert all categorical columns to integer
metadata_cat[colnames(metadata_cat)] <- lapply(metadata_cat[colnames(metadata_cat)], function(x) as.integer(factor(x)))
colnames(metadata_cat) <- paste0("metadata_", colnames(metadata_cat))                  
clusterings <- cbind(clusterings, metadata_cat)


### determine indices using clusterCrit
indices_df <- data.frame()
for(clust in colnames(clusterings)){
    indices_tmp <- intCriteria(traj=as.matrix(data), part=as.integer(clusterings[[clust]]), crit=c("Silhouette", "Calinski_Harabasz", "C_index", "Davies_Bouldin", "Dunn"))
    indices_tmp_df <- as.data.frame(t(unlist(indices_tmp)))
    indices_df <- rbind(indices_df, indices_tmp_df)
}
rownames(indices_df) <- colnames(clusterings)
                                 
                                 
### determine indices using AIC & BIC on top PC of PCA
AIC_sum <- rep(0L, ncol(clusterings))
BIC_sum <- rep(0L, ncol(clusterings))

for(i in 1:ncol(pca)){
    AIC_sum <- AIC_sum + unlist(apply(clusterings,2,function(x) AIC(lm(pca[,i]~as.factor(x)))))*pca_var[i,1]
    BIC_sum <- BIC_sum + unlist(apply(clusterings,2,function(x) BIC(lm(pca[,i]~as.factor(x)))))*pca_var[i,1]
}
                                      
indices_df$AIC <- AIC_sum
indices_df$BIC <- BIC_sum
                                 

### save results
write.csv(indices_df, file=result_path, row.names=TRUE)
