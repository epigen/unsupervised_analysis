#### load all clustering results and metadata and determine internal cluster indices ####

#### load libraries
library("clusterCrit")
library("stats")
set.seed(42)

# helper function for BIC calculation
do_BIC <- function(x) {
    # Check if there are at least two unique values in the column
    if(length(unique(x)) < 2){
        # Return Inf if there's only one unique value
        return(Inf)
    } else {
        # Perform BIC calculation if there are two or more unique values
        return(BIC(lm(data_mtx[,i] ~ as.factor(x))))
    }
}

### configurations

# input
metadata_path <- snakemake@input[["metadata"]]
clusterings_path <- snakemake@input[["clusterings"]]
pca_path <- snakemake@input[["pca"]]
pca_var_path <- snakemake@input[["pca_var"]]

# output
result_path <- snakemake@output[["internal_indices"]]

# parameters
internal_index <- as.character(snakemake@params['internal_index']) #"Silhouette"
sample_proportion <- as.numeric(snakemake@params['sample_proportion']) #0.1
metadata_of_interest <- unlist(c(snakemake@params['metadata_of_interest']))

### load data
metadata <- read.csv(file=file.path(metadata_path), row.names=1, header=TRUE)
clusterings <- read.csv(file=file.path(clusterings_path), row.names=1, header=TRUE)
pca_var <- read.csv(file=file.path(pca_var_path), row.names=1, header=TRUE)


# load PCs that explain >90% of the variance in the data
# Find the first row where the cumulative sum is greater than 0.9
cumulative_sum <- cumsum(pca_var[,1])
PCn <- which(cumulative_sum > 0.9)[1]
# Get the column classes of the full dataframe, select cols to load and set classes of the columns not to load to "NULL"
full_classes <- sapply(read.csv(file.path(pca_path), nrows = 1), class)
cols_to_load <- names(full_classes)[1:PCn+1]
classes <- ifelse(names(full_classes) %in% c("sample_name", cols_to_load), full_classes, "NULL")
# load the selected PCs only (slow)
pca <- read.csv(file.path(pca_path), colClasses = classes, row.names=1, header=TRUE)

# subset metadata to metadata_of_interest
if(length(metadata_of_interest)==0){
    metadata <- metadata[,1,drop=FALSE]
}else{
    metadata <- metadata[,metadata_of_interest,drop=FALSE]
}

# transform metadata
na_cols <- c()
for (col in colnames(metadata)){
    # if NA or less than 2 unique values -> remove column and move on
    if (any(is.na(metadata[[col]])) | length(unique(metadata[[col]]))<2){
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

### determine internal indices using clusterCrit
indices_df <- data.frame(matrix(ncol = 1, nrow = ncol(clusterings), dimnames = list(colnames(clusterings), internal_index)))

# prepare data for caluclations
# data_mtx <- as.matrix(data)
clusterings[colnames(clusterings)] <- lapply(clusterings[colnames(clusterings)], function(x) as.integer(x))
data_mtx <- as.matrix(pca)
data_mtx <- data_mtx[sample(nrow(data_mtx), ceiling(sample_proportion * nrow(data_mtx))), ]
clusterings <- clusterings[rownames(data_mtx),,drop=FALSE]

# calculate internal cluster index
if(internal_index %in% c("Silhouette", "Calinski_Harabasz", "C_index", "Davies_Bouldin", "Dunn")){
    for(clust in colnames(clusterings)){        
        indices_df[clust,internal_index] <- intCriteria(traj=data_mtx, part=clusterings[[clust]], crit=c(internal_index))
    }
} else if(internal_index=="AIC"){ # -> NOT USED
    ### determine indices using AIC on top PC of PCA
    AIC_sum <- rep(0L, ncol(clusterings))
    
    for(i in 1:ncol(data_mtx)){
        AIC_sum <- AIC_sum + unlist(apply(clusterings,2,function(x) AIC(lm(data_mtx[,i]~as.factor(x)))))*pca_var[i,1]
    }
    indices_df$AIC <- AIC_sum
} else if(internal_index=="BIC"){
    ### determine indices using BIC on top PC of PCA
    BIC_sum <- rep(0L, ncol(clusterings))
    
    for(i in 1:ncol(pca)){
#         BIC_sum <- BIC_sum + unlist(apply(clusterings,2,function(x) BIC(lm(data_mtx[,i]~as.factor(x)))))*pca_var[i,1] # crashed in case of only 1 cluster
        BIC_sum <- BIC_sum + unlist(apply(clusterings, 2, do_BIC)) * pca_var[i, 1]
    }
    indices_df$BIC <- BIC_sum
}

### save results
write.csv(indices_df, file=result_path, row.names=TRUE)
