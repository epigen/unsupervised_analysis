#### load libraries
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
#library("dendsort")

### configurations
set.seed(42)
ht_opt(fast_hclust = TRUE)

# inputs
data_path <- snakemake@input[["data"]]
metadata_path <- snakemake@input[["metadata"]]

# output
plot_path <- snakemake@output[["plot"]]

# parameters
samples_by_features <- as.integer(snakemake@params['samples_by_features']) #1
metric <- snakemake@params[["metric"]]# "pearson"
cluster_method <- snakemake@params[["cluster_method"]]# "complete"
metadata_col <- c(snakemake@config[["metadata_of_interest"]])[1]


result_dir <- file.path(dirname(plot_path))
# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

### load data
data <- read.csv(file=file.path(data_path), row.names=1, header=TRUE)
metadata <- read.csv(file=file.path(metadata_path), row.names=1, header=TRUE)

# check and fix orientation
if(samples_by_features==0){
    data <- t(data)
}

# scale data
data <- scale(data)

# remove NA columns & replace NA
data <- data[,colSums(is.na(data))<nrow(data)]
data[is.na(data)] <- 0

# prepare metadata
if(is.null(metadata_col)){ #|!(metadata_col %in% colnames(metadata))){
    metadata_col <- colnames(metadata)[1]
}

# check if metadata column is only NA and switch to the first that is not
if(all(is.na(metadata[[metadata_col]]))){
    for(col in colnames(metadata)){
        if(all(is.na(metadata[[col]]))){
            next
        }else{
            metadata_col <- col
            break
        }
    }
}

if (is.numeric(metadata[[metadata_col]]) & length(unique(metadata[[metadata_col]]))<=25){
    if(all(metadata[[metadata_col]] == round(metadata[[metadata_col]]))){
        metadata[metadata_col] <- as.factor(metadata[[metadata_col]])
    }
}
# if a metadata class is empty ("") fill with "unknown"
if (!any(is.na(metadata[[metadata_col]]))){
    if (any(metadata[[metadata_col]]=="")){
        metadata[metadata[[metadata_col]]=="", metadata_col] <- "unknown"
    }
}

# replace NA values
if (any(is.na(metadata[[metadata_col]]))){
    if (is.numeric(metadata[[metadata_col]])){
        metadata[is.na(metadata[[metadata_col]]),metadata_col] <- 0
    }else{
        metadata[is.na(metadata[[metadata_col]]),metadata_col] <- "NA"
    }
}


# plot specifications

# make colors for values
limit <- ceiling(max(abs(quantile(data, c(0.01, 0.99)))))
col_fun <- colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))

# make color mapping
if (is.numeric(metadata[[metadata_col]])){
    #check if divergent or sequential?
    row_annot <- HeatmapAnnotation(df=metadata[,metadata_col, drop=FALSE], which="row")
    
}else{
    n_cat <- length(unique(metadata[[metadata_col]]))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    colors <- sample(unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))),n_cat)
    names(colors) <- unique(metadata[[metadata_col]])
    colors_list <- list()
    colors_list[[metadata_col]] <- colors# put here the mapped colors
    
    row_annot <- HeatmapAnnotation(df=metadata[,metadata_col, drop=FALSE], which="row", col=colors_list)
}


# limit labeling to 100 in each dim for readability
show_row_names <- ifelse(nrow(data)>100, FALSE, TRUE)
show_column_names <- ifelse(ncol(data)>100, FALSE, TRUE)

# # alternative way to: determine distance, hierarchical clustering and order dendrograms
# # would be used in the parameters cluster_rows=row_dend and cluster_columns=col_dend as arguments
# if(metric %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
#     row_dend <- dendsort(hclust(dist(data, method=metric), method=cluster_method))
#     col_dend <- dendsort(hclust(dist(t(data), method=metric), method=cluster_method))
# } else if(metric %in% c("pearson", "spearman", "kendall")){
#     #rows
#     row_data <- cor(t(data), method = metric)
#     row_data[] <- sapply(row_data, function(x) 1-abs(x))
#     row_data <- as.dist(row_data)
#     row_dend <- dendsort(hclust(row_data, method=cluster_method))
    
#     #cols
#     col_data <- cor(data, method = metric)
#     col_data[] <- sapply(col_data, function(x) 1-abs(x))
#     col_data <- as.dist(col_data)
#     col_dend <- dendsort(hclust(col_data, method=cluster_method))
# } else {
#     print(paste0("metric not supported: ", metric))
# }

### make & save heatmap
# options(repr.plot.width=10, repr.plot.height=10)
png(filename=plot_path, width=20, height=20, units = "in", res=300)

Heatmap(data,
        name = "z-scores",
        column_title = paste0("Heatmap of data scaled by features (z-scores), hierarchically clustered using method ",cluster_method," with distance metric ",metric,", and colorscape limited to the top percentiles."),
        col = col_fun,
        left_annotation = row_annot,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        clustering_distance_rows = metric,
        clustering_method_rows = cluster_method,
        clustering_distance_columns = metric,
        clustering_method_columns = cluster_method,
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE,
        use_raster = TRUE,
        raster_quality = 9
       )

dev.off()
