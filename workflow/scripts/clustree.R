#### load libraries
library("clustree")
library("patchwork")

### configurations

# input
clustering_path <- snakemake@input[["metadata_clustering"]]
metadata_path <- snakemake@input[["metadata"]]

# output
plot_path <- snakemake@output[["plot"]]

# parameters
content <- snakemake@wildcards[["content"]] # wildcard, one of: metadata, features, default, custom
categorical_label_option <- as.character(snakemake@params[["categorical_label_option"]]) #"pure" or "majority"
numerical_aggregation_option <- as.character(snakemake@params[["numerical_aggregation_option"]]) # mean, median, max, min
custom_metadata <- c(snakemake@params[["custom_metadata"]]) # c("target")
count_filter <- as.numeric(snakemake@params[["count_filter"]]) # 0
prop_filter <- as.numeric(snakemake@params[["prop_filter"]]) # 0.1
layout <- as.character(snakemake@params[["layout"]])  # "tree" or "sugiyama"


if (content=="default" | content=="custom"){
    result_dir <- file.path(dirname(plot_path))
}else{
    result_dir <- file.path(plot_path)
}
# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}


# helper function for labeling catgorical metadata
categorical_labeler <- function(labels) {
    
    if (categorical_label_option == "majority"){
        label <- as.character(names(which.max(table(labels))))
    } else if (categorical_label_option == "pure"){
        if (length(unique(labels)) == 1) {
            label <- as.character(unique(labels))
        } else {
            label <- "mixed"
        }
    }
    return(label)
}

# helper function for aggregation of numerical metadata
numerical_aggregation <- function(values){
    if(numerical_aggregation_option=="mean"){
        label <- mean(values)
    } else if (numerical_aggregation_option=="median"){
        label <- median(values)
    } else if (numerical_aggregation_option=="min"){
        label <- min(values)
    } else if (numerical_aggregation_option=="max"){
        label <- max(values)
    }
    
    return(round(label,1))
}

plot_clustree <- function(data, col) {
    clustree_plot <- clustree(
      x = data,
      prefix = "X_",
      suffix = NULL,
      #metadata = features, # not used in case of dataframe as input
      count_filter = count_filter,
      prop_filter = prop_filter,
      layout = layout,
      use_core_edges = TRUE,
      highlight_core = FALSE, # check effect -> TODO
      node_colour = col, # depending on col
      node_colour_aggr = if(col=="X_") NULL else if(is.numeric(data[[col]])) "numerical_aggregation" else "categorical_labeler", # depending on col
      node_size = "size",
      node_size_aggr = NULL,
      node_size_range = c(4, 15),
      node_alpha = 1,
      node_alpha_aggr = NULL,
      node_text_size = 3,
      scale_node_text = FALSE,
      node_text_colour = "black",
      node_label = if(col=="X_") NULL else col, # depending on col
      node_label_aggr = if(col=="X_") NULL else if(is.numeric(data[[col]])) "numerical_aggregation" else "categorical_labeler", # depending on col
      node_label_size = 3,
      node_label_nudge = -0.2,
      edge_width = 1.5,
      edge_arrow = FALSE, # default: TRUE, but makes plot unnecessary busy
      edge_arrow_ends = "last",
      show_axis = FALSE,
      return = "plot"
    ) + 
      {if(col=="X_") guides(color="none")} +  # remove legend in case of color being clusterings
                         {if(is.factor(data[[col]]) & length(unique(data[[col]]))>25) guides(color="none")} +  # remove legend in case of color being categorical with more than 25 categories
                         scale_edge_color_continuous(low = "grey80", high = "firebrick2") +
                         {if(is.numeric(data[[col]])) scale_color_gradient2(midpoint=0, low="royalblue4", mid="grey80", high="firebrick2", space ="Lab")} + 
                         {if(is.numeric(data[[col]])) scale_fill_gradient2(midpoint=0, low="royalblue4", mid="grey80", high="firebrick2", space ="Lab")} +
                         labs(color=col) + guides(fill="none")+
                         geom_label(data = clustering_names, aes(x = Inf, y = index-1+0.25, label = clustering), hjust = 1, vjust = 0.5, size = 3, inherit.aes = FALSE) + # add clustering names
    {if(col=="X_") ggtitle(content) else ggtitle(col)} + theme(plot.title = element_text(size = 10))

    return(clustree_plot)
}


### load data
clusterings <- read.csv(file=file.path(clustering_path), row.names=1, header=TRUE)
metadata <- read.csv(file=file.path(metadata_path), row.names=1, header=TRUE)

### transform data

# sort (rows by clustering rownames and) columns by alphabet -> rows should be sorted correctly in all cases
metadata <- metadata[,sort(colnames(metadata)),drop=FALSE] # rownames(clusterings)

# transform metadata
if (content=="metadata"){
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
} else if (content=="custom"){# for custom plot with selected categorical metadata as clustering
    # subset to selected metadata & add prefix to colnames "metadata_"
    metadata <- metadata[,custom_metadata, drop=FALSE]
    colnames(metadata) <- paste0("metadata_", colnames(metadata))
    # convert to factor
    for (col in colnames(metadata)){
        metadata[col] <- as.factor(metadata[[col]])
    }
    
    # add to clusterings
    clusterings <- cbind(clusterings, metadata)
}

# sort by number of clusters
num_unique <- sapply(clusterings, function(x) length(unique(x)))
# Get the order of the columns sorted by the number of unique elements
order_num_unique <- order(num_unique)
# Sort the dataframe by the columns in this order
clusterings <- clusterings[, order_num_unique]

# save the mapping from column index to column name
clustering_names <- data.frame(clustering = colnames(clusterings), index = ncol(clusterings):1)

# rename columns to X_{n}
colnames(clusterings) <- paste0("X_", 1:ncol(clusterings))
                     
# plot specifications
width <- ceiling(length(unique(clusterings[,ncol(clusterings)]))/3) + 2 # per clustering 1/3 inch
height <- ncol(clusterings) + 2
                                          
# add metadata and features to clusterings
data <- cbind(clusterings, metadata)#, features)

                     
### clustree analysis
if (content=="default" | content=="custom"){
    # make default or custom plot without any metadata highlighted
    clustree_final <- plot_clustree(data, "X_")
    
    ggsave(basename(plot_path),
       plot = clustree_final,
       device = 'png',
       path = result_dir,
       scale = 1,
       dpi = 300,
       width = width,
       height = height,
       limitsize = FALSE,
      )
} else{
    # generate & save metadata/features plots
    for (col in colnames(metadata)){
        
        clustree_final <- plot_clustree(data, col)
        
        ggsave(paste0("clustree_",content,"_",col,".png"),
           plot = clustree_final,
           device = 'png',
           path = result_dir,
           scale = 1,
           dpi = 300,
           width = width,
           height = height,
           limitsize = FALSE,
          )
    }
}

