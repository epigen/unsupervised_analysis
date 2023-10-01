#### load libraries
library("ggplot2")
library("patchwork")
library("reshape2")
# library(dendextend)

### configurations

# input
index_paths <- snakemake@input

# output
plot_path <- snakemake@output[["plot"]]

result_dir <- file.path(dirname(plot_path))
# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

indices <- names(index_paths)
indices <- indices[indices!=""]

heatmaps <- list()

# loop through all proivided index paths and create (clustered) heatmap
for (idx in indices){
    
    scores <- read.csv(file.path(index_paths[[idx]]), row.names = 1)
    
    # scale internal indices, add seperator and reorder columns (max/sep/min)
    if(length(indices)==1){
        scores <- as.data.frame(scale(scores))
        scores$sep <- NA
        scores <- scores[,c("Silhouette", "Calinski_Harabasz", "Dunn", "sep", "C_index", "Davies_Bouldin", "BIC")]
    }

    height <- nrow(scores) * 0.5
    width <- ncol(scores)*0.5 + 5

    # Perform hierarchical clustering on the rows and columns of external indices
    if(nrow(scores)>2 & ncol(scores)>2 & length(indices)!=1){
        cluster_rows <- hclust(dist(scores))
        cluster_cols <- hclust(dist(t(scores)))
        scores <- scores[order.dendrogram(as.dendrogram(cluster_rows)), order.dendrogram(as.dendrogram(cluster_cols))]
    }

    # prepare plotting dataframe
    scores$clusterings <- rownames(scores)
    scores_plot <- melt(scores, id.vars = "clusterings", variable.name = "metadata")
    
    # retain MCDM ranking
    if(length(indices)==1){
        scores_plot$clusterings <- factor(scores_plot$clusterings, levels = rev(unique(scores_plot$clusterings)))
    }
    
    # Create a heatmap
    heatmaps[[idx]] <- ggplot(scores_plot, aes(x=metadata, y=clusterings, fill=value)) + 
        geom_tile(color = "white", size = 1) + 
        geom_text(aes(label=round(value, 2)), size=3) +
        theme_minimal() + 
        labs(title=gsub("_", " ", idx)) +
        {if(length(indices)==1) xlab("internal cluster indices (benefit/max || cost/min)") } +
        {if(length(indices)==1)  scale_x_discrete(labels = function(x) gsub(".*\\.", "", x), breaks = function(x) {x[x!="sep"]})  } +
        scale_fill_gradient2(midpoint=0, low="royalblue4", mid="white", high="firebrick2", space ="Lab") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels by 45 degrees
    
#     options(repr.plot.width=width, repr.plot.height=height)
#     print(heatmaps[[idx]])
}

# differentiate between internal (one heatmap) and external (panel of six heatmaps) indices
if(length(heatmaps)==1){
    heatmap_panel <- heatmaps[[1]]
    width_panel <- width
    height_panel <- height
}else{
    heatmap_panel <- wrap_plots(heatmaps, ncol=3)
    width_panel <- 3 * width
    height_panel <- ceiling(length(indices)/3) * height
}

### save plot
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# heatmap_panel

ggsave(basename(plot_path),
       plot = heatmap_panel,
       device = 'png',
       path = result_dir,
       scale = 1,
       dpi = 300,
       width = width_panel,
       height = height_panel,
       limitsize = FALSE,
      )
