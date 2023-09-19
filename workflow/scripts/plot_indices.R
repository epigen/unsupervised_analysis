#### load libraries
library("ggplot2")
library("patchwork")
library("reshape2")
# library(dendextend)

### configurations

# input
index_paths <- snakemake@input # list("AMI"='/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_AMI.csv', "ARI"='/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_ARI.csv', "FMI"='/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_FMI.csv')

# output
plot_path <- snakemake@output[["plot"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_indices.png"

# parameters
# content <- snakemake@wildcards[["content"]] 

result_dir <- file.path(dirname(plot_path))
# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

indices <- names(index_paths)
indices <- indices[indices!=""]

# print(indices)

heatmaps <- list()

# loop through all proivided index paths and create (clustered) heatmap
for (idx in indices){
#     print(index_paths[[idx]])
    scores <- read.csv(file.path(index_paths[[idx]]), row.names = 1)

    height <- nrow(scores) * 0.5
    width <- ncol(scores)*0.5 + 5

    # Perform hierarchical clustering on the rows and columns of the matrix
    if(nrow(scores)>2 & ncol(scores)>2){
        cluster_rows <- hclust(dist(scores))
        cluster_cols <- hclust(dist(t(scores)))
        scores <- scores[order.dendrogram(as.dendrogram(cluster_rows)), order.dendrogram(as.dendrogram(cluster_cols))]
    }

    scores$clusterings <- rownames(scores)
    scores_plot <- melt(scores, id.vars = "clusterings", variable.name = "metadata")

    # Create a heatmap
    heatmaps[[idx]] <- ggplot(scores_plot, aes(x=metadata, y=clusterings, fill=value)) + 
        geom_tile(color = "white", size = 1) + 
        geom_text(aes(label=round(value, 2)), size=3) +
        theme_minimal() + 
        labs(title=idx) +
        scale_fill_gradient2(midpoint=0, low="royalblue4", mid="white", high="firebrick2", space ="Lab") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels by 45 degrees
    
#     options(repr.plot.width=width, repr.plot.height=height)
#     print(heatmaps[[idx]])
}


# determine plot and panel size
# white borders?
# add dendrogram?

heatmap_panel <- wrap_plots(heatmaps, ncol=3)

width_panel <- 3 * width
height_panel <- ceiling(length(indices)/3) * height


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
