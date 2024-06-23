#### libraries
library("ggplot2")
library("patchwork")
library("data.table")

# utility function to adapt legend according to metadata
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 3, spaceLegend = 0, alpha = 1) {
    new_plot <- myPlot +
        guides(color = guide_legend(override.aes = list(size = pointSize, alpha = alpha))) +
        theme(legend.text  = element_text(size = textSize),
              legend.key.size = unit(spaceLegend, "lines"))
    return (new_plot)
}

### configurations

# input
data_path <- snakemake@input[["dimred_data"]]
axes_path <- snakemake@input[["dimred_axes"]]
metadata_path <- snakemake@input[["metadata"]]

# output
plot_path <- snakemake@output[["plot"]]

# parameters
size <- snakemake@params[["size"]]# 0.5
alpha <- snakemake@params[["alpha"]]# 1
coord_fixed_flag <- if(as.integer(snakemake@config[["coord_fixed"]])==1) TRUE else FALSE

dir.create(plot_path, recursive = TRUE)

### load data
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)[,1:2]
axes <- data.frame(fread(file.path(axes_path), header=TRUE), row.names=1)[1:2,]
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)

# plot specifications
width <- 4
height <- 3
shape <- if(nrow(data)>50000) '.' else 16

### make plots

for (col in sort(colnames(metadata))){
    print(col)
    
    # check if metadata column is only NA
    if(all(is.na(metadata[[col]]))){
        next
    }
    
    # convert to categorical if the substring "cluster" is present
    if (grepl("cluster", col)) {
        metadata[col] <- as.factor(metadata[[col]])
    }
    
    # convert metadata to categorical if less than 25 unique integer values
    if (grepl("metadata", basename(plot_path)) & is.numeric(metadata[[col]]) & length(unique(metadata[[col]]))<=25){
        if(all(metadata[[col]] == round(metadata[[col]]))){
            metadata[col] <- as.factor(metadata[[col]])
        }
    }
    
    # if a metadata class is empty ("") fill with "unknown"
    if (!any(is.na(metadata[[col]]))){
        if (any(metadata[[col]]=="")){
            metadata[metadata[[col]]=="", col] <- "unknown"
        }
    }
    
    # prepare data for plotting
    tmp_data <- cbind(data, metadata[col])
    
    # make 2D scatter plots
    if (!is.numeric(metadata[[col]])){
        # plot categorical data
        tmp_plot <- ggplot(tmp_data, aes_string(x=colnames(tmp_data)[1], y=colnames(tmp_data)[2])) +
        geom_point(aes_string(color=col), size=size, stroke=0, alpha=alpha, shape=shape) + 
        {if(coord_fixed_flag) coord_fixed()} +
        xlab(axes[1]) +    
        ylab(axes[2]) +
        ggtitle(col) +
        theme_linedraw() + 
        theme(plot.title = element_text(size = 10), legend.title = element_blank())
        
        if (length(unique(metadata[[col]]))>25){
            tmp_plot <- tmp_plot + theme(legend.position="none") #addSmallLegend(myPlot = tmp_plot, pointSize = 0, textSize = 0)
        }else if (length(unique(metadata[[col]]))>15){
            tmp_plot <- addSmallLegend(myPlot = tmp_plot, textSize = 4)
        }else{
            tmp_plot <- addSmallLegend(myPlot = tmp_plot, textSize = 10)
        }
    }else{
        # plot numerical data
        tmp_plot <- ggplot(tmp_data, aes_string(x=colnames(tmp_data)[1], y=colnames(tmp_data)[2])) +
        geom_point(aes_string(color=col), size=size, stroke=0, alpha=alpha, shape=shape) + 
        {if(coord_fixed_flag) coord_fixed()} +
        scale_color_gradient2(midpoint=0, low="royalblue4", mid="grey80", high="firebrick2", space ="Lab") +
        xlab(axes[1]) +    
        ylab(axes[2]) +
        ggtitle(col) +
        theme_linedraw() + 
        theme(plot.title = element_text(size = 10), legend.title = element_blank())
        
    }
    
    # save plot
    ggsave(paste0(col,".png"),
       plot = tmp_plot,
       device = 'png',
       path = plot_path,
       scale = 1,
       dpi = 300,
       width = width,
       height = height,
       limitsize = FALSE,
      )
}
