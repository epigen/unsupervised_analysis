#### load libraries
library("ggplot2")
library("patchwork")

# utility function to adapt legend according to metadata
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 3, spaceLegend = 0, alpha = 1) {
    new_plot <- myPlot +
        guides(color = guide_legend(override.aes = list(size = pointSize, alpha = alpha))) +
        theme(legend.text  = element_text(size = textSize),
              legend.key.size = unit(spaceLegend, "lines"))
    return (new_plot)
}

### configurations

# inputs
data_path <- snakemake@input[["dimred_data"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/PCA_data.csv"
axes_path <- snakemake@input[["dimred_axes"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/PCA_axes.csv"
metadata_path <- snakemake@input[["metadata"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/KOcall_NonTargeting/counts/CORRECTED_metadata.csv"

plot_path <- snakemake@output[["plot"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/plots/PCA_metadata.png"

size <- snakemake@params[["size"]]# 0.5
alpha <- snakemake@params[["alpha"]]# 1

result_dir <- file.path(dirname(plot_path))
# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

# parameters
# step <- snakemake@params[["step"]] #""


### load data
data <- read.csv(file=file.path(data_path), row.names=1, header=TRUE)[,1:2]
axes <- read.csv(file=file.path(axes_path), row.names=1, header=TRUE)[1:2,]
metadata <- read.csv(file=file.path(metadata_path), row.names=1, header=TRUE)

# plot specifications
n_col <- min(10, ncol(metadata))
width <- 5
height <- 3

width_panel <- n_col * width
height_panel <- ceiling(ncol(metadata)/n_col) * height

# col <- colnames(metadata)[1]
# col <- colnames(metadata)[2]

### make plots
scatter_plots <- list()

for (col in sort(colnames(metadata))){
    print(col)
    
    # check if metadata column is only NA
    if(all(is.na(metadata[[col]]))){
        next
    }
    
    # convert to categorical if less than 25 unique integer values
    if (is.numeric(metadata[[col]]) & length(unique(metadata[[col]]))<=25){
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
        geom_point(aes_string(color=col), size=size, stroke=0, alpha=alpha) + 
        coord_fixed() +
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
        geom_point(aes_string(color=col), size=size, stroke=0, alpha=alpha) + 
        coord_fixed() +
        scale_color_gradient2(midpoint=0, low="royalblue4", mid="grey80", high="firebrick2", space ="Lab") +
        xlab(axes[1]) +    
        ylab(axes[2]) +
        ggtitle(col) +
        theme_linedraw() + 
        theme(plot.title = element_text(size = 10), legend.title = element_blank())
        
    }
    
    scatter_plots[[col]] <- tmp_plot
}

scatter_plot_panel <- wrap_plots(scatter_plots, ncol = n_col)#, widths=4, heights=4)

### save plot
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# scatter_plot_panel


ggsave(basename(plot_path),
       plot = scatter_plot_panel,
       device = 'png',
       path = result_dir,
       scale = 1,
       dpi = 300,
       width = width_panel,
       height = height_panel,
       limitsize = FALSE,
      )
