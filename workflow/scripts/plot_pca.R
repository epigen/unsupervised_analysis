#### load libraries
library("ggplot2")
library("patchwork")
library("GGally")
library("ggrepel")

### configurations

# inputs
data_path <- snakemake@input[["dimred_data"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_data.csv"
var_path <- snakemake@input[["dimred_var"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_var.csv"
axes_path <- snakemake@input[["dimred_axes"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_axes.csv"
loadings_path <- snakemake@input[["dimred_loadings"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_loadings.csv"
metadata_path <- snakemake@input[["metadata"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/data/digits_labels.csv"

diagnostics_path <- snakemake@output[["diagnostics_plot"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/plots/PCA_diagnostics.png"
pairs_path <- snakemake@output[["pairs_plot"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/plots/PCA_pairs.png"
loadingsplot_path <- snakemake@output[["loadings_plot"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/plots/PCA_loadings.png"

pairs_size <- snakemake@config[["scatterplot2d"]][["size"]]/10 # 0.5
pairs_alpha <- snakemake@config[["scatterplot2d"]][["alpha"]]/2 # 1
metadata_col <- snakemake@config[["metadata_of_interest"]] # "target"

result_dir <- file.path(dirname(diagnostics_path))
# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

### load data
data <- read.csv(file=file.path(data_path), row.names=1, header=TRUE)
metadata <- read.csv(file=file.path(metadata_path), row.names=1, header=TRUE)

# prepare metadata
if(metadata_col==""){
    metadata_col <- colnames(metadata)[1]
}

data_axes <- read.csv(file=file.path(axes_path), row.names=1, header=TRUE)
data_loadings <- read.csv(file=file.path(loadings_path), row.names=1, header=TRUE)

data_var <- read.csv(file=file.path(var_path), row.names=1, header=TRUE)
colnames(data_var) <- c('var')
data_var$PC <- as.numeric(rownames(data_var))+1

### variance plot
print("Variance plots")
# plot specifications
n_col <- 2
width <- 5
height <- 3
point_size <- 0.25
line_size <- 0.1

width_panel <- n_col * width
height_panel <- 2 * height

# number of top 10% of PCs
top_n <- ceiling(nrow(data_var)*0.1)

# make plots
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)

pca_plots <- list()

pca_plots[["scree_all"]] <- ggplot(data_var, aes(x=PC,y=var, group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot of all Principal Components")+
theme_linedraw()+ theme(plot.title = element_text(size = 10))

pca_plots[["cum_all"]] <- ggplot(data_var, aes(x=PC,y=cumsum(var), group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Cumulative Explained Variance of all Principal Components")+
theme_linedraw()+ theme(plot.title = element_text(size = 10))

pca_plots[["scree_top"]] <- ggplot(data_var[1:top_n,], aes(x=PC,y=var, group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle(paste0("Scree Plot of top ", top_n," Principal Components"))+
theme_linedraw()+ theme(plot.title = element_text(size = 10))


pca_plots[["cum_top"]] <- ggplot(data_var[1:top_n,], aes(x=PC,y=cumsum(var), group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle(paste0("Cumulative Explained Variance of top ", top_n," Principal Components"))+
theme_linedraw()+ theme(plot.title = element_text(size = 10))


pca_plot_panel <- wrap_plots(pca_plots, ncol = n_col)

# save diangostics plot
# pca_plot_panel

ggsave(basename(diagnostics_path),
       plot = pca_plot_panel,
       device = 'png',
       path = result_dir,
       scale = 1,
       dpi = 300,
       width = width_panel,
       height = height_panel,
       limitsize = FALSE,
      )

### pairs plot
print("Pairs plot")
n_dim <- min(10, ncol(data))

# convert to categorical if less than 25 unique integer values
if (is.numeric(metadata[[metadata_col]]) & length(unique(metadata[[metadata_col]]))<=25){
    if(all(metadata[[metadata_col]] == round(metadata[[metadata_col]]))){
        metadata[metadata_col] <- as.factor(metadata[[metadata_col]])
    }
}
# if a metadata class is empty ("") fill with "unknown"
if (any(metadata[[metadata_col]]=="")){
    metadata[metadata[[metadata_col]]=="", metadata_col] <- "unknown"
}

# legend parameter according to data type
if (is.numeric(metadata[[metadata_col]])){
    legend <- NULL
}else{
    legend <- 1
}

# options(repr.plot.width=10, repr.plot.height=10)

# make pairs plot
pairs_plot <- ggpairs(
  data = data,
  mapping = ggplot2::aes(color = metadata[[metadata_col]]),
  columns = 1:n_dim,
  title = paste0("PCA pairs plot colored by ",metadata_col),
  upper = list(continuous = "density"),
  lower = list(continuous = wrap("points", alpha = pairs_alpha, size = pairs_size)),
  diag = list(continuous = "densityDiag"),
  params = NULL,
  xlab = NULL,
  ylab = NULL,
  axisLabels = c("show", "internal", "none"),
  columnLabels = data_axes[1:n_dim,'label'],
  labeller = "label_value",
  switch = NULL,
  showStrips = NULL,
  legend = legend,
  cardinality_threshold = 15,
  progress = NULL,
  proportions = NULL
)+  
theme(legend.position = "bottom") + 
labs(fill = metadata_col)

if (is.numeric(metadata[[metadata_col]])){
    pairs_plot <- pairs_plot + scale_color_gradient2(midpoint=0, low="royalblue4", mid="grey80", high="firebrick2", space ="Lab")
}

# save pairs plot
ggsave(basename(pairs_path),
       plot = pairs_plot,
       device = 'png',
       path = result_dir,
       scale = 1,
       dpi = 300,
       width = n_dim,
       height = n_dim,
       limitsize = FALSE,
      )


### loadings plot
print("Loadings plot")

# plot specifications
n_col <- min(5,n_dim)
height_panel <- ceiling(n_dim/n_col)*4
width_panel <- n_col * 4

loading_plots <- list()

for(i in 1:n_dim){
    tmp_x <- i
    tmp_y <- i+1
    
    # determine top 10 features per PC combination
    top_features <- rownames(data_loadings)[order(-(data_loadings[paste0("PC_",tmp_x)]^2 + data_loadings[paste0("PC_",tmp_y)]^2))][1:10]
    
    # subset data
    tmp_loadings <- data_loadings[top_features, c(paste0("PC_",tmp_x), paste0("PC_",tmp_y))]
    tmp_loadings$features <- rownames(tmp_loadings)
    text_var <- "features"
    
    # plot data
    loading_plots[[i]] <- ggplot(data=tmp_loadings, aes_string(x=paste0("PC_",tmp_x), y=paste0("PC_",tmp_y), label=text_var))+
    geom_segment(data=tmp_loadings, aes_string(x=0, y=0, xend=paste0("PC_",tmp_x), yend=paste0("PC_",tmp_y)), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_label_repel()+
    xlab(paste0("Principal Component ",tmp_x)) +
    ylab(paste0("Principal Component ",tmp_y)) +
    theme_linedraw()
    
#     print(tmp_plot)
}

loadings_plot_panel <- wrap_plots(loading_plots, ncol = n_col)

# save diangostics plot
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# loadings_plot_panel

ggsave(basename(loadingsplot_path),
       plot = loadings_plot_panel,
       device = 'png',
       path = result_dir,
       scale = 1,
       dpi = 300,
       width = width_panel,
       height = height_panel,
       limitsize = FALSE,
      )
