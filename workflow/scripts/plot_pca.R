#### load libraries
library("ggplot2")
library("patchwork")
library("GGally")

### configurations

# inputs
data_path <- snakemake@input[["dimred_data"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_data.csv"
var_path <- snakemake@input[["dimred_var"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_var.csv"
axes_path <- snakemake@input[["dimred_axes"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/PCA_axes.csv"
metadata_path <- snakemake@input[["metadata"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/data/digits_labels.csv"

diagnostics_path <- snakemake@output[["diagnostics_plot"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/plots/PCA_diagnostics.png"
pairs_path <- snakemake@output[["pairs_plot"]] # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/PCA/plots/PCA_pairs.png"

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
data_axes <- read.csv(file=file.path(axes_path), row.names=1, header=TRUE)

data_var <- read.csv(file=file.path(var_path), row.names=1, header=TRUE)
colnames(data_var) <- c('var')
data_var$PC <- as.numeric(rownames(data_var))+1

### diagnostics plot
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
n_dim <- min(10, ncol(data))

# convert to categorical if less than 25 unique integer values
if (is.numeric(metadata[[metadata_col]]) & length(unique(metadata[[metadata_col]]))<=25){
    if(all(metadata[[metadata_col]] == round(metadata[[metadata_col]]))){
        metadata[metadata_col] <- as.factor(metadata[[metadata_col]])
    }
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
