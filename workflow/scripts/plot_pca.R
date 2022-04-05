#### load libraries
library("ggplot2")
library("patchwork")

### configurations

# inputs
var_path <- snakemake@input[["var_data"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/PCA/PCA_var.csv"

plot_path <- snakemake@output[["diagnostics_plot"]] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/PCA/plots/PCA_diagnostics.png"

result_dir <- file.path(dirname(plot_path))
# make result directory if not exist
if (!dir.exists(result_dir)){
    dir.create(result_dir, recursive = TRUE)
}

### load data
data <- read.csv(file=file.path(var_path), row.names=1, header=TRUE)
colnames(data) <- c('var')
data$PC <- as.numeric(rownames(data))+1

# plot specifications
n_col <- 2
width <- 5
height <- 3
point_size <- 0.25
line_size <- 0.1

width_panel <- n_col * width
height_panel <- 2 * height

# number of top 10% of PCs
top_n <- ceiling(nrow(data)*0.1)

### make plots
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)

pca_plots <- list()

pca_plots[["scree_all"]] <- ggplot(data, aes(x=PC,y=var, group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot of all Principal Components")+
theme_linedraw()+ theme(plot.title = element_text(size = 10))

pca_plots[["cum_all"]] <- ggplot(data, aes(x=PC,y=cumsum(var), group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Cumulative Explained Variance of all Principal Components")+
theme_linedraw()+ theme(plot.title = element_text(size = 10))

pca_plots[["scree_top"]] <- ggplot(data[1:top_n,], aes(x=PC,y=var, group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle(paste0("Scree Plot of top ", top_n," Principal Components"))+
theme_linedraw()+ theme(plot.title = element_text(size = 10))


pca_plots[["cum_top"]] <- ggplot(data[1:top_n,], aes(x=PC,y=cumsum(var), group=1))+
  geom_point(size=point_size)+
  geom_line(size=line_size)+
    xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle(paste0("Cumulative Explained Variance of top ", top_n," Principal Components"))+
theme_linedraw()+ theme(plot.title = element_text(size = 10))


pca_plot_panel <- wrap_plots(pca_plots, ncol = n_col)

### save plot
# pca_plot_panel

ggsave(basename(plot_path),
       plot = pca_plot_panel,
       device = 'png',
       path = result_dir,
       scale = 1,
       dpi = 300,
       width = width_panel,
       height = height_panel,
       limitsize = FALSE,
      )

