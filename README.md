[![DOI](https://zenodo.org/badge/475465311.svg)](https://zenodo.org/badge/latestdoi/475465311)

# Unsupervised Analysis Workflow
A general purpose [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to perform unsupervised analyses (dimensionality reduction and cluster analysis) and visualizations of high-dimensional data.

This workflow adheres to the module specifications of [MR.PARETO](https://github.com/epigen/mr.pareto), an effort to augment research by modularizing (biomedical) data science. For more details and modules check out the project's repository. Please consider **starring** and sharing modules that are interesting or useful to you, this helps others to find and benefit from the effort and me to prioritize my efforts!

**If you use this workflow in a publication, please don't forget to give credit to the authors by citing it using this DOI [10.5281/zenodo.8405360](https://doi.org/10.5281/zenodo.8405360).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [scRNA-seq Analysis](#single-cell-RNA-sequencing-(scRNA-seq)-data-analysis)
  * [Links](#links)
  * [Resources](#resources)
  * [Publications](#publications)


# Authors
- [Stephan Reichl](https://github.com/sreichl)
- [Raphael Bednarsky](https://github.com/bednarsky)
- [Christoph Bock](https://github.com/chrbock)

# Software
This project wouldn't be possible without the following software and their dependencies

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| clusterCrit    | https://CRAN.R-project.org/package=clusterCrit    |
| clustree       | https://doi.org/10.1093/gigascience/giy083        |
| ComplexHeatmap | https://doi.org/10.1093/bioinformatics/btw313     |
| densMAP        | https://doi.org/10.1038/s41587-020-00801-7        |
| ggally         | https://CRAN.R-project.org/package=GGally         |
| ggplot2        | https://ggplot2.tidyverse.org/                    |
| ggrepel        | https://CRAN.R-project.org/package=ggrepel        |
| igraph         | https://doi.org/10.5281/zenodo.3630268            |
| leidenalg      | https://doi.org/10.5281/zenodo.1469356            |
| pandas         | https://doi.org/10.5281/zenodo.3509134            |
| patchwork      | https://CRAN.R-project.org/package=patchwork      |
| PCA            | https://doi.org/10.1080/14786440109462720         |
| plotly express | https://plot.ly                                   |
| pymcdm         | https://doi.org/10.1016/j.softx.2023.101368       |
| scikit-learn   | http://jmlr.org/papers/v12/pedregosa11a.html      |
| scipy          |  https://doi.org/10.1038/s41592-019-0686-2        |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |
| umap-learn     | https://doi.org/10.21105/joss.00861               |


# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

The outlined analyses were performed using the programming languages R (ver) [ref] and Python (ver) [ref] unless stated otherwise. We applied both linear and non-linear unsupervised analysis methods for dimensionality reduction on normalized data for downstream analyses (e.g., clustering) and to visualize emerging patterns in lower dimensional embeddings.

**Dimensionality Reduction**

**Principal Component Analysis (PCA)**
We used Principal Component Analysis (PCA) [ref] from scikit-learn (ver) [ref] as the linear approach. We visualized [n_components] principal components and kept [X/all] components for downstream analyses. For diagnostic purposes we visualized the variance explained of all and the top 10% of principal components (PCs) using elbow- and cumulative-variance-plots, sequential pair-wise PCs for up to 10 PCs using scatter-, and density-plots (colored by [metadata_of_interest]), and finally loadings plots showing the magnitude and direction of the 10 most influential features for each PC as lollipop plot and biplot for sequential combinations of PCs. The R packages ggally (ver) [ref] and ggrepel (ver) [ref] were used to improve the diagnostic visualizations.

**Uniform Manifold Approximation and Projection (UMAP)**
Uniform Manifold Approximation projection (UMAP) from umap-learn (ver) [ref] was used as the non-linear approach. The metric [metric] and number of neighbors [n_neighbors] were used for the generation of a shared k-nearest-neighbor graph. The graph was embedded in [n_components] dimensions with minimum distance parameter [min_dist].
(Optional) We used the density preserving regularization option, densMAP [ref], during the embedding step, with default parameters to account for varying local density of the data within its original high dimensional space.

**Hierarchically Clustered Heatmap**
Hierarchically clustered heatmaps of scaled data (z-score) were generated using the R package ComplexHeatmap (ver) [ref]. The distance metric [metric] and clustering method [clustering_method] were used to determine the hierarchical clustering of observations (rows) and features (columns), respectively. The heatmap was annotated with metadata [metadata_of_interest]. The values were colored by the top percentiles (0.01/0.99) of the data to avoid shifts in the coloring scheme caused by outliers.

**Visualization**
The R-packages ggplot2 (ver) [ref] and patchwork (ver) [ref] were used to generate all 2D visualizations colored by metadata [metadata], feature(s) [features_to_plot], and/or clustering results.
Interactive visualizations in self-contained HTML files of all 2D and 3D projections/embeddings were generated using plotly express (ver) [ref].

**Cluster Analysis**

**Leiden Clustering**
We applied the Leiden algorithm (ver) [ref] to the UMAP KNN graphs specified by the respective parameters (metric, n_neighbors). The adjacency matrix of the KNN graph was converted to a weighted undirected graph using igraph (ver) [ref]. The Leiden algorithm was then applied to this graph, using the specified partition type [partition_types], resolution [resolutions], and number of iterations [n_iterations]. All clustering results were visualized as described above as 2D and interactive 2D and 3D plots for all available embedings/projections.

**Clustification Approach**
We developed/employed an iterative clustering approach, termed Clustification, that merges clusters based on misclassification. The method was initialized with the clustering result that had the highest resolution (i.e., the most clusters). We then performed iterative classification using the cluster labels, to determine if the classifier can distinguish between clusters or if they should be merged. This involved a stratified 5-fold cross-validation and a Random Forest classifier with default parameters (e.g., 100 trees). The predicted labels were retained for each iteration. Clusters were merged based on a normalized confusion matrix built using the predicted labels. This matrix was made symmetric and upper triangular, resulting in a similarity graph, such that each edge weight ranges from 0 to 1, where 0 means that the classifier was able to distinguish all observations between the two respective clusters. The stopping criterion was set such that if the maximum edge weight was less than 2.5% (i.e., 0.025 – less than 5% of observations are misclassified between any two clusters), the process would stop and return the current cluster labels. Otherwise, the two clusters connected by the maximum edge weight were merged. This process was repeated until the stopping criterion was met.

**Clustree Analysis & Visualization**
We performed cluster analysis and visualization using the Clustree package (ver) [ref] with the parameters [count_filter], [prop_filter], and [layout]. The default analysis produced a standard Clustree visualization, ordered by the number of clusters and annotated accordingly. For the custom analysis, we extended the default behaviour by adding [metadata_of_interest] as additional "clusterings". Metadata and features, specified in the configuration, were highlighted on top of the clusterings using aggregation functions. For numerical data, we used the [numerical_aggregation_option] function , and for categorical data, we used the [categorical_label_option] function.

**Cluster Validation - External Indices**
We validated/analyzed the clustering results by comparing them with all categorical metadata using external cluster indices. The complementary indices used were Adjusted Mutual Information (AMI) [ref], Adjusted Rand Index (ARI) [ref], Fowlkes-Mallows Index (FMI) [ref], Homogeneity, Completeness, and V-Measure [ref] from scikit-learn (ver) [ref]. The indices were calculated for each clustering result and each categorical metadata, and visualized using hierarchically clustered heatmaps.

**Cluster Validation - Internal Indices & MCDM using TOPSIS**
We performed internal cluster validation using six complementary indices: Silhouette, Calinski-Harabasz, C-index, Dunn index, Davis-Bouldin Score from the clusterCrit package (ver) [ref], and a weighted Bayesian Information Criterion (BIC) approach as described in [Reichl 2018 - Chapter 4.2.2 - Internal Indices](https://repositum.tuwien.at/handle/20.500.12708/3488). Due to computational cost, PCA results representing 90% of variance explained were used as input, and only a random sample proportion of [sample_proportion] was used. These internal cluster indices are linear, using Euclidean distance metrics. To rank all clustering results and [metadata_of_interest] from best to worst, we applied the Multiple-criteria decision-making (MCDM) method TOPSIS from the the Python package pymcdm (ver) [ref] to the internal cluster indices, as described in [Reichl 2018 - Chapter 4.3.1 - The Favorite Approach](https://repositum.tuwien.at/handle/20.500.12708/3488).

**The analysis and visualizations described here were performed using a publicly available Snakemake [ver] (ref) workflow [10.5281/zenodo.8405360](https://doi.org/10.5281/zenodo.8405360).**


# Features
The workflow perfroms the following analyses on each dataset provided in the annotation file. A result folder "unsupervised_analysis" is generated containing a folder for each dataset.

## Dimensionality Reduction
> _"High-dimensional spaces are where intuition goes to die and dimensionality reduction becomes the antidote to the curse of dimensionality."_ from Anonymous
- Principal Component Anlaysis (PCA) keeping all components (.pickle and .CSV)
  - diagnostics (.PNG):
      - variance: scree-plot and cumulative explained variance-plot of all and top 10% principal components
      - pairs: sequential pair-wise PCs for up to 10 PCs using scatter- and density-plots colored by [metadata_of_interest]
      - loadings: showing the magnitude and direction of the 10 most influential features for each Principal Component combination (Biplot, but without the data)
      - loadings lolliplot: showing the magnitude of the 10 most influential features for each Principal Component
- Uniform Manifold Approximation & Projection (UMAP)
  - k-nearest-neighbor graph (.pickle): generated using the [n_neighbors] parameter together with the provided [metrics].
    - fix any pickle load issue by specifying Python version to 3.9 (in case you want to use the graph downstream)
  - low dimensional embedding (.pickle and .CSV): using the precomputed-knn graph from before, embeddings are parametrized using [min_dist] and [n_components]
  - densMAP (optional): local density preserving regularization as additional dimensionality reduction method (i.e., all UMAP parameter combinations and downstream visualizations apply)
  - diagnostics (.PNG): 2D embedding colored by PCA coordinates, vector quantization coordinates, approximated local dimension, neighborhood Jaccard index
  - connectivity (.PNG): graph/network-connectivity plot with edge-bundling (hammer algorithm variant)
- Hierarchically Clustered Heatmap (.PNG)
    - hierarchically clustered heatmaps of scaled data (z-score) with configured distance [metrics] and clustering methods ([hclust_methods]). All combinations are computed, and annotated with [metadata_of_interest].
- Visualization
  -  2D metadata and feature plots (.PNG) of the first 2 principal components and all 2D embeddings, respectively.
  -  interactive 2D and 3D visualizations as self contained HTML files of all projections/embeddings.
- Results directories for each dataset have the following structure:
  -  "method" (containing all the data as .pickle and/or .CSV files)
    -  plots (for all visualizations as .PNG files)

## Cluster Analysis
> _"The validation of clustering structures is the most difficult and frustrating part of cluster analysis. Without a strong effort in this direction, cluster analysis will remain a black art accessible only to those true believers who have experience and great courage."_ from _Algorithms for Clustering Data (1988)_ by Jain & Dubes

- Clustering
    - Leiden algorithm
        - Applied to the UMAP KNN graphs specified by the respective parameters (metric, n_neighbors).
        - All algorithm specific parameters are supported: [partition_types], [resolutions], and [n_iterations].
    - Clustification: an ML-based clustering approach that iteratively merges clusters based on misclassification
        0. User: Specify a clustering method [method].
        1. Chose the clustering with the most clusters as starting point (i.e., overclustered).
        2. Iterative classification using the cluster labels, to determine if the classifier can distinguish between clusters or if they should be merged.
            - Stratified 5-fold CV
            - RF with 100 trees (i.e., defaults)
            - Retain predicted labels
        3. Merging of clusters.
            - Build a normalized confusion matrix using the predicted labels.
            - Make it symmetric and upper triangle, resulting in a similarity graph.
            - Edge weight ranges from 0 to 1, where 0 means that the classifier was able to distinguish all observations between the two respective clusters.
            - Check stopping criterion: if maximum edge weight < 2.5% (i.e., 0.025 – less than 5% of observations are misclassified between any two clusters).
              - -> STOP and return current cluster labels
            - Merge the two clusters connected by the maximum edge weight.
        4. Back to 2. using the new labels.
- Clustree analysis and visualization
    - The following clustree specific parameters are supported: [count_filter], [prop_filter], and [layout].
    - default: produces the standard clustree visualization, ordered by number of clusters and annotated.
    - custom: extends default by adding [metadata_of_interest] as additional "clusterings".
    - metadata and features, specified in the config, are highlighted on top of the clusterings using aggregation functions
        - numeric: available aggregation functions: mean, median, max, min
        - categorical: available aggregation functions: "pure" or "majority"
- Cluster Validation
    - External cluster indices are determined comparing all clustering results with all categorical metadata
        - all complementary indices from sklearn are used: [AMI](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html#sklearn.metrics.adjusted_mutual_info_score), [ARI](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html#sklearn.metrics.adjusted_rand_score), [FMI](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.fowlkes_mallows_score.html#sklearn.metrics.fowlkes_mallows_score), [**Homogeneity** and **Completeness** and **V**-Measure](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.homogeneity_completeness_v_measure.html#sklearn.metrics.homogeneity_completeness_v_measure)
    - Internal cluster indices are determined for each clustering and [metadata_of_interest]
        - 6 complementary indices are used
            - 5 from the package [clusterCrit](https://rdrr.io/cran/clusterCrit/man/intCriteria.html): Silhouette, Calinski-Harabasz, C-index, Dunn index, Davis-Bouldin Score.
            - 1 weighted Bayesian Information Criterion (BIC) approach, previously described in [Reichl 2018 - Chapter 4.2.2 - Internal Indices](https://repositum.tuwien.at/handle/20.500.12708/3488) 
        - Due to the comutational cost PCA results, representing 90% of variance explained, are used for as input and a [sample_proportion] can be configured.
        - Caveat: internal cluster indices are linear i.e., using Euclidean distance metrics.
    - Multiple-criteria decision-making (MCDM) using TOPSIS for ranking clustering results
        - The MCDM method TOPSIS is applied to the internal cluster indices to rank all clustering results (and [metadata_of_interest]) from best to worst.
        - This approach has been described in [Reichl 2018 - Chapter 4.3.1 - The Favorite Approach](https://repositum.tuwien.at/handle/20.500.12708/3488)
        - Caveat: Silhouette score sometimes generates NA due to a known [bug](https://github.com/cran/clusterCrit/pull/1/commits/b37a5e361d0a12f9d3900089aa03e3947d0d4ef7). Clusterings with NA scores are removed before TOPSIS is applied.
- Visualization
    - all clustering results as 2D and interactive 2D & 3D plots for all available embedings/projections.
    - external cluster indices as hierarchically clustered heatmaps, aggregated in one panel.
    - internal cluster indices as one heatmap with clusterings (and [metadata_of_interest]) sorted by TOPSIS ranking from top to bottom and split cluster indices split by type (cost/benefit functions to be minimized/maximized).


# Usage
Here are some tips for the usage of this workflow:
- Start with minimal parameter combinations and without UMAP diagnostics and connectivity plots (they are computational expensive and slow).
- Heatmaps require **a lot** of memory, hence the memory allocation is solved dynamically based on retries. If a out-of-memory exception occurs the flag `--retries X` can be used to trigger automatic resubmission X time upon failure with X times the memory.
- Clustification performance scales with available cores, i.e., more cores faster internal parallelization of Random Forest training & testing.
- Cluster indices are extremely compute intense and scale linearly with every additional clustering result and specified metadata (can be skipped).

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
We provide a minimal example of the analysis of the [UCI ML hand-written digits datasets](https://archive.ics.uci.edu/ml/datasets/Optical+Recognition+of+Handwritten+Digits) imported from [sklearn](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_digits.html) in the [test folder](./test/):
- config
    - configuration: config/config.yaml
    - sample annotation: digits_unsupervised_analysis_annotation.csv
- data
    - dataset (1797 observations, 64 features): digits_data.csv
    - metadata (consisting of the ground truth label "target"): digits_labels.csv
- results will be generated in the configured subfolder `./test/results/`
- performance: on an HPC it took less than 5 minutes to complete a full run (with up to 32GB of memory per task)

# single-cell RNA sequencing (scRNA-seq) data analysis
Unsupervised analyses, dimensionality reduction and cluster analysis, are corner stones of scRNA-seq data analyses.
A full run on a [published](https://www.nature.com/articles/s41588-020-0636-z) scRNA-seq [cancer dataset](https://www.weizmann.ac.il/sites/3CA/colorectal) with 21,657 cells and 18,245 genes took 2.5 hours to complete (without heatmaps, with 32GB memory and 8 cores for clustification).
Below are configurations of the two most commonly used frameworks, [scanpy](https://scanpy.readthedocs.io/en/stable/index.html) (Python) and [Seurat](https://satijalab.org/seurat/) (R), and the original package's defaults as comparison and to facilitate reproducibility:

UMAP for dimensionality reduction
- [umap-learn](https://umap-learn.readthedocs.io/en/latest/api.html)
    - initialization: spectral
    - metric: Euclidean
    - neighbors: 15
    - min. distance: 0.1
- [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html#scanpy.pp.neighbors)
    - initialization: spectral
    - metric: Euclidean
    - neighbors: 15
    - min. distance: **0.5**
- [Seurat](https://satijalab.org/seurat/reference/runumap)
    - initialization: **PCA**
    - method: "uwot" (not umap-learn package)
    - metric: **Cosine (or Correlation)**
    - neighbors: **30**
    - min. distance: **0.3**

Leiden algorithm for clustering
- [leidenalg](https://leidenalg.readthedocs.io/en/stable/reference.html)
    - no defaults
- [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)
    - input: batch balanced UMAP KNN graph
    - partition type: RBConfigurationVertexPartition
    - resolution: 1
- [Seurat](https://github.com/satijalab/seurat/blob/763259d05991d40721dee99c9919ec6d4491d15e/R/clustering.R#L344)
    - input: SNN graph
    - partition type: RBConfigurationVertexPartition
    - resolution: 0.8
  

# Links
- [GitHub Repository](https://github.com/epigen/unsupervised_analysis/)
- [GitHub Page](https://epigen.github.io/unsupervised_analysis/)
- [Zenodo Repository](https://doi.org/10.5281/zenodo.8405360)
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/unsupervised_analysis)

# Resources
- Recommended compatible [MR.PARETO](https://github.com/epigen/mr.pareto) modules
  - for upstream processing:
    - [ATAC-seq Processing](https://github.com/epigen/atacseq_pipeline) to quantify  chromatin accessibility.
    - [scRNA-seq Data Processing & Visualization](https://github.com/epigen/scrnaseq_processing_seurat) for processing and preparing single cell data as input.
    - [Split, Filter, Normalize and Integrate Sequencing Data](https://github.com/epigen/spilterlize_integrate) process and preapre sequencing data as input.
    - [Perturbation Analysis using Mixscape from Seurat](https://github.com/epigen/mixscape_seurat) to identify perturbed cells from pooled (multimodal) CRISPR screens with sc/snRNA-seq read-out (scCRISPR-seq) as input.
- [Reichl, S. (2018). Mathematical methods in single cell RNA sequencing analysis with an emphasis on the validation of clustering results [Diploma Thesis, Technische Universität Wien]](https://doi.org/10.34726/hss.2018.49662)

# Publications
The following publications successfully used this module for their analyses.
- ...