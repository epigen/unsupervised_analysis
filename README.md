# Unsupervised Analysis Worfklow
A general purpose [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to perform selected unsupervised analyses and visualizations of high dimensional data.

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

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

# Authors
- [Stephan Reichl](https://github.com/sreichl)

# Software
This project wouldn't be possible without the following software and their dependencies

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| ggplot2        | https://ggplot2.tidyverse.org/                    |
| pandas         | https://doi.org/10.5281/zenodo.3509134            |
| patchwork      | https://CRAN.R-project.org/package=patchwork      |
| scikit-learn   | http://jmlr.org/papers/v12/pedregosa11a.html      |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |
| umap-learn     | https://doi.org/10.21105/joss.00861               |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

--- COMING SOON ---

# Features
The workflow perfroms the following analyses on each dataset provided in the annotation file. A result folder "unsupervised_analysis" is generated containing a folder for each dataset.
- Principal Component Anlaysis (PCA) keeping all components
  - diagnostics: scree-plot and cumulative explained variance-plot of all and top 10% principal components
- Uniform Manifold Approximation & Projection (UMAP)
  - graph: generated using the maximum n_neighorhood parameter together with the provided metrics
  - embedding: using the precomputed-knn graph from before, embeddings are parametrized using min_dist and n_components
  - diagnostics: 2D embedding colored by PCA coordinates, vector quantization coordinates, approximated local dimension, neighborhood Jaccard index
  - connectivity: graph/network-connectivity plot with edge-bundling (hammer algorithm variant)
- Visualization
  -  2D metadata plots of the first 2 principal components and all 2D embeddings, depending on the method
- Results directories for each dataset have the following structure:
  -  "method" (containing all the data as .pickle or .CSV files)
    -  plots (for all visualizations)


# Usage
Here are some tips for the usage of this workflow:
- Start with minimal parameter combinations and without UMAP diagnostics and connectivity plots (computational expensive and slow)
- 

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
--- COMING SOON ---
