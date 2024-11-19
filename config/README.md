# Configuration

You need one configuration file to configure the analyses and one annotation file describing the data to run the complete workflow. If in doubt read the comments in the config and/or try the default values. We provide a full example including data and configuration in `test/` as a starting point.

- project configuration (`config/config.yaml`): Different for every project and configures the analyses to be performed.
- sample annotation (annotation): CSV file consisting of four mandatory columns.
    -  name: A unique name for the dataset (tip: keep it short but descriptive).
    -  data: Path to the tabular data as a comma-separated table (CSV).
    -  metadata: Path to the metadata as a comma-separated table (CSV) with the first column being the index/identifier of each observation/sample and every other column metadata for the respective observation (either numeric or categorical, not mixed). **No NaN or empty values allowed, and no special characters (all except a-z, 0-9, `_`) in the index.**
    -  samples_by_features: Boolean indicator if the data matrix is observations/samples (rows) x features (columns): 0==no, 1==yes.

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.