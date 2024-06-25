# Configuration

You need one configuration file to configure the analyses and one annotation file describing the data to run the complete workflow. If in doubt read the comments in the config and/or try the default values. We provide a full example including data, configuration, results an report in `test/` as a starting point.

- project configuration (`config/config.yaml`): Different for every project and configures the analyses to be performed.
- sample annotation (annotation): CSV file consisting of four mandatory columns.
    -  name: A unique name of the dataset (tip: keep it short but descriptive).
    -  data: Path to the tabular data as comma separated table (CSV).
    -  metadata: Path to the metadata as comma separated table (CSV) with the first column being the index/identifier of each sample/observation and every other column metadata for the respective sample (either numeric or categorical, not mixed). **No NaN or empty values allowed.**
    -  samples_by_features: Boolean indicator if the data matrix is samples (rows) x features (columns) -> (0==no, 1==yes).
