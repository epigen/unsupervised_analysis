# Configuration

You need one configuration file for the analyses and one annotation file for the data to run the complete workflow. Always use absolute paths. If in doubt read the comments in the config and/or try the default values. We provide a full example including data, configuration, reuslts an report in .test/ as a starting point.

- project configuration (config/config.yaml): different for every project and configures the analyses to be performed
- sample annotation (sample_annotation): CSV file consisting of three columns
    -  name: name of the data set (tip: keep it short)
    -  data: absolute path to the tabular data as CSV
    -  metadata: absolute path to the metadata as CSV with the first column being the index/identifier of each sample/observation and every other coloumn metadata for the respective sample
    -  samples_by_features: 0 or 1 as boolean indicator if data matrix is samples (rows) x features (columns) -> (0==no, 1==yes)
