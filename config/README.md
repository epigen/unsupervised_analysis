# Configuration

You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. Always use absolute paths. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project and configures the analyses to be performed
- sample annotation (sample_annotation): CSV file consisting of three columns
    -  name: name of the data set (tip: keep it short)
    -  data: absolute path to the tabular data as CSV
    -  metadata: abosulte path to the metadata as CSV with the first column being the index/identifier of each sample and every other coloumn metadata for the respective sample
