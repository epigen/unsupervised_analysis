#### load all the internal indices of all clusterings and metadata and apply the MCDM/A method TOPSIS to rank them ####

#### libraries
import os
import numpy as np
import pandas as pd
from pymcdm.methods import TOPSIS
from pymcdm.helpers import rrankdata

#### configurations

# inputs
indices_path = os.path.join(snakemake.input["internal_indices"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/internal_indices.csv"

# outputs
indices_ranked_path = os.path.join(snakemake.output["internal_indices_ranked"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/internal_indices_ranked.csv"

# load the internal indices
indices = pd.read_csv(indices_path, index_col=0)

# perform TOPSIS

# give all scores equal weights as they are supposed to be complementary
weights = np.full(indices.shape[1], 1.0/indices.shape[1])
# set 1 for benefit (max) and -1 for cost (min) functions (hard coded for internal indices in order)
types = np.array([1, 1, -1, -1, 1, -1, -1])
# create object with defaults
topsis = TOPSIS()
# run TOPSIS
pref = topsis(indices.to_numpy(), weights, types)
# extract ranks
rank = rrankdata(pref)
print(rank)

indices_ranked = indices.iloc[pref.argsort()[::-1]]

# save ranked indices
indices_ranked.to_csv(indices_ranked_path)
