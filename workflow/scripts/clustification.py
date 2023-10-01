#### load a clustering result with high-resolution (i.e., overclustered) and iteratively merge clusters using a RF classifier ####

#### libraries
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import confusion_matrix
# from sklearn.metrics import f1_score, accuracy_score

# clustification function
def iterative_classification(data, labels, n_trees=100, max_iterations=100):

    for iteration in range(max_iterations):
        print(f"\nIteration {iteration+1}")
        # Train & test 5 times
        new_labels = np.zeros_like(labels)
        new_data = np.zeros_like(data)
        skf = StratifiedKFold(n_splits=5, shuffle=True)
        
        for train_index, test_index in skf.split(data, labels):
            X_train, X_test = data[train_index], data[test_index]
            y_train, y_test = labels[train_index], labels[test_index]
            clf = RandomForestClassifier(n_estimators=n_trees,
                                        random_state = 42,
                                        n_jobs = -1)
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            new_labels[test_index] = y_pred
            new_data[test_index] = X_test

        # determine symmetric confusiton matrix interpreted as weighted graph
        cm = confusion_matrix(y_true=labels, y_pred=new_labels, normalize='true')
        cm_symmetric = (cm + cm.T) / 2
        cm_symmetric = np.triu(cm_symmetric, 1)
        max_weight = np.max(cm_symmetric)
        print("max edge weight: {}".format(max_weight))
        
        # check stopping criterion
        if max_weight < 0.025:
            print("Misclassification threshold met, stopping.")
            break
        
        merge_indices = np.unravel_index(np.argmax(cm_symmetric, axis=None), cm_symmetric.shape)
        labels = np.array([merge_indices[0] if num == merge_indices[1] else num for num in labels])
        labels, _ = pd.factorize(labels)
        data = new_data
        print("#clusters: {}".format(len(set(labels))))
    return labels


#### configurations

# inputs
data_path = os.path.join(snakemake.input[0])
clusterings_path = os.path.join(snakemake.input[1])

# parameters
samples_by_features = int(snakemake.params['samples_by_features']) #0

# outputs
result_path = os.path.join(snakemake.output["clustering"])

# load data
# check data orientation to fit: samples/observations x features
if samples_by_features == 1:
    data = pd.read_csv(data_path, index_col=0)
else:
    data = pd.read_csv(data_path, index_col=0).T

    
clusterings = pd.read_csv(clusterings_path, index_col=0)
clustering_init = clusterings[clusterings.nunique().idxmax()].to_numpy().ravel()

# run clustification
clustering_new = iterative_classification(data.to_numpy(), clustering_init)#, n_trees=100, max_iterations=100)

# save clustering as CSV
pd.DataFrame({"clustification_clustering": clustering_new}, index=data.index).to_csv(result_path, index=True)
