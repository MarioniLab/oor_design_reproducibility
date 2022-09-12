### Run MELD ###
import graphtools as gt
import meld
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union
from numpy.typing import ArrayLike

def run_meld(X_red_dim: ArrayLike, 
            sample_labels: List[str], 
            conditions: List[str], 
            k: int = 15):
    '''
    Run MELD
    - X_red_dim: c x d matrix of dimensionality reduction to use for graph construction
    - sample_labels: assignment of cells to samples
    - conditions: vector of condition names
    '''
    # Make graph
    graph = gt.Graph(X_red_dim, knn=int(k))
    # Make MELD object
    meld_op = meld.MELD()
    meld_op.graph = graph
    # Compute density
    meld_fit = meld_op.transform(sample_labels=np.array(sample_labels))

    # Mean density per replicates
    mean_density = pd.DataFrame(
        np.zeros(shape=(meld_fit.shape[0], len(conditions))),
        index=meld_fit.index,
        columns=conditions,
    )

    for c in conditions:
        c_mean = meld_fit.loc[:, [c in x for x in meld_fit.columns]].mean(1)
        mean_density[c] = c_mean

    # From density to likelihood per condition
    likelihoods = meld.utils.normalize_densities(mean_density)
    likelihoods.columns = [col.split("_")[0] for col in likelihoods.columns]
    return(likelihoods)
