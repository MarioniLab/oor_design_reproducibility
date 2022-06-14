from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
import anndata
import pandas as pd
import numpy as np

from scipy.sparse import csr_matrix
from scipy import stats
from statsmodels.stats.multitest import multipletests


from sklearn.metrics import roc_curve, auc


def _check_inputs(adata: AnnData) -> bool:
    """
    Checks that all the required inputs are there
    - Cell-level confidence metric (stored in `adata.obs['confidence']`)
    - Common dimensions for controls and condition samples (stored in `adata.obsm['X_dims']`)
    - Assignment of cells to sample/donors (stored in `adata.obs['sample_id']`)
    - Assignment of cells to `adata.obsm['cell_nhoods']`, binary of dimensions cells x nhoods
    """

# --- Assign cells to populations of similar cells --- #


def make_cell_nhoods(adata: AnnData,
                     method: str,
                     #  X_dims: str,
                     neighbors_key: str = None,
                     clusters_key: str = None) -> None:
    """Assign cells to populations of similar cells in the phenotypic space defined by `X_dims`

   Params:
   ------
   - adata: AnnData object
    - X_dims: string of name of slot in `adata.obsm` storing the common phenotypic space to use for analysis
   - method: if `KNN`, compute milo neighbourhoods based on graph stored in `adata.obsp['neighbors_key']`, if `clusters` assign cells to clusters stored in `adata.obs['clusters_keys']`
   - neighbors_key: only used if `method = "KNN"`, which key in `adata.obsp` to use for assignment of cells to neighbourhoods
   - clusters_key: only used if `method = "clusters"`, which key in `adata.obsp` to use for assignment of cells to clusters

   Returns:
   -------
   None, adds in place `adata.obsm['cell_nhoods']`, a binary sparse matrix of dimensions cells x neighbourhoods
    """
    if method == 'clusters':
        _make_nhoods_clusters(adata, clusters_key=clusters_key)
    elif method == 'KNN':
        None
    else:
        raise ValueError("method must be either 'KNN' or 'clusters'")


def _make_nhoods_KNN(adata: AnnData, neighbors_key: str = None):
    """make milo neighbourhoods"""
    None


def _make_nhoods_clusters(adata: AnnData, clusters_key: str = None):
    """make neighbourhoods from cluster assignment"""
    try:
        clusters = adata.obs[clusters_key]
    except KeyError:
        raise ValueError('clusters_key is not a column in adata.obs')
    adata.obsm['cell_nhoods'] = csr_matrix(pd.get_dummies(clusters).values)


# --- Compute confidence metric per sample/donor --- #

def nhood_confidence(adata: AnnData,
                     confidence_col: str,
                     sample_col: str) -> None:
    """Aggregate cell-level confidence by cell neighourhood and sample
    Params:
    ------
    - adata: AnnData object
    - confidence_col: column in `adata.obs` storing cell-level confidence statistic. Higher values correspond to higher confidence in reference mapping.
    - sample_col: column in `adata.obs` storing assignment of cells to samples

    Returns:
    -------
    None, adds in place `adata.uns['nhood_adata']` of dimensions nhoods x samples storing cell counts in .X and confidence in `layers['confidence']` 
    """
    try:
        nhood_mat = adata.obsm["cell_nhoods"].copy()  #  cells x nhoods
    except:
        raise KeyError(
            "adata.obsm['cell_nhoods'] not found, please run make_cell_nhoods first")
    try:
        sample_dummies = pd.get_dummies(adata.obs[sample_col]).values
    except:
        raise KeyError("sample_col is not a column in adata.obs")

    # Aggregate confidence score by cell neighbourhood
    nhood_scores = np.array([]).reshape(nhood_mat.shape[1], 0)
    #  better ideas than a for loop?
    for i in np.arange(sample_dummies.shape[1]):
        ixs = sample_dummies[:, i] == 1
        nh_score = nhood_mat[ixs, :].T.dot(
            csr_matrix(adata.obs.loc[ixs, confidence_col]).T).toarray()
        nh_score = nh_score/np.array(nhood_mat[ixs, :].T.sum(1))
        nh_score[np.isnan(nh_score)] = 0
        nhood_scores = np.hstack([nhood_scores, nh_score])

    # Make anndata object for nhoods x samples if missing
    try:
        adata.uns['nhood_adata'].layers['confidence'] = nhood_scores
    except KeyError:
        _add_nhoods_adata(adata, sample_col)
        adata.uns['nhood_adata'].layers['confidence'] = nhood_scores


def _add_nhoods_adata(
    adata: AnnData,
    sample_col: str,
):
    '''
    - adata
    - sample_col: string, column in adata.obs that contains sample information 
    (what should be in the columns of the nhoodCount matrix)

    Returns: None
    Updated adata.uns slot to contain adata.uns["nhood_adata"], where:
    - adata.uns["nhood_adata"].obs_names are neighbourhoods
    - adata.uns["nhood_adata"].var_names are samples
    - adata.uns["nhood_adata"].X is the matrix counting the number of cells from each
    sample in each neighbourhood
    '''
    try:
        nhoods = adata.obsm["cell_nhoods"]
    except KeyError:
        raise KeyError(
            'Cannot find "cell_nhoods" slot in adata.obsm -- please run make_cell_nhoods'
        )
    #  Make nhood abundance matrix
    sample_dummies = pd.get_dummies(adata.obs[sample_col])
    all_samples = sample_dummies.columns
    sample_dummies = csr_matrix(sample_dummies.values)
    nhood_count_mat = adata.obsm["cell_nhoods"].T.dot(sample_dummies)
    nhood_var = pd.DataFrame(index=all_samples)
    nhood_adata = anndata.AnnData(X=nhood_count_mat, var=nhood_var)
    nhood_adata.uns["sample_col"] = sample_col
    # Save nhood index info if milo nhoods
    if "nhood_ixs_refined" in adata.obs.columns:
        nhood_adata.obs["index_cell"] = adata.obs_names[adata.obs["nhood_ixs_refined"] == 1]
        nhood_adata.obs["kth_distance"] = adata.obs.loc[adata.obs["nhood_ixs_refined"]
                                                        == 1, "nhood_kth_distance"].values
    adata.uns["nhood_adata"] = nhood_adata


# --- Make design matrix --- #


def make_design(adata: AnnData,
                categorical_covariates: List[str] = None,
                continuous_covariates: List[str] = None
                ) -> None:
    """Assign covariates used for testing to samples. Covariates need to be uniquely assigned to samples (i.e. only one value per sample).

    Params:
    ------
    - adata: AnnData object, storing nhoods x sample matrix in `adata.uns['nhood_adata']`
    - categorical_covariates: list of categorical covariates to assign to samples
    - categorical_covariates: list of continuous covariates to assign to samples

    Return:
    -------
    None, adds in place sample covariates to `adata.uns['nhood_adata'].var` and design matrix to `adata.uns['nhood_adata'].varm`
    """
    try:
        nhood_adata = adata.uns['nhood_adata'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'] not found -- please run nhood_confidence first")

    sample_col = nhood_adata.uns['sample_col']
    design_mat = pd.DataFrame()
    sample_obs = pd.DataFrame()
    if categorical_covariates is not None:
        #  Check categorical type
        covs_dtypes = adata.obs[categorical_covariates].dtypes.values
        assert not any([pd.api.types.is_numeric_dtype(x) for x in covs_dtypes]
                       ), 'some categorical_covariates are numeric. Add to continuous_covariates instead'
        # Get sample-level obs
        sample_obs_categorical = adata.obs[[sample_col] + categorical_covariates].drop_duplicates(
        ).reset_index(drop=True).set_index(sample_col)
        assert sample_obs_categorical.index.is_unique, 'Sample covariates are not unique.'
        # Make design matrix
        design_mat_categorical = pd.get_dummies(
            sample_obs_categorical[categorical_covariates])
        design_mat = pd.concat([design_mat, design_mat_categorical], axis=1)
        sample_obs = pd.concat([sample_obs, sample_obs_categorical], axis=1)
    if continuous_covariates is not None:
        #  Check continuout type
        covs_dtypes = adata.obs[continuous_covariates].dtypes.values
        assert all([pd.api.types.is_numeric_dtype(x) for x in covs_dtypes]
                   ), 'some continuous_covariates are not numeric. Add to categorical_covariates instead'
        # Get sample-level obs
        sample_obs_continuous = adata.obs[[sample_col] + continuous_covariates].drop_duplicates(
        ).reset_index(drop=True).set_index(sample_col)
        # Make design matrix
        design_mat_continuous = sample_obs_continuous.copy()
        assert design_mat_continuous.index.is_unique, 'Sample covariates are not unique.'
        design_mat = pd.concat([design_mat, design_mat_continuous], axis=1)
        sample_obs = pd.concat([sample_obs, sample_obs_continuous], axis=1)
    nhood_adata.var = sample_obs.loc[nhood_adata.var_names].copy()
    nhood_adata.varm['design_matrix'] = design_mat.loc[nhood_adata.var_names].copy()
    adata.uns['nhood_adata'] = nhood_adata.copy()
    # return(design_mat)


# --- Test for differences in uncertainties ---


def test_confidence(adata: AnnData,
                    test_covariate: str,
                    method: str,
                    ref_level: str = None) -> None:
    """Test for differences in confidence statistic between condition and controls

    Params:
    ------
    - adata: AnnData object
    - test_covariate: which covariate stores the condition of interest 
    - method: if 'AUROC' replicates are ignored, if 't-test' a difference in means is used 
    """
    try:
        nhood_adata = adata.uns['nhood_adata'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'] not found -- please run nhood_confidence first")

    try:
        confidence_mat = nhood_adata.layers['confidence'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'].layers['confidence'] not found -- please run nhood_confidence first")

    test_vec = _get_test_column(adata, test_covariate, ref_level=ref_level)

    if method == 'AUROC':
        AUROCs = np.apply_along_axis(
            lambda x: _nhood_AUROC(test_vec, x), 1, confidence_mat)
        nhood_adata.obs['confidence_test_statistic'] = AUROCs

    elif method == 't-test':
        group_index = test_vec == 1
        rest_index = test_vec == 0

        mean_group = confidence_mat[:, group_index].mean(1)
        var_group = confidence_mat[:, group_index].var(1)
        ns_group = sum(group_index)

        mean_rest = confidence_mat[:, rest_index].mean(1)
        var_rest = confidence_mat[:, rest_index].var(1)
        ns_rest = sum(rest_index)

        scores, pvals = stats.ttest_ind_from_stats(
            mean1=mean_group,
            std1=np.sqrt(var_group),
            nobs1=ns_group,
            mean2=mean_rest,
            std2=np.sqrt(var_rest),
            nobs2=ns_rest,
            equal_var=False,  # Welch's
        )
        pvals[np.isnan(pvals)] = 1

        # BH correction
        _, pvals_adj, _, _ = multipletests(
            pvals, alpha=0.05, method='fdr_bh'
        )

        nhood_adata.obs['confidence_test_statistic'] = scores
        nhood_adata.obs['confidence_test_pvals'] = pvals
        nhood_adata.obs['confidence_test_adj_pvals'] = pvals_adj
    else:
        raise ValueError("method must be either 'AUROC' or 't-test'")

    # Save params
    nhood_adata.uns['test_confidence'] = {
        'method': method, 'test_covariate': test_covariate, 'ref_level': ref_level}

    adata.uns['nhood_adata'] = nhood_adata.copy()


def _nhood_AUROC(test_vec, nhood_confidence):
    '''Compute AUROC for single neighbourhood'''
    fpr, tpr, _ = roc_curve(test_vec, nhood_confidence)
    AUROC = auc(fpr, tpr)
    return(AUROC)


def _get_test_column(adata, test_covariate, ref_level=None):
    '''Extract column for comparison in design matrix
    '''
    try:
        nhood_adata = adata.uns['nhood_adata'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'] not found -- please run nhood_confidence first")

    try:
        design_mat = nhood_adata.varm['design_matrix'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'].varm['design_matrix'] not found -- please run make_design first")

    if pd.api.types.is_numeric_dtype(nhood_adata.var[test_covariate].dtype):
        comp_column = test_covariate
    else:
        if ref_level is None:
            # Take last category
            ref_level = nhood_adata.var[test_covariate].cat.categories[-1]
        comp_column = test_covariate + "_" + ref_level
    test_vec = design_mat[comp_column].values
    return(test_vec)


def _get_confidence_stats(nhood_adata, func, test_covariate, cov_level):
    '''
    Get mean confidence for a covariate level
    '''
    ixs = nhood_adata.var[test_covariate] == cov_level
    if func == 'mean':
        s = nhood_adata[:, ixs].layers['confidence'].mean(1)
    elif func == 'var':
        s = nhood_adata[:, ixs].layers['confidence'].var(1)
    # nhood_adata.obs[f'mean_confidence_{test_covariate}{cov_level}'] = mean_conf
    return(s)
