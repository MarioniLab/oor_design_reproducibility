from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
import anndata
import pandas as pd
import numpy as np
import scanpy as sc

from scipy.sparse import csr_matrix
from scipy import stats
from statsmodels.stats.multitest import multipletests


from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LinearRegression, LogisticRegression
import statsmodels.api as sm

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

import multiprocessing


# def _check_inputs(adata: AnnData) -> bool:
#     """
#     Checks that all the required inputs are there
#     - Cell-level confidence metric (stored in `adata.obs['confidence']`)
#     - Common dimensions for controls and condition samples (stored in `adata.obsm['X_dims']`)
#     - Assignment of cells to sample/donors (stored in `adata.obs['sample_id']`)
#     - Assignment of cells to `adata.obsm['cell_nhoods']`, binary of dimensions cells x nhoods
#     """

# --- Assign cells to populations of similar cells --- #


def make_cell_nhoods(adata: AnnData,
                     method: str,
                     #  X_dims: str,
                     neighbors_key: str = None,
                     clusters_key: str = None,
                     key_added: str = None
                     ) -> None:
    """Assign cells to populations of similar cells in the phenotypic space defined by `X_dims`

   Params:
   ------
   - adata: AnnData object
    - X_dims: string of name of slot in `adata.obsm` storing the common phenotypic space to use for analysis
   - method: if `KNN`, compute milo neighbourhoods based on graph stored in `adata.obsp['neighbors_key']`, if `clusters` assign cells to clusters stored in `adata.obs['clusters_keys']`
   - neighbors_key: only used if `method = "KNN"`, which key in `adata.obsp` to use for assignment of cells to neighbourhoods
   - clusters_key: only used if `method = "clusters"`, which key in `adata.obsp` to use for assignment of cells to clusters
   - key_added: If not specified, the neighbourhoods are stored in .obsm['nhoods']. If specified, they are stored as .obsm[key_added+'_nhoods']

   Returns:
   -------
   None, adds in place nhood matric to `adata.obsm`, a binary sparse matrix of dimensions cells x neighbourhoods
    """
    if method == 'clusters':
        nhoods_mat = _make_nhoods_clusters(adata, clusters_key=clusters_key)
    elif method == 'KNN':
        None
    else:
        raise ValueError("method must be either 'KNN' or 'clusters'")
    if key_added is None:
        adata.obsm['nhoods'] = nhoods_mat
    else:
        adata.obsm[key_added + '_nhoods'] = nhoods_mat


def _make_nhoods_KNN(adata: AnnData, neighbors_key: str = None):
    """make milo neighbourhoods"""
    None


def _make_nhoods_clusters(adata: AnnData, clusters_key: str = None):
    """make neighbourhoods from cluster assignment"""
    try:
        clusters = adata.obs[clusters_key]
    except KeyError:
        raise ValueError('clusters_key is not a column in adata.obs')
    return(csr_matrix(pd.get_dummies(clusters).values))


# --- Compute confidence metric per sample/donor --- #

def nhood_confidence(adata: AnnData,
                     confidence_col: str,
                     sample_col: str,
                     nhoods_key: str = None,
                     impute_missing: bool = True,
                     min_cells: int = 1
                     ) -> None:
    """Aggregate cell-level confidence by cell neighourhood and sample
    Params:
    ------
    - adata: AnnData object
    - confidence_col: column in `adata.obs` storing cell-level confidence statistic. Higher values correspond to higher confidence in reference mapping.
    - sample_col: column in `adata.obs` storing assignment of cells to samples
    - nhoods_key: If not specified, nhood_confidence uses .obsm['nhoods'] (default storage place for make_cell_nhoods).
        If specified, nhood_confidence uses .obsm[key_added + '_nhoods']
    - impute_missing: boolean indicating whether confidence for missing samples (i.e. where there are no cells in the nhood) should be set to 0 (default: True)
    - min_cells: integer indicating the minimum number of cells per sample to compute mean confidence (default: 1)

    Returns:
    -------
    None, adds in place `adata.uns['nhood_adata']` of dimensions nhoods x samples storing cell counts in .X and confidence in `layers['confidence']`
    """
    if nhoods_key is None:
        nhoods_key = 'nhoods'
    else:
        nhoods_key = nhoods_key + '_nhoods'
    try:
        nhood_mat = adata.obsm[nhoods_key].copy()  #  cells x nhoods
    except:
        raise KeyError(
            f"adata.obsm['{nhoods_key}'] not found, please run make_cell_nhoods first")
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
        if impute_missing:
            nh_score[np.isnan(nh_score)] = 0
        nhood_scores = np.hstack([nhood_scores, nh_score])

    # Make anndata object for nhoods x samples if missing
    if 'nhood_adata' not in adata.uns:
        _add_nhoods_adata(adata, sample_col, nhoods_key)
        adata.uns['nhood_adata'].layers['confidence'] = nhood_scores
    elif adata.uns['nhood_adata'].uns['sample_col'] != sample_col or adata.uns['nhood_adata'].uns['nhoods_key'] != nhoods_key:
        _add_nhoods_adata(adata, sample_col)
        adata.uns['nhood_adata'].layers['confidence'] = nhood_scores
    else:
        adata.uns['nhood_adata'].layers['confidence'] = nhood_scores
    if min_cells > 1:
        adata.uns['nhood_adata'].layers['confidence'][adata.uns['nhood_adata'].X.toarray(
        ) < min_cells] = np.nan


def _add_nhoods_adata(
    adata: AnnData,
    sample_col: str,
    nhoods_key: str = 'nhoods'
):
    '''
    - adata
    - sample_col: string, column in adata.obs that contains sample information
    (what should be in the columns of the nhoodCount matrix)
    - nhoods_key: If not specified, nhood_confidence uses .obsm['nhoods'] (default storage place for make_cell_nhoods).
        If specified, nhood_confidence uses .obsm[key_added + '_nhoods']

    Returns: None
    Updated adata.uns slot to contain adata.uns["nhood_adata"], where:
    - adata.uns["nhood_adata"].obs_names are neighbourhoods
    - adata.uns["nhood_adata"].var_names are samples
    - adata.uns["nhood_adata"].X is the matrix counting the number of cells from each
    sample in each neighbourhood
    '''
    try:
        nhoods = adata.obsm[nhoods_key]
    except KeyError:
        raise KeyError(
            f'Cannot find {nhoods_key} slot in adata.obsm -- please run make_cell_nhoods'
        )
    #  Make nhood abundance matrix
    sample_dummies = pd.get_dummies(adata.obs[sample_col])
    all_samples = sample_dummies.columns
    sample_dummies = csr_matrix(sample_dummies.values)
    nhood_count_mat = adata.obsm[nhoods_key].T.dot(sample_dummies)
    nhood_var = pd.DataFrame(index=all_samples)
    nhood_adata = anndata.AnnData(X=nhood_count_mat, var=nhood_var)
    nhood_adata.uns["sample_col"] = sample_col
    # Save nhood index info if milo nhoods
    if "nhood_ixs_refined" in adata.obs.columns:
        nhood_adata.obs["index_cell"] = adata.obs_names[adata.obs["nhood_ixs_refined"] == 1]
        nhood_adata.obs["kth_distance"] = adata.obs.loc[adata.obs["nhood_ixs_refined"]
                                                        == 1, "nhood_kth_distance"].values
    nhood_adata.uns['nhoods_key'] = nhoods_key
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
        sample_obs_categorical = sample_obs_categorical.apply(
            lambda x: x.astype('category'))
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
    nhood_adata.varm['design_matrix'] = nhood_adata.varm['design_matrix'].astype(
        'float')
    adata.uns['nhood_adata'] = nhood_adata.copy()
    # return(design_mat)

# --- Filter nhoods ---


def highly_variable_nhoods(adata: AnnData):
    '''
    Find neighbourhoods with highly variable confidence scores

    Params:
    ------
    - adata: AnnData object

    Returns:
    -------
    None, adds in place adata.uns['nhood_adata'].obs['highly_variable']
    '''
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

    nhood_conf_var = np.nanvar(nhood_adata.layers['confidence'], axis=1)
    nhood_conf_mean = np.nanmean(nhood_adata.layers['confidence'], axis=1)

    lowess = sm.nonparametric.lowess
    z = lowess(nhood_conf_var, nhood_conf_mean)
    # Keep neighbourhoods with std higher than expected by the mean
    nhood_adata.obs['confidence_mean'] = nhood_conf_mean
    nhood_adata.obs['confidence_var'] = nhood_conf_var
    nhood_adata.obsm['confidence_var_lowess'] = z
    nhood_adata.obs['highly_variable'] = nhood_conf_var > z[:, 1]
    adata.uns['nhood_adata'] = nhood_adata.copy()

# --- Test for differences in uncertainties ---


def test_confidence(adata: AnnData,
                    test_covariate: str,
                    method: str,
                    ref_level: str = None,
                    alpha: float = 0.05,
                    use_highly_variable: bool = False
                    ) -> None:
    """Test for differences in confidence statistic between condition and controls

    Params:
    ------
    - adata: AnnData object
    - test_covariate: which covariate stores the condition of interest 
    - method: if 'AUROC' replicates are ignored, if 't-test' a difference in means is used,
    if 'OLS' also 
    - alpha: signif thresh for uncorrected pvals (used for BH correction in OLS and t-test)
    - use_highly_variable: boolean indicating whether the test should be performed only on highly variable nhoods 
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

    if use_highly_variable:
        keep_nhoods = nhood_adata.obs['highly_variable'].values
        confidence_mat = confidence_mat[keep_nhoods, :].copy()
    else:
        keep_nhoods = np.ones(nhood_adata.n_obs, dtype=bool)

    test_vec = _get_test_column(adata, test_covariate, ref_level=ref_level)

    nhood_adata.obs['confidence_test_statistic'] = np.nan
    nhood_adata.obs['confidence_test_pvals'] = np.nan
    nhood_adata.obs['confidence_test_adj_pvals'] = np.nan

    if method == 'AUROC':
        AUROCs = np.apply_along_axis(
            lambda x: _nhood_AUROC(test_vec, x), 1, confidence_mat)
        conf_stat = AUROCs
        pvalues = None
        fdr = None

    elif method == 'OLS':
        outs = np.apply_along_axis(
            lambda x: _nhood_OLS(test_vec, x), 1, confidence_mat)
        conf_stat = outs[:, 0]
        pvalues = outs[:, 1]
        pvalues[np.isnan(outs[:, 0])] = np.nan

    elif method == 't-test':
        group_index = test_vec == 1
        rest_index = test_vec == 0

        mean_group = np.nanmean(confidence_mat[:, group_index], axis=1)
        var_group = np.nanvar(confidence_mat[:, group_index], axis=1)
        ns_group = sum(group_index)

        mean_rest = np.nanmean(confidence_mat[:, rest_index], axis=1)
        var_rest = np.nanvar(confidence_mat[:, rest_index], axis=1)
        ns_rest = sum(rest_index)

        conf_stat, pvalues = stats.ttest_ind_from_stats(
            mean1=mean_group,
            std1=np.sqrt(var_group),
            nobs1=ns_group,
            mean2=mean_rest,
            std2=np.sqrt(var_rest),
            nobs2=ns_rest,
            equal_var=False,  # Welch's
        )
    else:
        raise ValueError("method must be 'AUROC', 'OLS' or 't-test'")

    nhood_adata.obs.loc[keep_nhoods, 'confidence_test_statistic'] = conf_stat

    # Multiple testing correction
    if pvalues is not None:
        nhood_adata.obs.loc[keep_nhoods,
                            'confidence_test_pvals'] = pvalues.copy()
        pvalues = nhood_adata.obs['confidence_test_pvals'].values.copy()
        keep_nhoods = ~np.isnan(pvalues)
        p_rank = stats.rankdata(pvalues[keep_nhoods])
        fdr = pvalues[keep_nhoods] * len(pvalues[keep_nhoods]) / p_rank
        fdr[fdr > 1] = 1

        nhood_adata.obs.loc[keep_nhoods,
                            "confidence_test_adj_pvals"] = fdr
    # Save params
    nhood_adata.uns['test_confidence'] = {
        'method': method,
        'test_covariate': test_covariate,
        'ref_level': np.where(ref_level is not None, ref_level, nhood_adata.var[test_covariate].cat.categories[0]).flatten()[-1]}

    adata.uns['nhood_adata'] = nhood_adata.copy()


def _nhood_AUROC(test_vec, nhood_confidence):
    '''Compute AUROC for single neighbourhood'''
    nan_mask = ~np.isnan(nhood_confidence)
    fpr, tpr, _ = roc_curve(test_vec[nan_mask], nhood_confidence[nan_mask])
    AUROC = auc(fpr, tpr)
    return(AUROC)


def _nhood_OLS(test_vec, nhood_confidence):
    '''Compute OLS regression for single neighbourhood'''
    nan_mask = ~np.isnan(nhood_confidence)
    y = nhood_confidence[nan_mask]
    X = test_vec[nan_mask]
    # X = test_vec[nan_mask].reshape(-1, 1)
    coef, _, _, pval, _ = stats.linregress(X, y)
    # ols = LinearRegression().fit(X, y)
    # coef = ols.coef_[0]
    # X = sm.add_constant(X)
    # pval = sm.OLS(y, X).fit().summary2().tables[1]['P>|t|'][0]
    return(coef, pval)


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


def augur_nhoods(adata,
                 test_covariate: str,
                 ref_level: str = None,
                 n_folds: int = 5,
                 sample_size: int = 20,
                 n_workers: int = 5
                 ) -> None:
    """Compute Augur

    Params:
    ------
    - adata: AnnData object
    - test_covariate: which covariate stores the condition of interest 
    - ref_level: which level of test_covariate to set to 1
    - n_folds: how many permutations to compute AUC on
    - sample_size: how many cells from each condition to downsample to 
        (crucial, otherwise AUC dependent on no of cells from each condition in nhood) 
    - n_workers: for multiprocessing 
    """
    try:
        nhood_adata = adata.uns['nhood_adata'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'] not found -- please run nhood_confidence first")

    Y = (adata.obs[test_covariate]
         == ref_level).astype('int').values

    if 'log1p' not in adata.uns.keys():
        adata.layers['counts'] = adata.X.copy()
        sc.pp.normalize_per_cell(adata)
        sc.pp.log1p(adata)

    X = adata.X.copy()

    # Make array w nhood x cell indexes (sample subsample of cells x nhood)
    subsamples = []
    for i in range(adata.obsm['nhoods'].shape[1]):
        pos_ixs = np.where((_get_cells_from_nhood(
            adata, i)) & (Y == 1))[0]
        neg_ixs = np.where((_get_cells_from_nhood(
            adata, i)) & (Y == 0))[0]

        pos_ixs = np.random.choice(pos_ixs, sample_size)
        neg_ixs = np.random.choice(neg_ixs, sample_size)
        subsamples.append(np.hstack([pos_ixs, neg_ixs]))
    nhoods_sample_ixs = np.vstack(subsamples)

    pool = multiprocessing.Pool(n_workers)
    inputs = [(Y[nhoods_sample_ixs[i, :]], X[nhoods_sample_ixs[i, :], :], 0.33)
              for i in range(nhoods_sample_ixs.shape[0])]

    # Do 5-fold CV dei poveri
    aucs = np.array(pool.map(_nhood_augur, inputs))
    for i in range(n_folds - 1):
        print(f'it {i+1}')
        aucs_i = np.array(pool.map(_nhood_augur, inputs))
        aucs = np.vstack([aucs, aucs_i])
    _get_nhood_adata(adata).obs['Augur_AUC'] = aucs.mean(0)


def _nhood_augur(params):
    Y_nhood, X_nhood, test_size = params

    # Feature selection
    hvgs = sc.pp.highly_variable_genes(anndata.AnnData(X_nhood), inplace=False)
    X_nhood = X_nhood[:, hvgs['highly_variable']]

    #  Split in train and test with subsampling
    X_train, X_test, Y_train, Y_test = train_test_split(
        X_nhood, Y_nhood, test_size=test_size)

    # Train classifier on train
    est = LogisticRegression(penalty='l2', random_state=None)
    est.fit(X_train, Y_train)

    ## Predict in test
    Y_pred = est.predict(X_test)

    #  Compute AUC
    nhood_auc = roc_auc_score(Y_test, Y_pred)
    return(nhood_auc)


def _get_nhood_adata(adata):
    return(adata.uns['nhood_adata'])


def _get_cells_from_nhood(adata, i):
    return((adata.obsm['nhoods'][:, i].toarray() == 1).flatten())
