# # Test that giving covariates of the wrong type triggers the right errors

# # Test that adata.uns['nhood_adata'].var contains only expected columns
# import itertools
# all_covariates = list(set(itertools.chain.from_iterable((continuous_covariates or [], categorical_covariates or []))))
# assert all(adata.uns['nhood_adata'].var.columns.isin(all_covariates))
