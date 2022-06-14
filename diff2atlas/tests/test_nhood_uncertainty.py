# Test that the size of adata.uns['nhood_adata'] is correct

# Test that the confidence is zero when the number of cells from a sample in a neighbourhood is zero

# Test that mean values are within the expected range

# def test_nhood_mean_range(adata):
#     annotate_nhoods_continuous(adata, anno_col='S_score')
#     assert adata.uns['nhood_adata'].obs['nhood_S_score'].max(
#     ) < adata.obs['S_score'].max()
#     assert adata.uns['nhood_adata'].obs['nhood_S_score'].min(
#     ) > adata.obs['S_score'].min()


# ## Test that value corresponds to mean

# def test_correct_mean(adata):
#     annotate_nhoods_continuous(adata, anno_col='S_score')
#     i = np.random.choice(np.arange(adata.uns['nhood_adata'].n_obs))
#     mean_val_nhood = adata.obs[adata.obsm['nhoods']
#                                [:, i].toarray() == 1]['S_score'].mean()
#     assert adata.uns['nhood_adata'].obs['nhood_S_score'][i] == pytest.approx(
#         mean_val_nhood, 0.0001)
