### Tutorial

We have uploaded pkl models for every cell lineage, download the pkl file and celltypist module, follow the following guidance:

##### 1.Import the modules


```python
import scanpy as sc
import celltypist
import time
import numpy as np
```

##### 2.Read in the query h5ad file


```python
adata_query=sc.read_h5ad('query.h5ad')
```

(optional)do these two process if data was not log-transformed before.


```python
# sc.pp.normalize_total(adata_query, target_sum = 1e4)
# sc.pp.log1p(adata_query)
```

    WARNING: adata.X seems to be already log-transformed.


##### 3.To predict the query adata from the model, use the `over_clustering` parameter if you have pre-clustered the data. The `majority_voting` option will perform voting within the clusters. 


```python
t_start = time.time()
predictions = celltypist.annotate(adata_query, model = './model_JHLineage_ref.pkl', over_clustering='leiden', majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")
```

    🔬 Input data has 168971 cells and 19019 genes
    🔗 Matching reference genes in the model
    🧬 2957 features used for prediction
    ⚖️ Scaling input data
    🖋️ Predicting labels
    ✅ Prediction done!
    👀 Detected a neighborhood graph in the input object, will run over-clustering on the basis of it
    ⛓️ Over-clustering input data with resolution set to 25
    🗳️ Majority voting the predictions
    ✅ Majority voting done!


    Time elapsed: 1580.080134153366 seconds


##### 4. Return the prediction to the adata_query object. Specifically, the `adata.obs` DataFrame will be updated with four new columns:

`predicted_labels`: The predicted cell type labels for each cell.
`over_clustering`: The clustering results if over-clustering was applied.
`majority_voting`: The results of the majority voting within clusters.
`conf_score`: The confidence scores for the predictions.


```python
adata_query= predictions.to_adata()
```


```python
adata_query
```


    AnnData object with n_obs × n_vars = 168971 × 19019
        obs: 'doublet_scores', 'predicted_doublets', 'project', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'S_score', 'G2M_score', 'phase', 'leiden', 'Lineage', 'predicted_labels', 'over_clustering', 'majority_voting', 'conf_score'
        var: 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection'
        uns: 'Lineage_colors', 'hvg', 'leiden', 'leiden_colors', 'log1p', 'neighbors', 'pca', 'project_colors', 'rank_genes_groups', 'umap'
        obsm: 'X_pca', 'X_umap'
        obsp: 'connectivities', 'distances'


##### 5.Visualization


```python
sc.pl.umap(adata_query,color=['Lineage','majority_voting','predicted_labels'],ncols=2,wspace=0.3)
```
