import anndata as ad
import pandas as pd
import scanpy as sc
import sys
from scipy.sparse import csr_matrix


SITE1_CELL = 16311
SITE2_CELL = 25171
SITE3_CELL = 32029
SITE4_CELL = 16750


adata = ad.read_h5ad('/workspace/JAMIE/data/citeseq_processed.h5ad')
adata.var_names_make_unique()

adata_GEX = adata[:, adata.var["feature_types"] == "GEX"].copy()
adata_ADT = adata[:, adata.var["feature_types"] == "ADT"].copy()
### step 1
sc.pp.normalize_total(adata_GEX, target_sum=1e4)
sc.pp.normalize_total(adata_ADT, target_sum=1e4)
### step 2
sc.pp.log1p(adata_GEX)
sc.pp.log1p(adata_ADT)
### step 3
sc.pp.highly_variable_genes(
    adata_GEX,
    n_top_genes=2000,
    subset=True
)

X_GEX = adata_GEX.X.toarray()
X_GEX[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL] = 0
adata_GEX.X = csr_matrix(X_GEX)

adata_GEX.write_h5ad('/workspace/JAMIE/data/citeseq_missing_gex.h5ad')
adata_ADT.write_h5ad('/workspace/JAMIE/data/citeseq_missing_adt.h5ad')