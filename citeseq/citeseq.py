import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, save_npz
from jamie import JAMIE


SITE1_CELL = 16311
SITE2_CELL = 25171
SITE3_CELL = 32029
SITE4_CELL = 16750
GEX = 2000


citeseq = ad.read_h5ad("/workspace/JAMIE/data/citeseq_processed.h5ad")
citeseq.var_names_make_unique()

adata_GEX = citeseq[:, citeseq.var["feature_types"] == "GEX"].copy()
adata_ADT = citeseq[:, citeseq.var["feature_types"] == "ADT"].copy()
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

adata_GEX = adata_GEX[-SITE4_CELL:]
adata_ADT = adata_ADT[-SITE4_CELL:]

data1 = adata_GEX.X.toarray()
data1[:, :GEX] = 0
data2 = adata_ADT.X.toarray()
corr = np.eye(data1.shape[0], data2.shape[0])

jm = JAMIE(min_epochs=500)
integrated_data = jm.fit_transform(dataset=[data1, data2], P=corr)

data1_imputed = jm.modal_predict(data2, 1)
data1_imputed_csr = csr_matrix(data1_imputed)

save_npz('/workspace/JAMIE/data/GEX_imputed_site4.npz', data1_imputed_csr)
