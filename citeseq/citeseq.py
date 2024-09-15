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

adata_GEX = ad.read_h5ad('/workspace/JAMIE/data/citeseq_missing_gex.h5ad')
adata_ADT = ad.read_h5ad('/workspace/JAMIE/data/citeseq_missing_adt.h5ad')

adata_GEX_train = adata_GEX[SITE1_CELL: SITE1_CELL + SITE2_CELL].copy()
adata_GEX_impute = adata_GEX[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL].copy()
adata_ADT_train = adata_ADT[SITE1_CELL: SITE1_CELL + SITE2_CELL].copy()
adata_ADT_impute = adata_ADT[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL].copy()

data1 = adata_GEX_train.X.toarray()
data2 = adata_ADT_train.X.toarray()
data3 = adata_GEX_impute.X.toarray()
data4 = adata_ADT_impute.X.toarray()
corr = np.eye(data1.shape[0], data2.shape[0])

jm = JAMIE(min_epochs=500)
integrated_data = jm.fit_transform(dataset=[data1, data2], P=corr)
jm.save_model('/workspace/JAMIE/ckpts/model.h5')

imputed = jm.modal_predict(data4, 1)
imputed_csr = csr_matrix(data1_imputed)

save_npz('/workspace/JAMIE/data/GEX_imputed_site4.npz', data1_imputed_csr)
