import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, save_npz
from jamie import JAMIE
import time

SITE1_CELL = 16311
SITE2_CELL = 25171
SITE3_CELL = 32029
SITE4_CELL = 16750
GEX = 2000

adata_GEX = ad.read_h5ad('/workspace/JAMIE/data/multiome_missing_gex.h5ad')
adata_ATAC = ad.read_h5ad('/workspace/JAMIE/data/multiome_missing_atac.h5ad')

adata_GEX_train = adata_GEX[:SITE1_CELL].copy()
adata_GEX_impute = adata_GEX[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL].copy()
adata_ATAC_train = adata_ATAC[:SITE1_CELL].copy()
adata_ATAC_impute = adata_ATAC[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL].copy()

data1 = adata_GEX_train.X.toarray()
data2 = adata_ATAC_train.X.toarray()
data3 = adata_GEX_impute.X.toarray()
data4 = adata_ATAC_impute.X.toarray()
corr = np.eye(data1.shape[0], data2.shape[0])

jm = JAMIE(min_epochs=500)
start_time = time.time()
integrated_data = jm.fit_transform(dataset=[data1, data2], P=corr)
end_time = time.time()
print(f"Model training took {end_time - start_time:.2f} seconds")
jm.save_model('/workspace/JAMIE/ckpts/model.h5')

start_time = time.time()
imputed = jm.modal_predict(data4, 1)
end_time = time.time()
print(f"Model predicting took {end_time - start_time:.2f} seconds")
print(imputed.shape)

imputed_csr = csr_matrix(imputed)
save_npz('/workspace/JAMIE/data/multiome_GEX_imputed_site3.npz', imputed_csr)
