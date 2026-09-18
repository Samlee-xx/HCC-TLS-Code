###### TLS deconvolution by cell2location ######

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import cell2location

np.random.seed(777)

# 薛老师的nature数据
adata_ref = sc.read_h5ad(
    r"E:\01-TLS\01-Omics\01-ST\02-New data\02-analysis-copy\02-cell2location\sc_trained.h5ad"
)

factor_names = list(adata_ref.uns["mod"]["factor_names"])

if "means_per_cluster_mu_fg" in adata_ref.varm:
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{x}" for x in factor_names]
    ].copy()
else:
    inf_aver = adata_ref.var[
        [f"means_per_cluster_mu_fg_{x}" for x in factor_names]
    ].copy()

inf_aver.columns = factor_names

print(inf_aver.shape)
print(inf_aver.columns)

###### TLS raw mini-bulk counts ######
# 需要最原始的count矩阵
tls_count = pd.read_csv(
    r"E:\01-TLS\01-Omics\01-ST\02-New data\02-analysis\03-Cell2location\TLS count ensemble_without batch.csv",
    index_col=0
)
# TLS的表型数据读入
tls_info = pd.read_csv(
    r"E:\01-TLS\01-Omics\01-ST\04-analysis final\08-Alternative_deconvolution\inputs\tls_metadata_461.csv"
)

tls_count = tls_count.loc[:, tls_info["TLS_ID"]]

print(tls_count.shape)
print(tls_info["sample"].value_counts())

adata_tls = ad.AnnData(tls_count.T)
adata_tls.var_names_make_unique()
adata_tls.obs_names_make_unique()

tls_info.index = tls_info["TLS_ID"]
tls_info = tls_info.loc[adata_tls.obs_names]
adata_tls.obs["sample"] = tls_info["sample"].astype(str)

common_gene = np.intersect1d(adata_tls.var_names, inf_aver.index)
adata_tls = adata_tls[:, common_gene].copy()
inf_aver = inf_aver.loc[common_gene, :].copy()

print(adata_tls.shape)

###### cell2location model ######

cell2location.models.Cell2location.setup_anndata(
    adata=adata_tls,
    batch_key="sample"
)

model = cell2location.models.Cell2location(
    adata_tls,
    cell_state_df=inf_aver,
    N_cells_per_location=30,
    detection_alpha=20
)

model.train(
    max_epochs=5000,
    batch_size=None,
    train_size=1
)

adata_tls = model.export_posterior(
    adata_tls,
    sample_kwargs={"num_samples": 1000, "batch_size": adata_tls.n_obs}
)

cell_abundance = adata_tls.obsm["means_cell_abundance_w_sf"]

if not isinstance(cell_abundance, pd.DataFrame):
    cell_abundance = pd.DataFrame(
        cell_abundance,
        index=adata_tls.obs_names,
        columns=factor_names
    )

print(cell_abundance.shape)
print(cell_abundance.head())
