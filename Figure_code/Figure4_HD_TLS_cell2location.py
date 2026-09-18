###### Cell2location deconvolution of Visium HD TLSs ######

from pathlib import Path

import anndata as ad
import cell2location
import numpy as np
import pandas as pd
import scanpy as sc


###### files ######

input_dir = Path(
    r"F:\TLS DATA\03-CODEX-ref_to_HD\03-analysis-final"
    r"\04-Correlation analysis\01-HD_TLS_region_count_sum"
)
output_dir = Path(
    r"F:\TLS Material\06-Code\result\Figure4_CODEX_HD_integration"
    r"\01-HD_TLS_cell2location"
)
output_dir.mkdir(parents=True, exist_ok=True)

count_file = input_dir / "HD_TLS_region_gene_count_sum.csv"
metadata_file = input_dir / "HD_TLS_region_metadata.csv"

# This trained reference was generated from PRJCA007744.
reference_file = Path(
    r"E:\01-TLS\01-Omics\01-ST\02-New data\02-analysis-copy"
    r"\02-cell2location\sc_trained.h5ad"
)


###### read the 123 Visium HD TLS regions ######

count_table = pd.read_csv(count_file)
tls_metadata = pd.read_csv(metadata_file).set_index("tls_id")

tls_ids = tls_metadata.index.tolist()

raw_count = count_table.set_index("gene_id")[tls_ids]
if raw_count.index.duplicated().any():
    raise ValueError("Duplicated gene_id values were found")

adata_tls = ad.AnnData(raw_count.T)
adata_tls.var = count_table.set_index("gene_id")[
    ["gene_name", "feature_type", "genome"]
].copy()
adata_tls.obs = tls_metadata.loc[adata_tls.obs_names].copy()

# All TLS mini-bulk profiles were fitted together as one Visium HD dataset.
adata_tls.obs["cell2location_batch"] = "Visium_HD"



###### read the PRJCA007744 reference signatures ######

adata_reference = sc.read_h5ad(reference_file, backed="r")
cell_types = list(adata_reference.uns["mod"]["factor_names"])

if "means_per_cluster_mu_fg" in adata_reference.varm:
    reference_signatures = adata_reference.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{x}" for x in cell_types]
    ].copy()
else:
    reference_signatures = adata_reference.var[
        [f"means_per_cluster_mu_fg_{x}" for x in cell_types]
    ].copy()

reference_signatures.columns = cell_types
adata_reference.file.close()

shared_genes = np.intersect1d(
    adata_tls.var_names,
    reference_signatures.index
)
if len(shared_genes) < 1000:
    raise ValueError("Too few genes overlap with the reference")

adata_tls = adata_tls[:, shared_genes].copy()
reference_signatures = reference_signatures.loc[shared_genes].copy()


###### cell2location ######

cell2location.models.Cell2location.setup_anndata(
    adata=adata_tls,
    batch_key="cell2location_batch"
)

model = cell2location.models.Cell2location(
    adata_tls,
    cell_state_df=reference_signatures,
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
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": adata_tls.n_obs
    }
)


###### save cell-type abundance and proportion ######

mean_abundance = pd.DataFrame(
    np.asarray(adata_tls.obsm["means_cell_abundance_w_sf"]),
    index=adata_tls.obs_names,
    columns=cell_types
)

mean_proportion = mean_abundance.div(
    mean_abundance.sum(axis=1),
    axis=0
)

mean_abundance.to_csv(
    output_dir / "HD_TLS_cell2location_mean_abundance.csv",
    index_label="tls_id"
)
mean_proportion.to_csv(
    output_dir / "HD_TLS_cell2location_mean_proportion.csv",
    index_label="tls_id"
)

model.save(str(output_dir / "model"), overwrite=True)
adata_tls.write(str(output_dir / "HD_TLS_cell2location.h5ad"))

qc = pd.DataFrame([{
    "n_HD_TLS": adata_tls.n_obs,
    "n_input_genes": count_table.shape[0],
    "n_shared_genes": len(shared_genes),
    "n_reference_cell_types": len(cell_types)
}])
qc.to_csv(output_dir / "HD_TLS_cell2location_QC.csv", index=False)
