###### Hierarchical clustering of 154 CODEX-defined TA-TLSs ######

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler


###### files ######

input_file = Path(
    r"F:\TLS DATA\03-CODEX-ref_to_HD\03-Analysis"
    r"\02-TLS-DBSCAN\02_Data_Output"
    r"\DBSCAN_TLS_summary_filtered_eps80_min50.csv"
)

initial_qc_file = Path(
    r"F:\TLS Material\06-Code"
    r"\Figure4_CODEX_initial_QC_154_TLS_ids.csv"
)

output_dir = Path(
    r"F:\TLS Material\06-Code\result\Figure4_CODEX_HD_integration"
    r"\00-CODEX_154_TLS_clustering"
)



###### read the DBSCAN CODEX TLSs ######

initial_qc_tls_ids = pd.read_csv(initial_qc_file)["tls_id"].tolist()
tls = pd.read_csv(input_file).set_index("tls_id")
tls = tls.loc[initial_qc_tls_ids].reset_index()
tls = tls.loc[
    (tls["pct_unassign"] <= 30)
    & (tls["pct_stroma"] <= 60)
    & (tls["pct_lymphoid"] >= 10)
].copy()


###### cell-type percentages among annotated TLS cells ######

defined_cells = tls["n_total"] - tls["n_unassign"]
if (defined_cells <= 0).any():
    raise ValueError("A TLS contains no annotated cells")

tls["dpct_B"] = tls["n_B_cell"] / defined_cells * 100
tls["dpct_CD4T"] = tls["n_CD4T"] / defined_cells * 100
tls["dpct_CD8T"] = tls["n_CD8T"] / defined_cells * 100
tls["dpct_T"] = tls["n_T"] / defined_cells * 100
tls["dpct_Neutrophil"] = tls["n_neutrophil"] / defined_cells * 100
tls["dpct_Macrophage"] = tls["n_macrophage"] / defined_cells * 100
tls["dpct_Endothelial_Lymphatic"] = (
    tls["n_endothelial"] + tls["n_lymphatic"]
) / defined_cells * 100
tls["dpct_Fibroblast"] = tls["n_fibroblast"] / defined_cells * 100
tls["dpct_Cholangiocyte"] = tls["n_bile_duct"] / defined_cells * 100


###### hierarchical clustering ######

clustering_features = [
    "dpct_B",
    "dpct_T",
    "dpct_Neutrophil",
    "dpct_Macrophage",
    "dpct_Endothelial_Lymphatic",
    "dpct_Fibroblast",
    "dpct_Cholangiocyte"
]

cell_composition = tls[clustering_features].astype(float)
if not np.isfinite(cell_composition.to_numpy()).all():
    raise ValueError("Non-finite cell-type percentages were found")

scaled_composition = StandardScaler().fit_transform(cell_composition)
cosine_distance = pdist(scaled_composition, metric="cosine")

finite_distance = cosine_distance[np.isfinite(cosine_distance)]
largest_distance = finite_distance.max() if finite_distance.size else 0
cosine_distance = np.nan_to_num(
    cosine_distance,
    nan=largest_distance,
    posinf=largest_distance,
    neginf=0
)

linkage_matrix = linkage(cosine_distance, method="complete")
tls["cluster"] = fcluster(
    linkage_matrix,
    t=2,
    criterion="maxclust"
).astype(str)


###### name the two TLS subtypes ######

bf_cluster = (
    tls.groupby("cluster")["dpct_B"]
    .mean()
    .idxmax()
)

tls["TLS_subtype"] = np.where(
    tls["cluster"].eq(bf_cluster),
    "BF_TLS",
    "TNC_TLS"
)

assert tls["TLS_subtype"].value_counts().to_dict() == {
    "TNC_TLS": 114,
    "BF_TLS": 40
}


###### save results ######

output_columns = [
    "tls_id",
    "batch",
    "slide",
    "fov",
    "tls_number",
    "dbscan_label",
    "n_total",
    "n_unassign",
    "cluster",
    "TLS_subtype",
    "dpct_B",
    "dpct_CD4T",
    "dpct_CD8T",
    "dpct_T",
    "dpct_Neutrophil",
    "dpct_Macrophage",
    "dpct_Endothelial_Lymphatic",
    "dpct_Fibroblast",
    "dpct_Cholangiocyte"
]

tls[output_columns].to_csv(
    output_dir / "CODEX_154_TLS_assignment.csv",
    index=False
)

tls.groupby("TLS_subtype")[clustering_features].mean().to_csv(
    output_dir / "CODEX_154_TLS_mean_cell_composition.csv"
)

pd.DataFrame(
    scaled_composition,
    index=tls["tls_id"],
    columns=clustering_features
).to_csv(
    output_dir / "CODEX_154_TLS_standardized_composition.csv",
    index_label="tls_id"
)

pd.DataFrame(
    linkage_matrix,
    columns=["cluster_1", "cluster_2", "distance", "n_TLS"]
).to_csv(
    output_dir / "CODEX_154_TLS_linkage_matrix.csv",
    index=False
)


###### dendrogram ######

plt.figure(figsize=(10, 4))
dendrogram(linkage_matrix, no_labels=True)
plt.xlabel("CODEX TA-TLSs")
plt.ylabel("Cosine distance")
plt.tight_layout()
plt.savefig(output_dir / "CODEX_154_TLS_dendrogram.pdf")
plt.close()
