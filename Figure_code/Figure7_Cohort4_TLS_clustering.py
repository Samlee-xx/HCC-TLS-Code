###### Cohort 4 TLS subtype clustering ######

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler


###### files ######

input_file = Path(
    r"F:\TLS DATA\05-ICI-CODEX\03-analysis\03-TLS-spatial-analysis"
    r"\01-TLS DBSCAN\03-result\01-tables\tls_summary_final.csv"
)

output_dir = Path(
    r"F:\TLS Material\06-Code\result\Figure6_Cohort4\01-TLS_clustering"
)
output_dir.mkdir(parents=True, exist_ok=True)


###### read the DBSCAN-defined TLSs ######

tls = pd.read_csv(input_file)

print(tls.shape)
print(tls["tls_id"].nunique())
print(tls["n_total"].describe())


###### cell-type percentages used for clustering ######

tls["dpct_B"] = tls["pct_B_cell"]
tls["dpct_CD4T"] = tls["pct_CD4T"]
tls["dpct_CD8T"] = tls["pct_CD8T"]
tls["dpct_T"] = tls["pct_T"]
tls["dpct_Neutrophil"] = tls["pct_Neutrophil"]
tls["dpct_Macrophage"] = tls["pct_Macrophage"]
tls["dpct_Endothelial_Lymphatic"] = tls["pct_Endothelial_Lymphatic"]
tls["dpct_Fibroblast"] = tls["pct_Fibroblast"]
tls["dpct_Cholangiocyte"] = tls["pct_Cholangiocyte"]


###### retain TLSs that passed cell-number ######

tls["pass_cell_number"] = tls["n_total"] >= 50
tls["used_for_clustering"] = (
    tls["pass_cell_number"]
)

tls_used = tls[tls["used_for_clustering"]].copy()


clustering_features = [
    "dpct_B",
    "dpct_T",
    "dpct_Neutrophil",
    "dpct_Macrophage",
    "dpct_Endothelial_Lymphatic",
    "dpct_Fibroblast",
    "dpct_Cholangiocyte",
]

print(tls_used[clustering_features].describe().round(2))

scaled_composition = StandardScaler().fit_transform(
    tls_used[clustering_features]
)

cosine_distance = pdist(scaled_composition, metric="cosine")
linkage_matrix = linkage(cosine_distance, method="complete")

tls_used["cluster"] = fcluster(
    linkage_matrix,
    t=2,
    criterion="maxclust",
)



print(tls_used.groupby("cluster")["dpct_B"].mean())

bf_cluster = tls_used.groupby("cluster")["dpct_B"].mean().idxmax()

tls_used["TLS_subtype"] = np.where(
    tls_used["cluster"] == bf_cluster,
    "BF-TLS",
    "TNC-TLS",
)

print(tls_used["TLS_subtype"].value_counts())
print(tls_used.groupby("TLS_subtype")[clustering_features].mean().round(2))


###### save the assignment ######

output_columns = [
    "tls_id", "batch", "slide", "fov", "tls_number", "n_total",
    "cluster", "TLS_subtype", "dpct_B", "dpct_CD4T", "dpct_CD8T",
    "dpct_T", "dpct_Neutrophil", "dpct_Macrophage",
    "dpct_Endothelial_Lymphatic", "dpct_Fibroblast",
    "dpct_Cholangiocyte",
]

tls_used[output_columns].to_csv(
    output_dir / "Cohort4_TLS_subtype_assignment.csv",
    index=False,
)

tls_used.groupby("TLS_subtype")[clustering_features].mean().to_csv(
    output_dir / "Cohort4_TLS_mean_cell_composition.csv"
)


###### simple plots for manual inspection ######

plt.figure(figsize=(9, 4))
dendrogram(linkage_matrix, no_labels=True)
plt.xlabel("TLS")
plt.ylabel("Cosine distance")
plt.tight_layout()
plt.savefig(output_dir / "Cohort4_TLS_dendrogram.pdf")
plt.close()

plot_order = np.argsort(tls_used["TLS_subtype"].to_numpy())

plt.figure(figsize=(7, 5))
plt.imshow(scaled_composition[plot_order, :], aspect="auto", cmap="RdBu_r")
plt.xticks(
    range(len(clustering_features)),
    [x.replace("dpct_", "") for x in clustering_features],
    rotation=45,
    ha="right",
)
plt.yticks([])
plt.colorbar(label="Z score")
plt.tight_layout()
plt.savefig(output_dir / "Cohort4_TLS_composition_heatmap.pdf")
plt.close()
