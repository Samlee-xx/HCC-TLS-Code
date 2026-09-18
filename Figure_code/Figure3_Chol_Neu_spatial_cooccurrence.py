from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


# Spatial co-occurrence of cholangiocyte-like cells and neutrophils


tls_file = Path(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\02-CN analysis\01-input\tls_delaunay_layer_cells_final.csv"
)
all_cell_file = Path(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\03-BileDuct analysis\01-input\final_codex_celltyped_obs.csv"
)
out_dir = Path(r"F:\TLS Material\06-Code\result\Figure2_Chol_Neu_spatial_cooccurrence")
out_dir.mkdir(parents=True, exist_ok=True)



tls_cell = pd.read_csv(tls_file, dtype={"slide": str, "fov": str})

print(tls_cell.shape)
print(tls_cell["tls_id"].nunique())
print(
    tls_cell[["tls_id", "tls_subtype"]]
    .drop_duplicates()["tls_subtype"]
    .value_counts()
)
print(tls_cell["celltype_l2"].value_counts())

chol = tls_cell.loc[
    tls_cell["celltype_l2"].eq("Cholangiocyte"),
    ["cellid", "slide", "fov", "tls_id", "tls_subtype", "centroid.0", "centroid.1"],
].copy()

chol["fov_key"] = chol["slide"] + "__" + chol["fov"]

print(chol.shape)
print(chol["tls_subtype"].value_counts())


# 2. Read all cells from the FOVs represented in the final TLS table
fov_keep = set(chol["fov_key"])
all_cell_part = []

for x in pd.read_csv(
    all_cell_file,
    usecols=[
        "cellid",
        "slide",
        "fov",
        "centroid-0",
        "centroid-1",
        "celltype_l2",
    ],
    dtype={"slide": str, "fov": str},
    chunksize=300000,
):
    x["fov_key"] = x["slide"] + "__" + x["fov"]
    x = x.loc[x["fov_key"].isin(fov_keep)]
    if len(x) > 0:
        all_cell_part.append(x)

all_cell = pd.concat(all_cell_part, ignore_index=True)

print(all_cell.shape)
print(all_cell["fov_key"].nunique())
print(all_cell["celltype_l2"].value_counts())


# 3. Calculate the neutrophil co-occurrence score at radii of 50-300 pixels
radii = np.arange(50, 301, 10)
score_raw = pd.DataFrame(
    np.nan,
    index=chol["cellid"],
    columns=[f"r{x}" for x in radii],
)
score_raw.index.name = "cellid"

chol_by_fov = chol.groupby("fov_key", sort=False)

for fov_key, fov_cell in all_cell.groupby("fov_key", sort=False):
    chol_fov = chol_by_fov.get_group(fov_key)

    xy_all = fov_cell[["centroid-0", "centroid-1"]].to_numpy(float)
    xy_chol = chol_fov[["centroid.0", "centroid.1"]].to_numpy(float)
    celltype = fov_cell["celltype_l2"].to_numpy()

    fov_neu_fraction = np.mean(celltype == "Neutrophil")
    if fov_neu_fraction == 0:
        continue

    tree = cKDTree(xy_all)

    for radius in radii:
        neighbor_index = tree.query_ball_point(xy_chol, r=radius)
        n_all = np.fromiter((len(i) for i in neighbor_index), dtype=int)
        n_neu = np.fromiter(
            (np.sum(celltype[i] == "Neutrophil") for i in neighbor_index),
            dtype=int,
        )

        local_neu_fraction = n_neu / n_all
        score_raw.loc[chol_fov["cellid"], f"r{radius}"] = (
            local_neu_fraction / fov_neu_fraction
        )


# 4. Retain complete, nonconstant radial profiles
complete = np.isfinite(score_raw).all(axis=1)
nonconstant = score_raw.std(axis=1) > 0
score_raw = score_raw.loc[complete & nonconstant].copy()

print(score_raw.shape)


# 5. log2 transform, standardize each radius across cells, and perform K-means
score_log = np.log2(score_raw + 0.05)
score_for_cluster = StandardScaler().fit_transform(score_log)

kmeans = KMeans(n_clusters=4, random_state=20260619, n_init=50)
raw_cluster = kmeans.fit_predict(score_for_cluster)


# 6. Name the two clusters with the higher short-range enrichment as NeuProx-Chols
# Here, short range is defined as 50-100 pixels, consistent with the 100-pixel
# neighborhood used in the preceding spatial analyses.
short_columns = [f"r{x}" for x in range(50, 101, 10)]
cluster_mean = score_raw.assign(raw_cluster=raw_cluster).groupby("raw_cluster").mean()
cluster_order = cluster_mean[short_columns].mean(axis=1).sort_values(ascending=False).index
cluster_name = {x: f"Cluster{i + 1}" for i, x in enumerate(cluster_order)}

assignment = chol.set_index("cellid").loc[score_raw.index].copy()
assignment["cluster"] = pd.Series(raw_cluster, index=score_raw.index).map(cluster_name)
assignment["chol_group"] = np.where(
    assignment["cluster"].isin(["Cluster1", "Cluster2"]),
    "NeuProx-Chols",
    "NeuDist-Chols",
)

cluster_mean.index = cluster_mean.index.map(cluster_name)
cluster_mean = cluster_mean.loc[["Cluster1", "Cluster2", "Cluster3", "Cluster4"]]
cluster_mean.index.name = "cluster"

print(assignment["cluster"].value_counts().sort_index())
print(pd.crosstab(assignment["cluster"], assignment["chol_group"]))
print(cluster_mean[short_columns].mean(axis=1))

# Row-wise Z scores are made only for displaying the radial heatmap.
score_row_z = score_raw.sub(score_raw.mean(axis=1), axis=0)
score_row_z = score_row_z.div(score_raw.std(axis=1), axis=0)


# 7. Calculate the NeuProx-Chol percentage in each TLS
tls_percentage = (
    assignment.groupby(["tls_id", "tls_subtype"])
    .size()
    .rename("n_clusterable_chol")
    .reset_index()
)

neu_prox_count = (
    assignment.loc[assignment["chol_group"].eq("NeuProx-Chols")]
    .groupby("tls_id")
    .size()
    .rename("n_neu_prox_chol")
    .reset_index()
)

tls_percentage = tls_percentage.merge(neu_prox_count, on="tls_id", how="left")
tls_percentage["n_neu_prox_chol"] = tls_percentage["n_neu_prox_chol"].fillna(0).astype(int)
tls_percentage["neu_prox_percentage"] = (
    100 * tls_percentage["n_neu_prox_chol"] / tls_percentage["n_clusterable_chol"]
)

# Methods: only TLSs with more than 10 eligible, clusterable cholangiocyte-like cells. To improve the robustness of our analyses.
tls_percentage = tls_percentage.loc[tls_percentage["n_clusterable_chol"] > 10].copy()

print(pd.crosstab(tls_percentage["tls_subtype"], columns="n_TLS"))
print(tls_percentage.groupby("tls_subtype")["neu_prox_percentage"].describe())

bf = tls_percentage.loc[
    tls_percentage["tls_subtype"].eq("BF-TLS"), "neu_prox_percentage"
]
tnc = tls_percentage.loc[
    tls_percentage["tls_subtype"].eq("TNC-TLS"), "neu_prox_percentage"
]

test = mannwhitneyu(tnc, bf, alternative="two-sided")
print("Mann-Whitney U =", test.statistic)
print("P =", test.pvalue)


# 8. Save the tables used for plotting
assignment.reset_index().to_csv(
    out_dir / "chol_neu_cell_assignment.csv.gz", index=False
)
score_raw.to_csv(out_dir / "chol_neu_score_raw.csv.gz")
score_row_z.to_csv(out_dir / "chol_neu_score_row_z.csv.gz")
cluster_mean.to_csv(out_dir / "chol_neu_cluster_mean_raw.csv")
tls_percentage.to_csv(out_dir / "chol_neu_TLS_percentage.csv", index=False)
