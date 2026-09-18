###### cellular neighborhood analysis in TA-TLS ######

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import MiniBatchKMeans


###### read the final TA-TLS cells used in the analysis ######

cells = pd.read_csv(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\02-CN analysis\01-input\tls_delaunay_layer_cells_final.csv",
    usecols=[
        "cellid", "slide", "fov", "centroid.0", "centroid.1",
        "tls_id", "tls_subtype", "celltype_l2",
    ],
    dtype={
        "cellid": str,
        "slide": str,
        "fov": str,
        "tls_id": str,
        "tls_subtype": str,
        "celltype_l2": str,
    },
)

cells["fov_region"] = cells["fov"].str.extract(
    r"(?:-|_)(L|T|N)(?:-|_)",
    expand=False,
)


cells["centroid.0"] = pd.to_numeric(cells["centroid.0"], errors="coerce")
cells["centroid.1"] = pd.to_numeric(cells["centroid.1"], errors="coerce")
cells = cells.dropna(subset=["centroid.0", "centroid.1"]).copy()

print(cells.shape)
print(cells["tls_id"].nunique())
print(pd.crosstab(cells["fov_region"], columns="n_cells"))
print(cells["tls_subtype"].value_counts())


###### cell types used to describe each neighborhood ######

celltype_order = [
    "B",
    "CD4T",
    "CD8T",
    "T-B",
    "Cholangiocyte",
    "DC",
    "Endothelial",
    "Fibroblast",
    "HEV",
    "Lymphatic",
    "Macrophage",
    "Neutrophil",
    "Unassigned",
]

cells = cells[cells["celltype_l2"].isin(celltype_order)].copy()
cells["celltype_code"] = pd.Categorical(
    cells["celltype_l2"],
    categories=celltype_order,
).codes

print(cells["celltype_l2"].value_counts())


###### neighbor composition within a 100-pixel radius ######

all_neighborhood = []

for tls_id, one_tls in cells.groupby("tls_id", sort=False):

    one_tls = one_tls.reset_index(drop=True).copy()
    xy = one_tls[["centroid.0", "centroid.1"]].to_numpy(dtype=float)
    celltype_code = one_tls["celltype_code"].to_numpy(dtype=int)

    tree = cKDTree(xy)
    neighbor_index = tree.query_ball_point(xy, r=100)

    neighbor_count = np.zeros(
        (one_tls.shape[0], len(celltype_order)),
        dtype=np.int32,
    )

    for i, index in enumerate(neighbor_index):

        index = np.asarray(index, dtype=int)
        index = index[index != i]

        if index.shape[0] == 0:
            continue

        neighbor_count[i, :] = np.bincount(
            celltype_code[index],
            minlength=len(celltype_order),
        )

    n_neighbors = neighbor_count.sum(axis=1)
    has_neighbor = n_neighbors > 0

    if has_neighbor.sum() == 0:
        continue

    neighbor_proportion = (
        neighbor_count[has_neighbor, :]
        / n_neighbors[has_neighbor, None]
    )

    one_result = one_tls.loc[
        has_neighbor,
        [
            "cellid", "slide", "fov", "fov_region", "tls_id",
            "tls_subtype", "celltype_l2",
        ],
    ].reset_index(drop=True)

    one_result["n_neighbors"] = n_neighbors[has_neighbor]

    one_result = pd.concat(
        [
            one_result,
            pd.DataFrame(neighbor_proportion, columns=celltype_order),
        ],
        axis=1,
    )

    all_neighborhood.append(one_result)


neighborhood = pd.concat(all_neighborhood, ignore_index=True)

print(neighborhood.shape)
print(neighborhood["tls_id"].nunique())
print(cells.shape[0] - neighborhood.shape[0])
print(neighborhood["n_neighbors"].describe())


###### locked K = 9 CN clustering ######

kmeans = MiniBatchKMeans(
    n_clusters=9,
    random_state=20260618,
    n_init=20,
    batch_size=20000,
)

neighborhood["CN"] = (
    kmeans.fit_predict(neighborhood[celltype_order].to_numpy(dtype=float)) + 1
)
neighborhood["CN"] = "CN" + neighborhood["CN"].astype(str)

print(neighborhood["CN"].value_counts().sort_index())


###### inspect mean composition and annotate the nine CNs ######

cn_mean = (
    neighborhood
    .groupby("CN")[celltype_order]
    .mean()
    .sort_index()
)

print(cn_mean.round(3))

cn_name = {
    "CN1": "CN1_LEC_HEV",
    "CN2": "CN2_B",
    "CN3": "CN3_Neutrophil_Chol",
    "CN4": "CN4_CD8T",
    "CN5": "CN5_CD4T",
    "CN6": "CN6_B_CD4T",
    "CN7": "CN7_Macrophage",
    "CN8": "CN8_CD4T_CD8T",
    "CN9": "CN9_Fibro_Endo",
}

neighborhood["CN_name"] = neighborhood["CN"].map(cn_name)
cn_mean.insert(0, "CN_name", cn_mean.index.map(cn_name))

print(neighborhood["CN_name"].value_counts().sort_index())

###### save the cell assignment and the mean CN composition ######

output_path = Path(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\02-CN analysis\03-result\01-tables"
)

neighborhood[
    [
        "cellid", "slide", "fov", "fov_region", "tls_id",
        "tls_subtype", "celltype_l2", "n_neighbors", "CN", "CN_name",
    ]
].to_csv(
    output_path / "CN_cell_assignment_K9_final_TA_TLS.csv.gz",
    index=False,
    compression="gzip",
)

cn_mean.to_csv(
    output_path / "CN_mean_composition_K9_final_TA_TLS.csv",
)
