###### cell-cell proximity analysis ######

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


###### TA-TLS cells ######

tls_cells = pd.read_csv(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\02-CN analysis\01-input\tls_delaunay_layer_cells_final.csv",
    usecols=[
        "tls_id", "tls_subtype", "centroid.0", "centroid.1", "celltype_l2",
    ],
    dtype={
        "tls_id": str,
        "tls_subtype": str,
        "celltype_l2": str,
    },
)

tls_cells["centroid.0"] = pd.to_numeric(
    tls_cells["centroid.0"], errors="coerce"
)
tls_cells["centroid.1"] = pd.to_numeric(
    tls_cells["centroid.1"], errors="coerce"
)
tls_cells = tls_cells.dropna(
    subset=["centroid.0", "centroid.1", "celltype_l2"]
).copy()

print(tls_cells.shape)
print(tls_cells["tls_id"].nunique())
print(tls_cells["tls_subtype"].value_counts())


###### cell types ######

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

tls_cells = tls_cells[tls_cells["celltype_l2"].isin(celltype_order)].copy()
tls_cells["celltype_code"] = pd.Categorical(
    tls_cells["celltype_l2"],
    categories=celltype_order,
).codes

print(tls_cells["celltype_l2"].value_counts())


###### proximity score in each TLS ######

np.random.seed(123)

pair_rows = []
chol_neu_rows = []

for tls_id, one_tls in tls_cells.groupby("tls_id", sort=True):

    one_tls = one_tls.reset_index(drop=True).copy()
    n_cells = one_tls.shape[0]

    if n_cells < 5: # to improve robust of our analysis, which mentioned in our methods
        continue

    tls_subtype = one_tls.loc[0, "tls_subtype"]
    xy = one_tls[["centroid.0", "centroid.1"]].to_numpy(dtype=float)
    celltype_code = one_tls["celltype_code"].to_numpy(dtype=int)

    tree = cKDTree(xy)
    neighbor_index = tree.query_ball_point(xy, r=100)

    neighbor_count = np.zeros(
        (n_cells, len(celltype_order)),
        dtype=np.int32,
    )

    for i, index in enumerate(neighbor_index):
        neighbor_count[i, :] = np.bincount(
            celltype_code[np.asarray(index, dtype=int)],
            minlength=len(celltype_order),
        )

    sample_size = int(round(n_cells * 0.20))
    sample_size = max(2, min(sample_size, n_cells))

    for center_celltype in celltype_order:

        center_code = celltype_order.index(center_celltype)
        center_index = np.where(celltype_code == center_code)[0]

        if center_index.shape[0] == 0:
            continue

        background_index = np.random.choice(
            n_cells,
            size=sample_size,
            replace=False,
        )

        background_count = neighbor_count[background_index, :]
        background_mean = background_count.mean(axis=0)
        background_sd = background_count.std(axis=0, ddof=1)
        observed_mean = neighbor_count[center_index, :].mean(axis=0)

        for neighbor_celltype in celltype_order:

            neighbor_code = celltype_order.index(neighbor_celltype)
            one_sd = background_sd[neighbor_code]

            if not np.isfinite(one_sd) or one_sd == 0:
                z_score = np.nan
            else:
                z_score = (
                    observed_mean[neighbor_code]
                    - background_mean[neighbor_code]
                ) / one_sd

            pair_rows.append(
                {
                    "tls_id": tls_id,
                    "tls_subtype": tls_subtype,
                    "center": center_celltype,
                    "neighbor": neighbor_celltype,
                    "z_score": z_score,
                    "observed_mean": observed_mean[neighbor_code],
                    "background_mean": background_mean[neighbor_code],
                    "background_sd": one_sd,
                    "center_n_cells": center_index.shape[0],
                    "tls_n_cells": n_cells,
                    "background_n_cells": sample_size,
                }
            )

    chol_code = celltype_order.index("Cholangiocyte")
    neu_code = celltype_order.index("Neutrophil")
    chol_index = np.where(celltype_code == chol_code)[0]
    neu_index = np.where(celltype_code == neu_code)[0]

    if chol_index.shape[0] > 0:

        neu_near_chol = set()

        for i in chol_index:
            one_neu = [
                j for j in neighbor_index[i]
                if celltype_code[j] == neu_code
            ]
            neu_near_chol.update(one_neu)

        n_neu_near_chol = len(neu_near_chol)

        chol_neu_rows.append(
            {
                "tls_id": tls_id,
                "tls_subtype": tls_subtype,
                "tls_n_cells": n_cells,
                "n_chol": chol_index.shape[0],
                "n_neu": neu_index.shape[0],
                "n_neu_within_100px_of_chol": n_neu_near_chol,
                "percentage_of_all_cells": 100 * n_neu_near_chol / n_cells,
                "percentage_of_all_neutrophils": (
                    100 * n_neu_near_chol / neu_index.shape[0]
                    if neu_index.shape[0] > 0
                    else np.nan
                ),
            }
        )


pair_by_tls = pd.DataFrame(pair_rows)
chol_neu_by_tls = pd.DataFrame(chol_neu_rows)

print(pair_by_tls.shape)
print(pair_by_tls["tls_id"].nunique())
print(chol_neu_by_tls.shape)
print(pd.crosstab(chol_neu_by_tls["tls_subtype"], columns="n_TLS"))


###### average TLS-level scores for the network ######

pair_summary = (
    pair_by_tls
    .groupby(["center", "neighbor"], as_index=False)
    .agg(
        mean_z_score=("z_score", "mean"),
        n_TLS=("z_score", "count"),
    )
)

print(
    pair_summary[
        (pair_summary["neighbor"] == "Cholangiocyte")
        & (pair_summary["center"] != "Cholangiocyte")
    ][["center", "mean_z_score", "n_TLS"]]
)


###### save the two tables used for plotting ######

output_path = Path(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\03-BileDuct analysis\03-result\01-tables"
)
output_path.mkdir(parents=True, exist_ok=True)

pair_by_tls.to_csv(
    output_path / "cell_cell_proximity_by_TLS_simple.csv",
    index=False,
)

chol_neu_by_tls.to_csv(
    output_path / "chol_neu_proximity_by_TLS_simple.csv",
    index=False,
)
