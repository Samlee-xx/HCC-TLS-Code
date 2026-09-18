###### TLS boundary detection by DBSCAN ######

from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.colors import to_hex
from sklearn.cluster import DBSCAN

import ark.settings as settings
from ark.utils.plot_utils import cohort_cluster_plot


###### read the final cell type annotation ######

cells = pd.read_csv(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\01-TLS subtype\01-DBSCAN\01-input\final_codex_celltyped_obs.csv",
    usecols=[
        "cellid", "slide", "fov", "label", "centroid-0", "centroid-1",
        "celltype_l2",
    ],
    dtype={"cellid": str, "slide": str, "fov": str},
)

cells = cells.rename(
    columns={
        "centroid-0": "centroid.0",
        "centroid-1": "centroid.1",
    }
)

cells["label"] = pd.to_numeric(cells["label"], errors="coerce")
cells["centroid.0"] = pd.to_numeric(cells["centroid.0"], errors="coerce")
cells["centroid.1"] = pd.to_numeric(cells["centroid.1"], errors="coerce")

print(cells.shape)
print(cells["celltype_l2"].value_counts())


###### cells used for DBSCAN ######

lymphoid = ["B", "CD4T", "CD8T", "T-B"]
myeloid = ["Macrophage", "DC", "Neutrophil"]
stroma = ["Fibroblast", "Endothelial", "HEV", "Lymphatic"]
other = ["Cholangiocyte", "Unassigned"]

tls_celltypes = lymphoid + myeloid + stroma + other

dbscan_cells = cells[cells["celltype_l2"].isin(tls_celltypes)].copy()

# 排除掉背景的FOV，这些是前期探索纳入的背景，后期没有用此进行分析，故排除
dbscan_cells = dbscan_cells[~dbscan_cells["fov"].str.endswith("B", na=False)].copy()
dbscan_cells = dbscan_cells.dropna(
    subset=["label", "centroid.0", "centroid.1"]
).copy()

dbscan_cells["label"] = dbscan_cells["label"].astype(int)
dbscan_cells["dbscan_label"] = -1

print(dbscan_cells.shape)
print(dbscan_cells["celltype_l2"].value_counts())


###### run DBSCAN in each FOV ######

for (slide, fov), one_fov in dbscan_cells.groupby(["slide", "fov"], sort=False):

    if one_fov.shape[0] < 50:
        continue

    xy = one_fov[["centroid.0", "centroid.1"]].to_numpy()
    label = DBSCAN(eps=80, min_samples=50).fit_predict(xy)
    dbscan_cells.loc[one_fov.index, "dbscan_label"] = label


###### name candidate TLS regions ######

dbscan_cells["tls_id"] = "Noise"
dbscan_cells["tls_number"] = pd.NA

clustered = dbscan_cells[dbscan_cells["dbscan_label"] >= 0]

for (slide, fov), one_fov in clustered.groupby(["slide", "fov"], sort=False):

    cluster_order = one_fov["dbscan_label"].value_counts().index.tolist()

    for tls_number, cluster_number in enumerate(cluster_order, start=1):
        index = one_fov.index[one_fov["dbscan_label"] == cluster_number]
        dbscan_cells.loc[index, "tls_number"] = tls_number
        dbscan_cells.loc[index, "tls_id"] = f"{slide}_{fov}_TLS{tls_number}"

dbscan_cells["tls_number"] = dbscan_cells["tls_number"].astype("Int64")
candidate_cells = dbscan_cells[dbscan_cells["dbscan_label"] >= 0].copy()

print(candidate_cells["tls_id"].nunique())
print(candidate_cells.groupby("tls_id").size().describe())


###### candidate TLS composition ######

tls_summary = (
    candidate_cells
    .groupby(["tls_id", "slide", "fov", "tls_number", "dbscan_label"])
    .size()
    .reset_index(name="n_cells")
)

celltype_count = pd.crosstab(
    candidate_cells["tls_id"],
    candidate_cells["celltype_l2"],
).reset_index()

tls_summary = tls_summary.merge(celltype_count, on="tls_id", how="left")

for celltype in tls_celltypes:
    if celltype not in tls_summary.columns:
        tls_summary[celltype] = 0

tls_summary["n_lymphoid"] = tls_summary[lymphoid].sum(axis=1)
tls_summary["n_stroma"] = tls_summary[stroma].sum(axis=1)

tls_summary["pct_Unassigned"] = (
    tls_summary["Unassigned"] / tls_summary["n_cells"] * 100
)
tls_summary["pct_lymphoid"] = (
    tls_summary["n_lymphoid"] / tls_summary["n_cells"] * 100
)
tls_summary["pct_stroma"] = (
    tls_summary["n_stroma"] / tls_summary["n_cells"] * 100
)

tls_summary["keep"] = (
    (tls_summary["n_cells"] > 50)
    & (tls_summary["pct_Unassigned"] <= 30)
    & (tls_summary["pct_stroma"] <= 60)
    & (tls_summary["pct_lymphoid"] >= 10)
)

print(tls_summary["keep"].value_counts())

keep_id = tls_summary.loc[tls_summary["keep"], "tls_id"]
tls_summary_keep = tls_summary[tls_summary["keep"]].copy()
tls_cells_keep = candidate_cells[candidate_cells["tls_id"].isin(keep_id)].copy()

print(tls_summary_keep.shape)
print(tls_cells_keep.shape)


