###### Delaunay reconstruction of TLS layers #####

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, QhullError, cKDTree


###### final TLS cells used in the analysis ######

tls_cells = pd.read_csv(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\01-TLS subtype\03-TLS Delaunay\01-input\dbscan_tls_cells_final.csv.gz",
    usecols=[
        "cellid", "slide", "fov", "label", "centroid.0", "centroid.1",
        "tls_id", "tls_number", "celltype_l2",
    ],
    dtype={"cellid": str, "slide": str, "fov": str, "tls_id": str},
)

tls_subtype = pd.read_csv(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\01-TLS subtype\03-TLS Delaunay\01-input\tls_subtype_assignment_final.csv",
    usecols=["tls_id", "tls_subtype"],
    dtype={"tls_id": str},
)

# The final upstream TLS table is used here; no FOV name is removed directly.
tls_cells = tls_cells.merge(tls_subtype, on="tls_id", how="inner")

print(tls_cells.shape)
print(tls_cells["tls_id"].nunique())
print(tls_cells["tls_subtype"].value_counts())


###### Delaunay graph and layer expansion ######

all_layers = []
layer_summary = []

for tls_id, one_tls in tls_cells.groupby("tls_id", sort=False):

    one_tls = one_tls.copy().reset_index(drop=True)
    one_tls["delaunay_level"] = np.nan
    one_tls["delaunay_layer"] = pd.NA

    xy = one_tls[["centroid.0", "centroid.1"]].to_numpy(dtype=float)

    try:
        tri = Delaunay(xy)
    except QhullError:
        all_layers.append(one_tls)
        continue

    edge = np.vstack(
        [
            tri.simplices[:, [0, 1]],
            tri.simplices[:, [0, 2]],
            tri.simplices[:, [1, 2]],
        ]
    )

    edge = np.unique(np.sort(edge, axis=1), axis=0)
    distance = np.linalg.norm(xy[edge[:, 0]] - xy[edge[:, 1]], axis=1)
    edge = edge[distance < 40]

    if edge.shape[0] == 0:
        all_layers.append(one_tls)
        continue

    graph_nodes = np.unique(edge.ravel())
    graph = {int(node): set() for node in graph_nodes}

    for node_a, node_b in edge:
        graph[int(node_a)].add(int(node_b))
        graph[int(node_b)].add(int(node_a))

    tree = cKDTree(xy[graph_nodes])
    neighbor_number = np.array(
        [
            len(x) - 1
            for x in tree.query_ball_point(xy[graph_nodes], r=45)
        ]
    )

    boundary = graph_nodes[neighbor_number < 15]

    if boundary.shape[0] == 0:
        all_layers.append(one_tls)
        continue

    level_by_node = {int(node): 1 for node in boundary}
    visited = set(int(node) for node in boundary)
    frontier = set(int(node) for node in boundary)
    level = 1

    while len(frontier) > 0:

        next_frontier = set()

        for node in frontier:
            next_frontier.update(graph[node])

        next_frontier = next_frontier - visited

        if len(next_frontier) == 0:
            break

        level += 1

        for node in next_frontier:
            level_by_node[node] = level

        visited.update(next_frontier)
        frontier = next_frontier

    for node, node_level in level_by_node.items():
        one_tls.loc[node, "delaunay_level"] = node_level

    max_level = int(one_tls["delaunay_level"].max())

    if max_level <= 3:
        one_tls.loc[one_tls["delaunay_level"] == 1, "delaunay_layer"] = "Periphery"
        one_tls.loc[one_tls["delaunay_level"] == 2, "delaunay_layer"] = "Transition"
        one_tls.loc[one_tls["delaunay_level"] == 3, "delaunay_layer"] = "Core"
    else:
        one_third = max_level // 3

        one_tls.loc[
            one_tls["delaunay_level"] <= one_third,
            "delaunay_layer",
        ] = "Periphery"

        one_tls.loc[
            one_tls["delaunay_level"].between(one_third + 1, 2 * one_third),
            "delaunay_layer",
        ] = "Transition"

        one_tls.loc[
            one_tls["delaunay_level"] > 2 * one_third,
            "delaunay_layer",
        ] = "Core"

    all_layers.append(one_tls)

    layer_summary.append(
        {
            "tls_id": tls_id,
            "n_cells": one_tls.shape[0],
            "max_level": max_level,
            "n_labeled": one_tls["delaunay_layer"].notna().sum(),
        }
    )

tls_layers = pd.concat(all_layers, ignore_index=True)
layer_summary = pd.DataFrame(layer_summary)

print(tls_layers["delaunay_layer"].value_counts(dropna=False))
print(layer_summary["max_level"].describe())
print((layer_summary["n_cells"] - layer_summary["n_labeled"]).sum())


###### draw the three spatial layers ######

plot_path = Path(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\01-TLS subtype\03-TLS Delaunay\03-result\02-figures\tls_delaunay_layer_plots_simple"
)
plot_path.mkdir(parents=True, exist_ok=True)

layer_color = {
    "Periphery": "#377EB8",
    "Transition": "#FF7F00",
    "Core": "#E41A1C",
}

for tls_id, one_tls in tls_layers.groupby("tls_id", sort=False):

    fig, ax = plt.subplots(figsize=(4, 4))

    unlabeled = one_tls[one_tls["delaunay_layer"].isna()]
    ax.scatter(
        unlabeled["centroid.0"],
        unlabeled["centroid.1"],
        s=2,
        color="lightgrey",
    )

    for layer in ["Periphery", "Transition", "Core"]:
        one_layer = one_tls[one_tls["delaunay_layer"] == layer]
        ax.scatter(
            one_layer["centroid.0"],
            one_layer["centroid.1"],
            s=2,
            color=layer_color[layer],
            label=layer,
        )

    ax.set_title(tls_id)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend(frameon=False, markerscale=3)
    fig.tight_layout()
    fig.savefig(plot_path / f"{tls_id}.png", dpi=200)
    plt.close(fig)


