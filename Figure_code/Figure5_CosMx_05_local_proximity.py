###### Local proximity enrichment across TLS layers ######

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.spatial import cKDTree


###### files ######

all_cell_file = Path(
    r"F:\TLS DATA\07-Cosmx\02-analysis-2"
    r"\06-updated_TLS_label_mechanism\06-QC_sensitivity\data"
    r"\01_all_cells_updated_labels_classifier_key_genes.parquet"
)
tls_cell_file = Path(
    r"F:\TLS DATA\07-Cosmx\02-analysis-2"
    r"\06-updated_TLS_label_mechanism\06-QC_sensitivity\data"
    r"\02_tls_cells_updated_TLS_labels_key_genes.parquet"
)
tls_id_file = Path(r"F:\TLS Material\06-Code\Figure5_CosMx_152_TLS.csv")
result_dir = Path(
    r"F:\TLS Material\06-Code\result\Figure5_CosMx\05_FigureS10E"
)
result_dir.mkdir(parents=True, exist_ok=True)


###### cells and definitions ######

tls_ids = pd.read_csv(tls_id_file)["tls_id"]
all_cells = pd.read_parquet(
    all_cell_file,
    columns=[
        "cell_barcode", "sample", "x_um", "y_um",
        "updated_celltype", "CCL5",
    ],
)
tls_cells = pd.read_parquet(
    tls_cell_file,
    columns=[
        "cell_barcode", "sample", "x_um", "y_um", "tls_id",
        "updated_celltype", "delaunay_layer", "CCR1", "CXCR2",
    ],
)
inside = tls_cells.loc[tls_cells["tls_id"].isin(tls_ids)].copy()

axes = {
    "CCL5-high T / CCR1+ Neu": {
        "source": all_cells["updated_celltype"].eq("T cell") & all_cells["CCL5"].ge(3),
        "receiver": inside["updated_celltype"].eq("Neutrophil") & inside["CCR1"].ge(1),
    },
    "Chol-like / CXCR2+ Neu": {
        "source": all_cells["updated_celltype"].eq("Cholangiocyte"),
        "receiver": inside["updated_celltype"].eq("Neutrophil") & inside["CXCR2"].ge(1),
    },
}


###### nearest source in the same section ######

receiver_tables = []

for axis, definition in axes.items():
    sources = all_cells.loc[
        definition["source"],
        ["sample", "x_um", "y_um"],
    ]
    receivers = inside.loc[
        definition["receiver"],
        [
            "cell_barcode", "sample", "tls_id",
            "delaunay_layer", "x_um", "y_um",
        ],
    ].copy()
    receivers["nearest_source_distance_um"] = np.inf

    for sample, receiver_index in receivers.groupby("sample").groups.items():
        source_xy = sources.loc[
            sources["sample"].eq(sample), ["x_um", "y_um"]
        ].to_numpy()
        if len(source_xy) == 0:
            continue

        receiver_xy = receivers.loc[receiver_index, ["x_um", "y_um"]].to_numpy()
        distance, _ = cKDTree(source_xy).query(receiver_xy, k=1)
        receivers.loc[receiver_index, "nearest_source_distance_um"] = distance

    receivers["within_100um"] = receivers["nearest_source_distance_um"].le(100)
    receivers["axis"] = axis
    receiver_tables.append(receivers)

receiver_data = pd.concat(receiver_tables, ignore_index=True)
receiver_data.to_csv(
    result_dir / "FigureS10E_receiver_level_proximity.csv", index=False
)


###### pooled relative observed/expected enrichment ######

layers = ["Periphery", "Transition", "Core"]
rows = []

for axis, one in receiver_data.groupby("axis", sort=False):
    baseline_n = len(one)
    baseline_contact = int(one["within_100um"].sum())
    baseline_fraction = baseline_contact / baseline_n

    for layer in layers:
        layer_cells = one.loc[one["delaunay_layer"].eq(layer)]
        n_receptor = len(layer_cells)
        n_contact = int(layer_cells["within_100um"].sum())
        layer_fraction = n_contact / n_receptor

        rows.append(
            {
                "axis": axis,
                "layer": layer,
                "n_receptor": n_receptor,
                "n_contact": n_contact,
                "baseline_receptors": baseline_n,
                "baseline_contacts": baseline_contact,
                "ROE": layer_fraction / baseline_fraction,
            }
        )

roe = pd.DataFrame(rows)
roe.to_csv(result_dir / "FigureS10E_local_proximity_ROE.csv", index=False)

expected_baselines = {
    "CCL5-high T / CCR1+ Neu": (681, 445),
    "Chol-like / CXCR2+ Neu": (193, 88),
}
for axis, expected in expected_baselines.items():
    one = roe.loc[roe["axis"].eq(axis)].iloc[0]
    assert (one["baseline_receptors"], one["baseline_contacts"]) == expected


###### heatmap ######

matrix = (
    roe.pivot(index="layer", columns="axis", values="ROE")
    .reindex(index=layers, columns=list(axes))
)

fig, ax = plt.subplots(figsize=(3.4, 4.2))
image = ax.imshow(
    matrix.to_numpy(),
    cmap="RdBu_r",
    norm=TwoSlopeNorm(vmin=0.75, vcenter=1, vmax=1.25),
    aspect="equal",
)

ax.set_xticks(range(2), matrix.columns, rotation=90)
ax.set_yticks(range(3), matrix.index)
ax.yaxis.tick_right()
ax.tick_params(length=0)
ax.set_title("Local proximity\nenrichment", loc="left")

for edge in np.arange(-0.5, 3, 1):
    ax.axhline(edge, color="white", linewidth=1)
for edge in np.arange(-0.5, 2, 1):
    ax.axvline(edge, color="white", linewidth=1)

colorbar = fig.colorbar(
    image,
    ax=ax,
    orientation="horizontal",
    fraction=0.08,
    pad=0.42,
)
colorbar.set_label("ROE")
colorbar.set_ticks([0.75, 1.25], labels=["Low", "High"])
colorbar.ax.invert_xaxis()

fig.tight_layout()
fig.savefig(result_dir / "FigureS10E_local_proximity_ROE.pdf", bbox_inches="tight")
fig.savefig(
    result_dir / "FigureS10E_local_proximity_ROE.png",
    dpi=400,
    bbox_inches="tight",
)
plt.close(fig)
