###### Receptor expression along the TLS boundary-to-core axis ######

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel

###### files ######

cell_file = Path(
    r"F:\TLS Material\06-Code\result\Figure5_CosMx\02_TLS_identification"
    r"\CosMx_TLS_cell_assignments_152.parquet"
)
tls_id_file = Path(r"F:\TLS Material\06-Code\Figure5_CosMx_152_TLS.csv")
result_dir = Path(
    r"F:\TLS Material\06-Code\result\Figure5_CosMx\04_Figure5T"
)
result_dir.mkdir(parents=True, exist_ok=True)


###### TLS-level receptor abundance ######

tls_ids = pd.read_csv(tls_id_file)["tls_id"]
columns = [
    "cell_barcode",
    "tls_id",
    "inside_tls",
    "updated_celltype",
    "signed_distance_um",
    "nCount_RNA",
    "CXCR2",
    "CCR1",
]
cells = pd.read_parquet(cell_file, columns=columns)
cells = cells.loc[
    cells["inside_tls"]
    & cells["tls_id"].isin(tls_ids)
    & cells["updated_celltype"].eq("Neutrophil")
].copy()

cells["distance_band"] = pd.cut(
    cells["signed_distance_um"],
    [-0.001, 50, 100, np.inf],
    right=False,
    labels=["0–50 µm", "50–100 µm", "≥100 µm"],
)

counts = (
    cells.groupby(["tls_id", "distance_band"], observed=True)
    .agg(
        n_neutrophils=("cell_barcode", "size"),
        total_neutrophil_RNA=("nCount_RNA", "sum"),
        CXCR2_RNA=("CXCR2", "sum"),
        CCR1_RNA=("CCR1", "sum"),
    )
    .reset_index()
)

grid = pd.MultiIndex.from_product(
    [tls_ids, ["0–50 µm", "50–100 µm", "≥100 µm"]],
    names=["tls_id", "distance_band"],
).to_frame(index=False)
counts = grid.merge(counts, how="left", on=["tls_id", "distance_band"])
counts[["n_neutrophils", "total_neutrophil_RNA", "CXCR2_RNA", "CCR1_RNA"]] = (
    counts[["n_neutrophils", "total_neutrophil_RNA", "CXCR2_RNA", "CCR1_RNA"]]
    .fillna(0)
)

for receptor in ["CXCR2", "CCR1"]:
    counts[f"{receptor}_RNA_per10k"] = (
        counts[f"{receptor}_RNA"]
        / counts["total_neutrophil_RNA"].replace(0, np.nan)
        * 10000
    )

band_coverage = counts.pivot(
    index="tls_id", columns="distance_band", values="n_neutrophils"
)
three_band_tls = band_coverage.index[(band_coverage > 0).all(axis=1)]


counts_116 = counts.loc[counts["tls_id"].isin(three_band_tls)].copy()
counts_116.to_csv(
    result_dir / "Figure5T_TLS_level_receptor_abundance.csv", index=False
)


###### mean profiles and paired tests ######

band_order = ["0–50 µm", "50–100 µm", "≥100 µm"]
profile_rows = []
test_rows = []

for receptor in ["CXCR2", "CCR1"]:
    value = f"{receptor}_RNA_per10k"
    wide = counts_116.pivot(index="tls_id", columns="distance_band", values=value)
    wide = wide[band_order]

    means = wide.mean(axis=0)
    relative = means / means.iloc[0]

    for order, band in enumerate(band_order, start=1):
        profile_rows.append(
            {
                "receptor": receptor,
                "distance_band": band,
                "band_order": order,
                "n_TLS": len(wide),
                "mean_RNA_per10k": means[band],
                "relative_RNA": relative[band],
            }
        )

    comparisons = [
        ("50–100 µm", "0–50 µm"),
        ("≥100 µm", "0–50 µm"),
        ("≥100 µm", "50–100 µm"),
    ]
    receptor_tests = []

    for band_1, band_2 in comparisons:
        test = ttest_rel(wide[band_1], wide[band_2])
        receptor_tests.append(
            {
                "receptor": receptor,
                "comparison": f"{band_1} vs {band_2}",
                "n_TLS": len(wide),
                "t": test.statistic,
                "p_raw": test.pvalue,
            }
        )

    raw_p = np.array([row["p_raw"] for row in receptor_tests])
    order = np.argsort(raw_p)
    adjusted = np.empty(3)
    ranked = raw_p[order] * 3 / np.arange(1, 4)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1)

    for row, p_adjusted in zip(receptor_tests, adjusted):
        row["p_BH_within_receptor"] = p_adjusted
        test_rows.append(row)

profiles = pd.DataFrame(profile_rows)
tests = pd.DataFrame(test_rows)
profiles.to_csv(result_dir / "Figure5T_mean_profiles.csv", index=False)
tests.to_csv(result_dir / "Figure5T_paired_tests.csv", index=False)


###### plot ######

fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0))

for ax, receptor, color in zip(
    axes,
    ["CXCR2", "CCR1"],
    ["#327CAD", "#BE3F50"],
):
    one = profiles.loc[profiles["receptor"].eq(receptor)].sort_values("band_order")
    endpoint = tests.loc[
        tests["receptor"].eq(receptor)
        & tests["comparison"].eq("≥100 µm vs 0–50 µm")
    ].iloc[0]

    ax.axhline(1, color="grey60", linestyle="--", linewidth=1)
    ax.plot(
        range(3), one["relative_RNA"],
        color=color, marker="o", linewidth=1.5, markersize=5
    )
    ax.set_xticks(range(3), band_order)
    ax.set_title(f"Neutrophil {receptor} RNA")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim((0.30, 1.20) if receptor == "CXCR2" else (0.90, 2.05))

    p_label = (
        f"{endpoint.p_raw:.6f}" if receptor == "CXCR2"
        else f"{endpoint.p_raw:.3f}"
    )
    ax.text(
        0.48, 0.80,
        f"≥100 vs 0–50 µm\n$P$ = {p_label}",
        transform=ax.transAxes,
    )

axes[0].set_ylabel("Relative RNA abundance")
fig.supxlabel("Distance inward from the TLS boundary")
fig.tight_layout()
fig.savefig(result_dir / "Figure5T_receptor_depth.pdf", bbox_inches="tight")
fig.savefig(result_dir / "Figure5T_receptor_depth.png", dpi=400, bbox_inches="tight")
plt.close(fig)
