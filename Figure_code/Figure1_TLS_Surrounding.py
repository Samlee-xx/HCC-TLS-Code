###### Cell composition and enrichment analysis of TA-TLSs ######

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata, wilcoxon
from statsmodels.stats.multitest import multipletests


###### paths and parameters ######

input_dir = Path(
    r"F:\TLS DATA\01-CODEX\03-Analysis-final\04-inside-outside TLS\01-input"
)
output_dir = Path(r"F:\TLS Material\06-Code\Figure1I_source_data")
output_dir.mkdir(parents=True, exist_ok=True)

cutoffs = [100, 200, 300, 400, 500]
main_cutoff = 200
full_fov_tls_fraction = 0.95

celltype_order = [
    "B", "CD4T", "CD8T", "T-B", "DC", "Neutrophil",
    "Macrophage", "Endothelial", "HEV", "Lymphatic",
    "Fibroblast", "Cholangiocyte", "Tumor", "Unassigned"
]


###### select the 667 clustered TA-TLSs ######

tls_info = pd.read_csv(
    input_dir / "tls_subtype_assignment_final.csv",
    dtype={"tls_id": str, "slide": str, "fov": str}
)

tls_ids = set(tls_info["tls_id"])

tls_cells = pd.read_csv(
    input_dir / "dbscan_tls_cells_final.csv.gz",
    usecols=["cellid", "slide", "fov", "tls_id"],
    dtype={"cellid": str, "slide": str, "fov": str, "tls_id": str}
)

tls_cells = tls_cells[tls_cells["tls_id"].isin(tls_ids)].copy()
tls_cells["fov_key"] = tls_cells["slide"] + "__" + tls_cells["fov"]

tls_cellids = set(tls_cells["cellid"])
tls_fovs = set(tls_cells["fov_key"])

print(tls_cells["tls_id"].nunique())
print(len(tls_fovs))


###### read cells from the selected FOVs ######

cell_list = []

for one_chunk in pd.read_csv(
    input_dir / "final_codex_celltyped_obs.csv",
    usecols=[
        "cellid", "slide", "fov", "centroid-0", "centroid-1",
        "celltype_l2"
    ],
    dtype={
        "cellid": str,
        "slide": str,
        "fov": str,
        "celltype_l2": str
    },
    chunksize=300000,
    low_memory=False
):
    one_chunk["fov_key"] = one_chunk["slide"] + "__" + one_chunk["fov"]
    one_chunk = one_chunk[one_chunk["fov_key"].isin(tls_fovs)]
    cell_list.append(one_chunk)

cells = pd.concat(cell_list, ignore_index=True)

cells = cells.rename(columns={
    "centroid-0": "x",
    "centroid-1": "y",
    "celltype_l2": "celltype"
})

cells["x"] = pd.to_numeric(cells["x"], errors="coerce")
cells["y"] = pd.to_numeric(cells["y"], errors="coerce")
cells = cells.dropna(subset=["x", "y"]).copy()
cells["celltype"] = cells["celltype"].replace({"Unassign": "Unassigned"})
cells["is_tls"] = cells["cellid"].isin(tls_cellids)


fov_coverage = (
    cells.groupby(["slide", "fov", "fov_key"], as_index=False)
    .agg(n_all_cells=("cellid", "size"), n_tls_cells=("is_tls", "sum"))
)

fov_coverage["tls_cell_fraction"] = (
    fov_coverage["n_tls_cells"] / fov_coverage["n_all_cells"]
)
fov_coverage["excluded"] = (
    fov_coverage["tls_cell_fraction"] >= full_fov_tls_fraction
)

excluded_fovs = set(fov_coverage.loc[fov_coverage["excluded"], "fov_key"])
cells = cells[~cells["fov_key"].isin(excluded_fovs)].copy()
tls_cells = tls_cells[~tls_cells["fov_key"].isin(excluded_fovs)].copy()

fov_coverage.to_csv(
    output_dir / "Figure1I_FOV_inclusion.csv",
    index=False
)

##### distance from each non-TLS cell to the nearest TLS cell ###

cells["nearest_tls_distance"] = np.nan

for fov_key, one_fov in cells.groupby("fov_key"):
    tls_index = one_fov.index[one_fov["is_tls"]]
    other_index = one_fov.index[~one_fov["is_tls"]]

    tls_xy = cells.loc[tls_index, ["x", "y"]].to_numpy()
    other_xy = cells.loc[other_index, ["x", "y"]].to_numpy()

    tree = cKDTree(tls_xy)
    distance, _ = tree.query(other_xy, k=1)
    cells.loc[other_index, "nearest_tls_distance"] = distance


### percentage, enrichment and cutoff sensitivity ###

all_percentage = []
all_percentage_stats = []
all_enrichment = []
all_enrichment_stats = []

fov_table = cells[["slide", "fov", "fov_key"]].drop_duplicates()

for cutoff in cutoffs:
    cells["region"] = "Other"
    cells.loc[cells["is_tls"], "region"] = "TLS"
    cells.loc[
        (~cells["is_tls"]) & (cells["nearest_tls_distance"] <= cutoff),
        "region"
    ] = "TLS-surrounding"

    local_cells = cells[cells["region"].isin(["TLS", "TLS-surrounding"])].copy()
    local_cells = local_cells[local_cells["celltype"] != "Artifact"]
    local_cells = local_cells[local_cells["celltype"].isin(celltype_order)]

    celltype_count = (
        local_cells.groupby(["slide", "fov", "fov_key", "region", "celltype"])
        .size()
        .reset_index(name="n_cells")
    )

    region_total = (
        local_cells.groupby(["slide", "fov", "fov_key", "region"])
        .size()
        .reset_index(name="n_region_cells")
    )

    complete_rows = []
    for _, one_fov in fov_table.iterrows():
        for region in ["TLS", "TLS-surrounding"]:
            for celltype in celltype_order:
                complete_rows.append({
                    "slide": one_fov["slide"],
                    "fov": one_fov["fov"],
                    "fov_key": one_fov["fov_key"],
                    "region": region,
                    "celltype": celltype
                })

    percentage = pd.DataFrame(complete_rows)
    percentage = percentage.merge(
        celltype_count,
        on=["slide", "fov", "fov_key", "region", "celltype"],
        how="left"
    )
    percentage = percentage.merge(
        region_total,
        on=["slide", "fov", "fov_key", "region"],
        how="left"
    )

    percentage["n_cells"] = percentage["n_cells"].fillna(0).astype(int)
    percentage["n_region_cells"] = percentage["n_region_cells"].fillna(0).astype(int)
    percentage["percentage"] = np.where(
        percentage["n_region_cells"] > 0,
        percentage["n_cells"] / percentage["n_region_cells"] * 100,
        np.nan
    )
    percentage["cutoff_px"] = cutoff
    all_percentage.append(percentage)

    percentage_stat = []

    for celltype in celltype_order:
        one = percentage[percentage["celltype"] == celltype]
        wide = one.pivot(
            index="fov_key",
            columns="region",
            values="percentage"
        ).dropna(subset=["TLS", "TLS-surrounding"])

        difference = wide["TLS"] - wide["TLS-surrounding"]
        nonzero_difference = difference[difference != 0].to_numpy()

        if len(nonzero_difference) == 0:
            p_value = 1.0
            r_rb = 0.0
        else:
            p_value = wilcoxon(difference, alternative="two-sided").pvalue
            ranks = rankdata(abs(nonzero_difference), method="average")
            w_positive = ranks[nonzero_difference > 0].sum()
            w_negative = ranks[nonzero_difference < 0].sum()
            r_rb = (w_positive - w_negative) / (w_positive + w_negative)

        percentage_stat.append({
            "cutoff_px": cutoff,
            "celltype": celltype,
            "n_fov_pairs": len(wide),
            "TLS_mean_percentage": wide["TLS"].mean(),
            "surrounding_mean_percentage": wide["TLS-surrounding"].mean(),
            "TLS_median_percentage": wide["TLS"].median(),
            "surrounding_median_percentage": wide["TLS-surrounding"].median(),
            "p_value": p_value,
            "r_rb": r_rb
        })

    percentage_stat = pd.DataFrame(percentage_stat)
    percentage_stat["p_adj_BH"] = multipletests(
        percentage_stat["p_value"], method="fdr_bh"
    )[1]
    all_percentage_stats.append(percentage_stat)

    count_wide = percentage.pivot_table(
        index=["slide", "fov", "fov_key", "celltype"],
        columns="region",
        values="n_cells",
        aggfunc="first"
    ).reset_index().fillna(0)

    count_wide = count_wide.rename(columns={
        "TLS": "n_tls_celltype",
        "TLS-surrounding": "n_surrounding_celltype"
    })

    total_wide = region_total.pivot_table(
        index=["slide", "fov", "fov_key"],
        columns="region",
        values="n_region_cells",
        aggfunc="first"
    ).reset_index().fillna(0)

    total_wide = total_wide.rename(columns={
        "TLS": "n_tls_total",
        "TLS-surrounding": "n_surrounding_total"
    })

    enrichment = count_wide.merge(
        total_wide,
        on=["slide", "fov", "fov_key"],
        how="left"
    )

    enrichment["n_celltype_total"] = (
        enrichment["n_tls_celltype"] + enrichment["n_surrounding_celltype"]
    )
    enrichment["n_local_total"] = (
        enrichment["n_tls_total"] + enrichment["n_surrounding_total"]
    )

    enrichment["observed_tls_fraction"] = (
        enrichment["n_tls_celltype"] / enrichment["n_celltype_total"]
    )
    enrichment["expected_tls_fraction"] = (
        enrichment["n_tls_total"] / enrichment["n_local_total"]
    )
    enrichment["enrichment_score"] = (
        enrichment["observed_tls_fraction"] /
        enrichment["expected_tls_fraction"]
    )

    enrichment = enrichment[
        (enrichment["n_celltype_total"] > 50) &
        (enrichment["n_surrounding_total"] > 0)
    ].copy()
    enrichment["cutoff_px"] = cutoff
    all_enrichment.append(enrichment)

    enrichment_stat = []

    for celltype in celltype_order:
        one = enrichment.loc[
            enrichment["celltype"] == celltype,
            "enrichment_score"
        ].dropna()

        if len(one) == 0:
            p_value = np.nan
        elif np.allclose(one, 1):
            p_value = 1.0
        else:
            p_value = wilcoxon(one - 1, alternative="two-sided").pvalue

        enrichment_stat.append({
            "cutoff_px": cutoff,
            "celltype": celltype,
            "n_fov": len(one),
            "mean_enrichment": one.mean(),
            "median_enrichment": one.median(),
            "p_value_against_1": p_value
        })

    enrichment_stat = pd.DataFrame(enrichment_stat)
    enrichment_stat["p_adj_BH"] = multipletests(
        enrichment_stat["p_value_against_1"].fillna(1),
        method="fdr_bh"
    )[1]
    all_enrichment_stats.append(enrichment_stat)

