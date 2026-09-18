###### CosMx TLS identification ######

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates
from skimage.feature import peak_local_max
from skimage.morphology import remove_small_objects
from skimage.segmentation import expand_labels, watershed


###### files ######

cell_file = Path(
    r"F:\TLS DATA\07-Cosmx\02-analysis-2"
    r"\06-updated_TLS_label_mechanism\06-QC_sensitivity\data"
    r"\01_all_cells_updated_labels_classifier_key_genes.parquet"
)
tls_id_file = Path(r"F:\TLS Material\06-Code\Figure5_CosMx_152_TLS.csv")
result_dir = Path(
    r"F:\TLS Material\06-Code\result\Figure5_CosMx\02_TLS_identification"
)


###### parameters ######

bin_um = 20.0
sigma_bins = 2.0
peak_quantile = 0.95
grow_quantile = 0.75
peak_min_distance_bins = 5
expand_distance_bins = 2

lymphoid_types = {"T cell", "B cell", "GC-B", "Plasma cell"}
b_lineage_types = {"B cell", "GC-B", "Plasma cell"}

min_total_cells = 50
min_lymphoid_cells = 20
min_b_lineage_cells = 10
min_lymphoid_percent = 15.0
min_area_um2 = 4000.0
max_area_um2 = 2_000_000.0


###### lymphoid-density watershed ######

def make_raster(one):
    x = one["x_um"].to_numpy()
    y = one["y_um"].to_numpy()

    x0 = np.floor(x.min() / bin_um) * bin_um - bin_um
    y0 = np.floor(y.min() / bin_um) * bin_um - bin_um
    nx = int(np.ceil((x.max() - x0) / bin_um)) + 2
    ny = int(np.ceil((y.max() - y0) / bin_um)) + 2

    ix = np.clip(((x - x0) // bin_um).astype(int), 0, nx - 1)
    iy = np.clip(((y - y0) // bin_um).astype(int), 0, ny - 1)

    total = np.zeros((ny, nx), dtype=np.float32)
    lymphoid = np.zeros_like(total)
    np.add.at(total, (iy, ix), 1)

    is_lymphoid = one["updated_celltype"].isin(lymphoid_types).to_numpy()
    np.add.at(lymphoid, (iy[is_lymphoid], ix[is_lymphoid]), 1)
    return x0, y0, ix, iy, total, lymphoid


def segment_tls(total, lymphoid):
    smooth_total = gaussian_filter(total, sigma_bins)
    smooth_lymphoid = gaussian_filter(lymphoid, sigma_bins)
    tissue = smooth_total > 0.20

    score = (
        smooth_lymphoid / (smooth_total + 0.10)
    ) * np.log1p(smooth_lymphoid)

    peak_cutoff = np.quantile(score[tissue], peak_quantile)
    grow_cutoff = np.quantile(score[tissue], grow_quantile)

    peaks = peak_local_max(
        score,
        min_distance=peak_min_distance_bins,
        threshold_abs=peak_cutoff,
        exclude_border=False,
        labels=tissue.astype(np.uint8),
    )

    markers = np.zeros_like(total, dtype=np.int32)
    for number, (row, column) in enumerate(peaks, start=1):
        markers[row, column] = number

    mask = tissue & (score >= grow_cutoff) & (smooth_lymphoid >= 0.05)
    mask = remove_small_objects(mask, min_size=4)
    labels = watershed(-score, markers=markers, mask=mask)
    labels = expand_labels(labels, distance=expand_distance_bins)
    labels[~tissue] = 0
    return labels, tissue


###### identify TLSs section by section ######

cells = pd.read_parquet(cell_file)
tls_ids = pd.read_csv(tls_id_file)["tls_id"]
tls_id_set = set(tls_ids)

assignment_tables = []
region_tables = []

for section_number, (sample, one) in enumerate(
    cells.groupby("sample", sort=True), start=1
):
    one = one.copy()
    x0, y0, ix, iy, total, lymphoid = make_raster(one)
    labels, tissue = segment_tls(total, lymphoid)
    one["region"] = labels[iy, ix]

    tissue_edge_distance = distance_transform_edt(tissue) * bin_um
    rows = []

    for region, region_cells in one.loc[one["region"] > 0].groupby("region"):
        region_mask = labels == region
        is_lymphoid = region_cells["updated_celltype"].isin(lymphoid_types)
        is_b_lineage = region_cells["updated_celltype"].isin(b_lineage_types)

        row = {
            "region": int(region),
            "n_cells": len(region_cells),
            "n_lymphoid": int(is_lymphoid.sum()),
            "n_B_lineage": int(is_b_lineage.sum()),
            "lymphoid_percent": 100 * is_lymphoid.mean(),
            "area_um2": int(region_mask.sum()) * bin_um**2,
            "centroid_x_um": region_cells["x_um"].mean(),
            "centroid_y_um": region_cells["y_um"].mean(),
            "max_interior_depth_um": distance_transform_edt(region_mask).max() * bin_um,
            "edge_fraction_lt40um": np.mean(tissue_edge_distance[region_mask] < 40),
        }
        rows.append(row)

    regions = pd.DataFrame(rows)
    regions["meets_tls_definition"] = (
        (regions["n_cells"] >= min_total_cells)
        & (regions["n_lymphoid"] >= min_lymphoid_cells)
        & (regions["n_B_lineage"] >= min_b_lineage_cells)
        & (regions["lymphoid_percent"] >= min_lymphoid_percent)
        & (regions["area_um2"] >= min_area_um2)
        & (regions["area_um2"] <= max_area_um2)
    )

    valid_regions = regions.loc[regions["meets_tls_definition"]].sort_values(
        ["centroid_y_um", "centroid_x_um"]
    )
    region_to_tls = {
        region: f"COSMX_S{section_number:02d}_TLS{number:03d}"
        for number, region in enumerate(valid_regions["region"], start=1)
    }

    one["tls_id"] = one["region"].map(region_to_tls)
    one["inside_tls"] = one["tls_id"].isin(tls_ids)
    one["signed_distance_um"] = np.nan

    for region, tls_id in region_to_tls.items():
        if tls_id not in tls_id_set:
            continue

        region_mask = labels == region
        raster_rows, raster_columns = np.where(region_mask)
        row0 = max(0, raster_rows.min() - 2)
        row1 = min(labels.shape[0], raster_rows.max() + 3)
        column0 = max(0, raster_columns.min() - 2)
        column1 = min(labels.shape[1], raster_columns.max() + 3)

        mask = labels[row0:row1, column0:column1] == region
        distance_inside = distance_transform_edt(mask)
        distance_outside = distance_transform_edt(~mask)
        signed_raster = np.where(
            mask,
            distance_inside - 0.5,
            -(distance_outside - 0.5),
        ) * bin_um

        cell_index = one.index[one["region"].eq(region)]
        region_cells = one.loc[cell_index]
        raster_coordinates = np.vstack(
            [
                (region_cells["y_um"].to_numpy() - y0) / bin_um - 0.5 - row0,
                (region_cells["x_um"].to_numpy() - x0) / bin_um - 0.5 - column0,
            ]
        )
        distance = map_coordinates(
            signed_raster,
            raster_coordinates,
            order=1,
            mode="nearest",
        )
        one.loc[cell_index, "signed_distance_um"] = np.maximum(distance, 0)

    assignment_tables.append(
        one.loc[
            one["inside_tls"],
            [
                "cell_barcode", "sample", "fov", "tls_id", "inside_tls",
                "signed_distance_um", "updated_celltype", "nCount_RNA",
                "CXCR2", "CCR1",
            ],
        ]
    )

    regions["sample"] = sample
    regions["tls_id"] = regions["region"].map(region_to_tls)
    region_tables.append(regions.loc[regions["tls_id"].isin(tls_ids)])

assignments = pd.concat(assignment_tables, ignore_index=True)
tls_summary = pd.concat(region_tables, ignore_index=True)

assignments.to_parquet(
    result_dir / "CosMx_TLS_cell_assignments_152.parquet", index=False
)
tls_summary.to_csv(result_dir / "CosMx_TLS_regions_152.csv", index=False)
