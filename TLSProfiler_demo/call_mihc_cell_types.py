from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from skimage.filters import threshold_otsu


MARKERS = ["CD45", "CD20", "CD3", "CD66B"]
LINEAGE_MARKERS = ["CD20", "CD3", "CD66B"]
CELL_SIZE_MIN = 150.0
CELL_SIZE_MAX = 1500.0
TRIM_LOW = 1.0
TRIM_HIGH = 99.8
MINIMUM_OTSU_CELLS = 100
MODE_OFFSET = 0.5
CD66B_DOMINANCE = 1.5


def trimmed_otsu(values: np.ndarray) -> tuple[float, int, int, str]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    n_original = len(values)
    if n_original == 0:
        return math.nan, 0, 0, "no_finite_values"
    trimmed = values
    if n_original > MINIMUM_OTSU_CELLS:
        low, high = np.percentile(values, [TRIM_LOW, TRIM_HIGH])
        candidate = values[(values >= low) & (values <= high)]
        if len(candidate) >= 2:
            trimmed = candidate
    if np.unique(trimmed).size < 2:
        return float(np.median(trimmed)), n_original, len(trimmed), "median_constant_fallback"
    return float(threshold_otsu(trimmed)), n_original, len(trimmed), "trimmed_otsu"


def first_major_mode(values: np.ndarray) -> tuple[float, float, int, str]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) < MINIMUM_OTSU_CELLS or np.unique(values).size < 3:
        return math.nan, math.nan, len(values), "insufficient_positive_values"
    low, high = np.percentile(values, [TRIM_LOW, TRIM_HIGH])
    values = values[(values >= low) & (values <= high)]
    if len(values) < MINIMUM_OTSU_CELLS or high <= low:
        return math.nan, math.nan, len(values), "insufficient_trimmed_positive_values"
    edges = np.linspace(low, high, 257)
    counts, _ = np.histogram(values, bins=edges)
    density = gaussian_filter1d(counts.astype(float), sigma=3.0, mode="nearest")
    centers = (edges[:-1] + edges[1:]) / 2
    peaks, properties = find_peaks(
        density,
        prominence=max(float(density.max()) * 0.02, 1e-12),
        distance=8,
    )
    qualifying = [
        int(peak)
        for i, peak in enumerate(peaks)
        if density[peak] >= density.max() * 0.25
        and properties["prominences"][i] >= density.max() * 0.10
    ]
    peak = min(qualifying, key=lambda i: centers[i]) if qualifying else int(np.argmax(density))
    mode = float(centers[peak])
    floor = min(mode + MODE_OFFSET, float(np.percentile(values, 95)))
    return mode, floor, len(values), "smoothed_histogram_first_major_mode_plus_offset"


def estimate_thresholds(cells: pd.DataFrame, slide: str):
    cd45, n_values, n_trimmed, method = trimmed_otsu(cells["CD45"].to_numpy())
    lineage_pool = cells.loc[cells["CD45"] > cd45]
    lineage_fallback = len(lineage_pool) < MINIMUM_OTSU_CELLS
    if lineage_fallback:
        lineage_pool = cells
    thresholds = {"CD45": cd45}
    rows = [{
        "fov": slide,
        "marker": "CD45",
        "threshold": cd45,
        "positive_pool_otsu": cd45,
        "first_major_positive_mode": np.nan,
        "first_major_mode_plus_offset_floor": np.nan,
        "threshold_method": method,
        "threshold_pool": "all_size_filtered_cells",
        "n_values_for_threshold": n_values,
        "n_values_after_trim": n_trimmed,
        "n_positive_values_for_lineage_threshold": np.nan,
        "lineage_pool_fallback_to_all_cells": False,
    }]
    for marker in LINEAGE_MARKERS:
        positive = lineage_pool.loc[lineage_pool[marker] > 0, marker].to_numpy(float)
        if len(positive) >= MINIMUM_OTSU_CELLS:
            pool = positive
            pool_name = "CD45_positive_strictly_positive_marker_values"
        else:
            pool = lineage_pool[marker].to_numpy(float)
            pool_name = "CD45_positive_all_values_fallback"
        otsu, n_values, n_trimmed, otsu_method = trimmed_otsu(pool)
        mode, floor, n_positive, mode_method = first_major_mode(positive)
        valid = [value for value in (otsu, floor) if np.isfinite(value)]
        threshold = float(max(valid)) if valid else math.nan
        thresholds[marker] = threshold
        rows.append({
            "fov": slide,
            "marker": marker,
            "threshold": threshold,
            "positive_pool_otsu": otsu,
            "first_major_positive_mode": mode,
            "first_major_mode_plus_offset_floor": floor,
            "threshold_method": f"max({otsu_method}, {mode_method})",
            "threshold_pool": pool_name,
            "n_values_for_threshold": n_values,
            "n_values_after_trim": n_trimmed,
            "n_positive_values_for_lineage_threshold": n_positive,
            "lineage_pool_fallback_to_all_cells": lineage_fallback,
        })
    threshold_table = pd.DataFrame(rows)
    for marker in MARKERS:
        threshold_table.loc[threshold_table["marker"].eq(marker), "target_positive_pct_all_retained_cells"] = (
            cells[marker].gt(thresholds[marker]).mean() * 100
        )
    return thresholds, threshold_table


def assign_cell_types(cells: pd.DataFrame, thresholds: dict[str, float]) -> pd.DataFrame:
    out = cells.copy()
    for marker in MARKERS:
        out[f"positive_{marker}"] = out[marker] > thresholds[marker]
        out[f"fold_{marker}"] = out[marker] / max(abs(thresholds[marker]), 1e-8)
    immune = out["positive_CD45"].to_numpy(bool)
    p20 = out["positive_CD20"].to_numpy(bool)
    p3 = out["positive_CD3"].to_numpy(bool)
    p66 = out["positive_CD66B"].to_numpy(bool)
    f20 = out["fold_CD20"].to_numpy(float)
    f3 = out["fold_CD3"].to_numpy(float)
    f66 = out["fold_CD66B"].to_numpy(float)
    cell_type = np.full(len(out), "Other_cell", dtype=object)
    cell_type[immune] = "Other_immune"
    neutrophil = immune & p66 & ((~p20 & ~p3) | (f66 >= CD66B_DOMINANCE * np.maximum(f20, f3)))
    cell_type[neutrophil] = "Neutrophil"
    remaining = immune & ~neutrophil
    cell_type[remaining & p20 & ~p3] = "B_cell"
    cell_type[remaining & ~p20 & p3] = "T_cell"
    both = remaining & p20 & p3
    cell_type[both & (f20 >= f3)] = "B_cell"
    cell_type[both & (f3 > f20)] = "T_cell"
    out["predicted_is_immune"] = immune
    out["predicted_cell_type"] = cell_type
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    data = pd.read_csv(args.input, low_memory=False)
    required = ["fov", "label", "cell_size", "centroid-0", "centroid-1", *MARKERS]
    missing = sorted(set(required) - set(data.columns))
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    for column in required[1:]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data[
        data["cell_size"].gt(CELL_SIZE_MIN)
        & data["cell_size"].lt(CELL_SIZE_MAX)
    ].dropna(subset=required[1:])

    assignments = []
    threshold_tables = []
    for slide, one in data.groupby("fov", sort=False):
        thresholds, threshold_table = estimate_thresholds(one, str(slide))
        assignments.append(assign_cell_types(one, thresholds))
        threshold_tables.append(threshold_table)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(assignments, ignore_index=True).to_csv(
        args.output_dir / "cell_type_assignments.csv.gz", index=False, compression="gzip"
    )
    pd.concat(threshold_tables, ignore_index=True).to_csv(
        args.output_dir / "marker_thresholds.csv", index=False
    )


if __name__ == "__main__":
    main()
