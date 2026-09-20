from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from shapely.geometry import LineString, MultiPoint
from shapely.ops import polygonize, unary_union
from sklearn.cluster import DBSCAN

try:
    from shapely import contains_xy
except ImportError:
    from shapely.vectorized import contains as contains_xy


EPS = 100.0
MIN_SAMPLES = 50
ALPHA_RADIUS = 150.0
MIN_BT_PCT = 10.0
MIN_POLYGON_CELLS = 100
DEFAULT_MINIMUM_TLS_CELLS = 300


def alpha_shape(points: np.ndarray, radius: float):
    if len(points) < 4:
        return MultiPoint(points).convex_hull, "convex_hull"
    triangles = Delaunay(points)
    edges = set()
    for simplex in triangles.simplices:
        a, b, c = points[simplex]
        ab = np.linalg.norm(a - b)
        bc = np.linalg.norm(b - c)
        ca = np.linalg.norm(c - a)
        area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2
        if area == 0:
            continue
        circumradius = ab * bc * ca / (4 * area)
        if circumradius <= radius:
            for i, j in ((simplex[0], simplex[1]), (simplex[1], simplex[2]), (simplex[2], simplex[0])):
                edges.add(tuple(sorted((int(i), int(j)))))
    if not edges:
        return MultiPoint(points).convex_hull, "convex_hull"
    lines = [LineString([points[i], points[j]]) for i, j in edges]
    shape = unary_union(list(polygonize(lines)))
    if shape.is_empty:
        return MultiPoint(points).convex_hull, "convex_hull"
    if shape.geom_type == "MultiPolygon":
        shape = max(shape.geoms, key=lambda polygon: polygon.area)
    return shape, "alpha_shape"


def probability_tnc(b_pct: float, t_pct: float, neutrophil_pct: float, model: dict) -> float:
    coefficient = model["coefficients"]
    eta = model["intercept"]
    eta += coefficient["B_pct"] * math.log1p(b_pct)
    eta += coefficient["T_pct"] * math.log1p(t_pct)
    eta += coefficient["Neutrophil_pct"] * math.log1p(neutrophil_pct)
    return 1 / (1 + math.exp(-eta))


def process_slide(one: pd.DataFrame, model: dict, minimum_tls_cells: int = DEFAULT_MINIMUM_TLS_CELLS):
    immune = one[one["predicted_cell_type"].isin(["B_cell", "T_cell", "Neutrophil", "Other_immune"])].copy()
    if len(immune) < MIN_SAMPLES:
        return [], [], []
    coordinates = immune[["centroid-1", "centroid-0"]].to_numpy(float)
    immune["dbscan_label"] = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES, algorithm="kd_tree").fit_predict(coordinates)
    candidates = []
    polygons = []
    memberships = []
    slide = str(one["fov"].iloc[0])
    all_x = one["centroid-1"].to_numpy(float)
    all_y = one["centroid-0"].to_numpy(float)

    tls_number = 0
    for label, cluster in immune[immune["dbscan_label"].ge(0)].groupby("dbscan_label", sort=True):
        counts = cluster["predicted_cell_type"].value_counts()
        immune_n = len(cluster)
        bt_pct = 100 * (counts.get("B_cell", 0) + counts.get("T_cell", 0)) / immune_n
        if immune_n <= MIN_SAMPLES or bt_pct < MIN_BT_PCT:
            continue
        tls_number += 1
        tls_id = f"{slide}_TLS{tls_number}"
        polygon, boundary_method = alpha_shape(
            cluster[["centroid-1", "centroid-0"]].to_numpy(float), ALPHA_RADIUS
        )
        inside = contains_xy(polygon, all_x, all_y)
        members = one.loc[inside].copy()
        member_counts = members["predicted_cell_type"].value_counts()
        denominator = len(members)
        b_n = int(member_counts.get("B_cell", 0))
        t_n = int(member_counts.get("T_cell", 0))
        neu_n = int(member_counts.get("Neutrophil", 0))
        pass_operational = denominator >= MIN_POLYGON_CELLS
        if denominator:
            b_pct, t_pct, neu_pct = [100 * value / denominator for value in (b_n, t_n, neu_n)]
        else:
            b_pct = t_pct = neu_pct = math.nan
        probability = probability_tnc(b_pct, t_pct, neu_pct, model) if pass_operational else math.nan
        subtype = (
            "TNC-TLS" if probability >= model["probability_cutoff"] else "BF-TLS"
        ) if pass_operational else "not_classified"
        candidates.append({
            "slide": slide,
            "tls_id": tls_id,
            "tls_number": tls_number,
            "dbscan_label": int(label),
            "dbscan_immune_n": immune_n,
            "dbscan_B_n": int(counts.get("B_cell", 0)),
            "dbscan_T_n": int(counts.get("T_cell", 0)),
            "dbscan_Neu_n": int(counts.get("Neutrophil", 0)),
            "dbscan_Other_immune_n": int(counts.get("Other_immune", 0)),
            "dbscan_BT_pct": bt_pct,
            "boundary_method": boundary_method,
            "tls_area_pixel2": polygon.area,
            "denominator_n": denominator,
            "B_n": b_n,
            "T_n": t_n,
            "Neutrophil_n": neu_n,
            "B_pct": b_pct,
            "T_pct": t_pct,
            "Neutrophil_pct": neu_pct,
            "probability_TNC": probability,
            "probability_cutoff": model["probability_cutoff"],
            "predicted_tls_subtype": subtype,
            "pass_tls_filter": pass_operational,
            "minimum_tls_cells": minimum_tls_cells,
            "retained_for_output": pass_operational and denominator >= minimum_tls_cells,
        })
        polygons.append({"slide": slide, "tls_id": tls_id, "wkt": polygon.wkt})
        if denominator:
            members = members[["fov", "label", "centroid-0", "centroid-1", "predicted_cell_type"]].copy()
            members["tls_id"] = tls_id
            memberships.append(members)
    return candidates, polygons, memberships


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-types", type=Path, required=True)
    parser.add_argument("--model", type=Path, default=Path(__file__).parent / "model" / "tlsprofiler_log1p_BTN.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--minimum-tls-cells",
        type=int,
        default=DEFAULT_MINIMUM_TLS_CELLS,
        help="Minimum number of cells within a reconstructed TLS polygon retained in the final output.",
    )
    args = parser.parse_args()
    if args.minimum_tls_cells < MIN_POLYGON_CELLS:
        raise ValueError(f"--minimum-tls-cells must be at least {MIN_POLYGON_CELLS}.")
    model = json.loads(args.model.read_text(encoding="utf-8"))
    cells = pd.read_csv(args.cell_types, low_memory=False)
    required = ["fov", "label", "centroid-0", "centroid-1", "predicted_cell_type"]
    missing = sorted(set(required) - set(cells.columns))
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    candidates = []
    polygons = []
    memberships = []
    for _, one in cells.groupby("fov", sort=False):
        one_candidates, one_polygons, one_memberships = process_slide(
            one, model, minimum_tls_cells=args.minimum_tls_cells
        )
        candidates.extend(one_candidates)
        polygons.extend(one_polygons)
        memberships.extend(one_memberships)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    candidate_table = pd.DataFrame(candidates)
    if candidate_table.empty:
        candidate_table = pd.DataFrame(columns=[
            "slide", "tls_id", "tls_number", "dbscan_label", "dbscan_immune_n",
            "denominator_n", "B_n", "T_n", "Neutrophil_n", "B_pct", "T_pct",
            "Neutrophil_pct", "probability_TNC", "predicted_tls_subtype",
            "pass_tls_filter", "minimum_tls_cells", "retained_for_output",
        ])
    candidate_table.to_csv(args.output_dir / "tls_candidates_and_subtypes.csv", index=False)
    retained = candidate_table[
        candidate_table["retained_for_output"].astype(bool)
    ].copy()
    retained["display_tls_number"] = retained.groupby("slide", sort=False).cumcount() + 1
    retained["display_tls_id"] = (
        retained["slide"].astype(str) + "_TLS" + retained["display_tls_number"].astype(str)
    )
    leading_columns = ["display_tls_number", "display_tls_id"]
    retained = retained[leading_columns + [column for column in retained.columns if column not in leading_columns]]
    retained.to_csv(args.output_dir / "retained_tls_and_subtypes.csv", index=False)
    pd.DataFrame(polygons).to_csv(args.output_dir / "tls_polygons.csv", index=False)
    if memberships:
        pd.concat(memberships, ignore_index=True).to_csv(
            args.output_dir / "cell_tls_membership.csv.gz", index=False, compression="gzip"
        )


if __name__ == "__main__":
    main()
