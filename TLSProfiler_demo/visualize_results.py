from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.patches import Patch
from shapely import wkt
from skimage.draw import polygon as rasterize_polygon
from skimage.segmentation import find_boundaries


CHANNELS = ["DAPI", "CD45", "CD20", "CD3", "CD66B"]
MARKERS = ["CD45", "CD20", "CD3", "CD66B"]
CHANNEL_COLORS = {
    "DAPI": "#4D7CFE",
    "CD45": "#F2F2F2",
    "CD20": "#E5368A",
    "CD3": "#00A6A6",
    "CD66B": "#F2B134",
}
CELL_TYPE_COLORS = {
    "B_cell": "#3B82C4",
    "T_cell": "#F28E45",
    "Neutrophil": "#D94B45",
    "Other_immune": "#85899A",
    "Other_cell": "#D6D7DA",
}
CELL_TYPE_LABELS = {
    "B_cell": "B cell",
    "T_cell": "T cell",
    "Neutrophil": "Neutrophil",
    "Other_immune": "Other immune",
    "Other_cell": "Other cell",
}
TLS_COLORS = {
    "BF-TLS": "#2E7DAD",
    "TNC-TLS": "#C81D25",
}


def read_image(path: Path) -> np.ndarray:
    image = np.asarray(tifffile.imread(path), dtype=np.float32).squeeze()
    if image.ndim != 2:
        raise ValueError(f"Expected a two-dimensional image: {path}")
    return image


def normalize_display(image: np.ndarray, low_pct: float = 1.0, high_pct: float = 99.8) -> np.ndarray:
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        return np.zeros(image.shape, dtype=np.float32)
    low, high = np.percentile(finite, [low_pct, high_pct])
    if high <= low:
        return np.zeros(image.shape, dtype=np.float32)
    return np.clip((image - low) / (high - low), 0, 1).astype(np.float32)


def hex_rgb(value: str) -> np.ndarray:
    value = value.lstrip("#")
    return np.array([int(value[i : i + 2], 16) for i in (0, 2, 4)], dtype=np.float32) / 255


def composite_image(images: dict[str, np.ndarray]) -> np.ndarray:
    rgb = np.zeros((*next(iter(images.values())).shape, 3), dtype=np.float32)
    weights = {"DAPI": 0.70, "CD45": 0.34, "CD20": 0.92, "CD3": 0.88, "CD66B": 0.92}
    for marker in CHANNELS:
        rgb += normalize_display(images[marker])[:, :, None] * hex_rgb(CHANNEL_COLORS[marker]) * weights[marker]
    return np.clip(rgb, 0, 1) ** 0.82


def style_image_axis(axis, title: str) -> None:
    axis.set_title(title, fontsize=11, pad=6)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)


def save_input_channels(channel_paths: dict[str, Path], output: Path) -> None:
    images = {marker: read_image(channel_paths[marker]) for marker in CHANNELS}
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2), constrained_layout=True)
    for axis, marker in zip(axes.flat[:5], CHANNELS):
        shown = normalize_display(images[marker])
        axis.imshow(shown, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        style_image_axis(axis, marker)
    axes.flat[5].imshow(composite_image(images), interpolation="nearest")
    style_image_axis(axes.flat[5], "Five-channel composite")
    fig.suptitle("Input image channels", fontsize=15, fontweight="semibold")
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def save_thresholds(assignments_path: Path, thresholds_path: Path, output: Path) -> None:
    cells = pd.read_csv(assignments_path, low_memory=False)
    thresholds = pd.read_csv(thresholds_path).set_index("marker")
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.8))
    for axis, marker in zip(axes.flat, MARKERS):
        values = pd.to_numeric(cells[marker], errors="coerce").dropna().to_numpy(float)
        threshold = float(thresholds.loc[marker, "threshold"])
        low, high = np.percentile(values, [0.5, 99.5])
        low = min(low, threshold)
        high = max(high, threshold)
        if high <= low:
            high = low + 1
        bins = np.linspace(low, high, 75)
        counts, edges = np.histogram(np.clip(values, low, high), bins=bins, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        widths = np.diff(edges)
        colors = np.where(centers > threshold, "#D95F59", "#B7C3CF")
        axis.bar(centers, counts, width=widths, color=colors, edgecolor="none", alpha=0.92)
        axis.axvline(threshold, color="#A51C30", linewidth=2)
        positive = 100 * float(np.mean(values > threshold))
        axis.text(
            0.98,
            0.95,
            f"Threshold = {threshold:.3f}\nPositive = {positive:.1f}%\nn = {len(values):,}",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#D5D7DA", "boxstyle": "square,pad=0.35"},
        )
        axis.set_title(marker, fontsize=12, fontweight="semibold")
        axis.set_xlabel("Arcsinh intensity")
        axis.set_ylabel("Density")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#ECEDEF", linewidth=0.7)
    fig.suptitle("Marker thresholding", fontsize=15, fontweight="semibold")
    fig.subplots_adjust(left=0.07, right=0.985, top=0.90, bottom=0.11, hspace=0.37, wspace=0.14)
    fig.text(
        0.5,
        0.025,
        "Distributions show retained cells; lineage-marker thresholds are estimated within the predicted CD45-positive pool.",
        ha="center",
        fontsize=9,
        color="#555555",
    )
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def cell_type_rgb(mask: np.ndarray, assignments: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    max_label = int(mask.max())
    lookup = np.zeros((max_label + 1, 3), dtype=np.float32)
    assigned_lookup = np.zeros(max_label + 1, dtype=bool)
    labels = pd.to_numeric(assignments["label"], errors="coerce").fillna(-1).astype(int).to_numpy()
    cell_types = assignments["predicted_cell_type"].astype(str).to_numpy()
    valid = (labels > 0) & (labels <= max_label)
    for label, cell_type in zip(labels[valid], cell_types[valid]):
        lookup[label] = hex_rgb(CELL_TYPE_COLORS.get(cell_type, "#D6D7DA"))
        assigned_lookup[label] = True
    return lookup[np.clip(mask, 0, max_label)], assigned_lookup[np.clip(mask, 0, max_label)]


def save_cell_types(
    dapi_path: Path,
    mask_path: Path,
    assignments_path: Path,
    output: Path,
) -> None:
    dapi = normalize_display(read_image(dapi_path))
    mask = np.asarray(tifffile.imread(mask_path)).squeeze().astype(np.int32)
    assignments = pd.read_csv(assignments_path, low_memory=False)
    type_rgb, assigned = cell_type_rgb(mask, assignments)
    boundaries = find_boundaries(mask, mode="outer")
    map_rgb = np.ones((*mask.shape, 3), dtype=np.float32)
    map_rgb[assigned] = type_rgb[assigned]
    map_rgb[boundaries] = 0.18
    dapi_rgb = np.stack([dapi * 0.20, dapi * 0.34, dapi * 0.78], axis=-1)
    overlay = dapi_rgb.copy()
    overlay[assigned] = 0.30 * dapi_rgb[assigned] + 0.70 * type_rgb[assigned]
    overlay[boundaries] *= 0.55

    order = ["B_cell", "T_cell", "Neutrophil", "Other_immune", "Other_cell"]
    counts = assignments["predicted_cell_type"].value_counts()
    fig = plt.figure(figsize=(16, 6.2), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.72])
    ax_map = fig.add_subplot(grid[0, 0])
    ax_overlay = fig.add_subplot(grid[0, 1])
    ax_bar = fig.add_subplot(grid[0, 2])
    ax_map.imshow(map_rgb, interpolation="nearest")
    style_image_axis(ax_map, "Cell-type assignment")
    ax_overlay.imshow(overlay, interpolation="nearest")
    style_image_axis(ax_overlay, "Cell types over DAPI")

    values = np.array([int(counts.get(name, 0)) for name in order])
    y = np.arange(len(order))
    ax_bar.barh(y, values, color=[CELL_TYPE_COLORS[name] for name in order], height=0.65)
    ax_bar.set_yticks(y, [CELL_TYPE_LABELS[name] for name in order])
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Cells")
    ax_bar.set_title(f"Retained cells (n = {values.sum():,})", fontsize=11)
    ax_bar.spines[["top", "right", "left"]].set_visible(False)
    ax_bar.grid(axis="x", color="#ECEDEF", linewidth=0.7)
    for index, value in enumerate(values):
        ax_bar.text(value, index, f"  {value:,}", va="center", fontsize=9)

    handles = [Patch(color=CELL_TYPE_COLORS[name], label=CELL_TYPE_LABELS[name]) for name in order]
    fig.legend(handles=handles, loc="outside lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("Cell-type identification", fontsize=15, fontweight="semibold")
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def polygon_parts(geometry):
    if geometry.geom_type == "Polygon":
        return [geometry]
    if geometry.geom_type == "MultiPolygon":
        return list(geometry.geoms)
    return []


def label_position(x: float, y: float, occupied: list[tuple[float, float]], shape: tuple[int, int]):
    height, width = shape
    offsets = [
        (0, 0),
        (0.10 * width, 0.07 * height),
        (-0.10 * width, -0.07 * height),
        (0.10 * width, -0.07 * height),
        (-0.10 * width, 0.07 * height),
        (0, 0.12 * height),
        (0, -0.12 * height),
    ]
    minimum_x = max(320.0, 0.08 * width)
    minimum_y = max(260.0, 0.06 * height)
    margin_x = max(220.0, 0.055 * width)
    margin_y = max(170.0, 0.04 * height)
    for x_offset, y_offset in offsets:
        candidate_x = float(np.clip(x + x_offset, margin_x, width - margin_x))
        candidate_y = float(np.clip(y + y_offset, margin_y, height - margin_y))
        collision = any(
            abs(candidate_x - used_x) < minimum_x and abs(candidate_y - used_y) < minimum_y
            for used_x, used_y in occupied
        )
        if not collision:
            occupied.append((candidate_x, candidate_y))
            return candidate_x, candidate_y
    occupied.append((x, y))
    return x, y


def save_tls_regions(
    channel_paths: dict[str, Path],
    tls_table_path: Path,
    polygons_path: Path,
    output: Path,
) -> None:
    images = {marker: read_image(channel_paths[marker]) for marker in CHANNELS}
    composite = composite_image(images)
    dapi = normalize_display(images["DAPI"])
    dapi_rgb = np.stack([dapi * 0.20, dapi * 0.34, dapi * 0.78], axis=-1)
    tls = pd.read_csv(tls_table_path) if tls_table_path.exists() and tls_table_path.stat().st_size else pd.DataFrame()
    try:
        polygons = pd.read_csv(polygons_path)
    except pd.errors.EmptyDataError:
        polygons = pd.DataFrame(columns=["tls_id", "wkt"])
    records = tls.merge(polygons[["tls_id", "wkt"]], on="tls_id", how="inner") if not tls.empty else tls

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.2), constrained_layout=True)
    axes[0].imshow(composite, interpolation="nearest")
    axes[1].imshow(dapi_rgb, interpolation="nearest")
    style_image_axis(axes[0], "TLS boundaries on five-channel composite")
    style_image_axis(axes[1], "TLS region masks on DAPI")
    filled = np.zeros((*dapi.shape, 3), dtype=np.float32)
    filled_alpha = np.zeros(dapi.shape, dtype=np.float32)
    label_records = []

    for _, row in records.iterrows():
        geometry = wkt.loads(row["wkt"])
        subtype = str(row["predicted_tls_subtype"])
        color_name = TLS_COLORS[subtype]
        color = hex_rgb(color_name)
        for part in polygon_parts(geometry):
            xy = np.asarray(part.exterior.coords)
            for axis in axes:
                axis.plot(xy[:, 0], xy[:, 1], color=color_name, linewidth=2.4)
            rr, cc = rasterize_polygon(xy[:, 1], xy[:, 0], shape=dapi.shape)
            filled[rr, cc] = color
            filled_alpha[rr, cc] = 0.72
            for interior in part.interiors:
                hole = np.asarray(interior.coords)
                rr_hole, cc_hole = rasterize_polygon(hole[:, 1], hole[:, 0], shape=dapi.shape)
                filled_alpha[rr_hole, cc_hole] = 0
        center = geometry.representative_point()
        tls_number = int(row["display_tls_number"])
        probability = row.get("probability_TNC", np.nan)
        probability_text = f"pTNC={float(probability):.3f}" if pd.notna(probability) else "pTNC=NA"
        label = f"TLS {tls_number}\n{subtype}\n{probability_text}"
        label_records.append((center.x, center.y, label, color_name))

    occupied = []
    for center_x, center_y, label, color_name in label_records:
        label_x, label_y = label_position(center_x, center_y, occupied, dapi.shape)
        for axis in axes:
            axis.annotate(
                label,
                xy=(center_x, center_y),
                xytext=(label_x, label_y),
                textcoords="data",
                color="white",
                fontsize=7.5,
                ha="center",
                va="center",
                bbox={"facecolor": color_name, "edgecolor": "white", "alpha": 0.86, "boxstyle": "square,pad=0.25"},
                arrowprops={"arrowstyle": "-", "color": color_name, "linewidth": 1.0}
                if (label_x, label_y) != (center_x, center_y)
                else None,
            )

    mask_overlay = dapi_rgb * (1 - filled_alpha[:, :, None]) + filled * filled_alpha[:, :, None]
    axes[1].images[0].set_data(mask_overlay)
    if records.empty:
        for axis in axes:
            axis.text(
                0.5,
                0.5,
                "No TLS met the minimum cell threshold",
                transform=axis.transAxes,
                ha="center",
                va="center",
            )

    handles = [
        Patch(color=TLS_COLORS["BF-TLS"], label="BF-TLS"),
        Patch(color=TLS_COLORS["TNC-TLS"], label="TNC-TLS"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Retained TLS regions and subtype inference", fontsize=15, fontweight="semibold")
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def generate_visualizations(
    channel_paths: dict[str, Path],
    mask_path: Path,
    assignments_path: Path,
    thresholds_path: Path,
    tls_table_path: Path,
    polygons_path: Path,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "input_channels": output_dir / "01-input_channels.png",
        "marker_thresholds": output_dir / "02-marker_thresholds.png",
        "cell_types": output_dir / "03-cell_types.png",
        "tls_regions": output_dir / "04-tls_regions.png",
    }
    save_input_channels(channel_paths, outputs["input_channels"])
    save_thresholds(assignments_path, thresholds_path, outputs["marker_thresholds"])
    save_cell_types(channel_paths["DAPI"], mask_path, assignments_path, outputs["cell_types"])
    save_tls_regions(channel_paths, tls_table_path, polygons_path, outputs["tls_regions"])
    return outputs
