from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from segment_and_quantify import DEFAULT_DEEPCELL_IMAGE, build_cell_table
from visualize_results import generate_visualizations


ROOT = Path(__file__).resolve().parent


def run_pipeline(
    sample_id: str,
    channel_paths: dict[str, Path],
    output_dir: Path,
    mask_path: Path | None = None,
    cofactor: float = 0.01,
    image_mpp: float = 0.5,
    minimum_tls_cells: int = 300,
    segmentation_mode: str = "docker",
    docker_image: str = DEFAULT_DEEPCELL_IMAGE,
) -> dict:
    segmentation_dir = output_dir / "01-segmentation"
    phenotype_dir = output_dir / "02-cell_types"
    tls_dir = output_dir / "03-tls"
    visualization_dir = output_dir / "04-visualizations"
    cell_table, mask_output = build_cell_table(
        sample_id,
        channel_paths,
        segmentation_dir,
        mask_path,
        cofactor,
        image_mpp,
        segmentation_mode,
        docker_image,
    )
    subprocess.run([
        sys.executable, str(ROOT / "call_mihc_cell_types.py"),
        "--input", str(cell_table), "--output-dir", str(phenotype_dir)
    ], check=True)
    subprocess.run([
        sys.executable, str(ROOT / "detect_and_classify_tls.py"),
        "--cell-types", str(phenotype_dir / "cell_type_assignments.csv.gz"),
        "--model", str(ROOT / "model" / "tlsprofiler_log1p_BTN.json"),
        "--output-dir", str(tls_dir),
        "--minimum-tls-cells", str(minimum_tls_cells),
    ], check=True)
    candidate_tls_path = tls_dir / "tls_candidates_and_subtypes.csv"
    tls_path = tls_dir / "retained_tls_and_subtypes.csv"
    candidates = (
        pd.read_csv(candidate_tls_path)
        if candidate_tls_path.exists() and candidate_tls_path.stat().st_size
        else pd.DataFrame()
    )
    retained = pd.read_csv(tls_path) if tls_path.exists() and tls_path.stat().st_size else pd.DataFrame()
    visualizations = generate_visualizations(
        channel_paths=channel_paths,
        mask_path=mask_output,
        assignments_path=phenotype_dir / "cell_type_assignments.csv.gz",
        thresholds_path=phenotype_dir / "marker_thresholds.csv",
        tls_table_path=tls_path,
        polygons_path=tls_dir / "tls_polygons.csv",
        output_dir=visualization_dir,
    )
    retained_subtypes = retained["predicted_tls_subtype"].dropna().unique().tolist() if not retained.empty else []
    if not retained_subtypes:
        sample_subtype = "non-TLS"
    elif len(retained_subtypes) == 1:
        sample_subtype = retained_subtypes[0]
    else:
        sample_subtype = "mixed"
    summary = {
        "sample_id": sample_id,
        "minimum_tls_cells": int(minimum_tls_cells),
        "candidate_tls_count": int(len(candidates)),
        "retained_tls_count": int(len(retained)),
        "BF_TLS_count": int((retained.get("predicted_tls_subtype") == "BF-TLS").sum()) if not retained.empty else 0,
        "TNC_TLS_count": int((retained.get("predicted_tls_subtype") == "TNC-TLS").sum()) if not retained.empty else 0,
        "sample_subtype": sample_subtype,
        "mask": str(mask_output),
        "cell_table": str(cell_table),
        "tls_table": str(tls_path),
        "candidate_tls_table": str(candidate_tls_path),
        "visualizations": {name: str(path) for name, path in visualizations.items()},
    }
    (output_dir / "TLSProfiler_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--dapi", type=Path, required=True)
    parser.add_argument("--cd45", type=Path, required=True)
    parser.add_argument("--cd20", type=Path, required=True)
    parser.add_argument("--cd3", type=Path, required=True)
    parser.add_argument("--cd66b", type=Path, required=True)
    parser.add_argument("--mask", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--cofactor",
        type=float,
        default=0.01,
        help="ARK-compatible arcsinh cofactor applied after per-cell mean intensity (default: 0.01).",
    )
    parser.add_argument("--image-mpp", type=float, default=0.5)
    parser.add_argument(
        "--minimum-tls-cells",
        type=int,
        default=300,
        help="Minimum polygon cell count retained in the final TLS results (default: 300).",
    )
    parser.add_argument("--segmentation-mode", choices=["docker", "python"], default="docker")
    parser.add_argument(
        "--deepcell-image",
        default=DEFAULT_DEEPCELL_IMAGE,
        help="DeepCell Docker image used when --segmentation-mode docker is selected.",
    )
    args = parser.parse_args()
    summary = run_pipeline(
        sample_id=args.sample_id,
        channel_paths={
            "DAPI": args.dapi,
            "CD45": args.cd45,
            "CD20": args.cd20,
            "CD3": args.cd3,
            "CD66B": args.cd66b,
        },
        output_dir=args.output_dir,
        mask_path=args.mask,
        cofactor=args.cofactor,
        image_mpp=args.image_mpp,
        minimum_tls_cells=args.minimum_tls_cells,
        segmentation_mode=args.segmentation_mode,
        docker_image=args.deepcell_image,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
