from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage
from skimage.measure import regionprops_table


DEFAULT_DEEPCELL_IMAGE = (
    "vanvalenlab/deepcell-applications@"
    "sha256:fbfba01a02827ef785cda974e56407d331a35daa2038601c15a4ea680ef951ce"
)


def read_channel(path: Path) -> np.ndarray:
    image = np.asarray(tifffile.imread(path), dtype=np.float32).squeeze()
    if image.ndim != 2:
        raise ValueError(f"Expected a two-dimensional TIFF: {path}")
    return image


def scale_for_segmentation(image: np.ndarray) -> np.ndarray:
    low, high = np.percentile(image[np.isfinite(image)], [1, 99.8])
    if high <= low:
        return np.zeros_like(image, dtype=np.float32)
    return np.clip((image - low) / (high - low), 0, 1).astype(np.float32)


def segment_cells_python(dapi: np.ndarray, membrane: np.ndarray, image_mpp: float) -> np.ndarray:
    try:
        from deepcell.applications import Mesmer
    except ImportError as error:
        raise RuntimeError("DeepCell is not installed. Supply --mask or install requirements.txt") from error
    model_input = np.stack([scale_for_segmentation(dapi), scale_for_segmentation(membrane)], axis=-1)[None, ...]
    prediction = Mesmer().predict(model_input, image_mpp=image_mpp, compartment="whole-cell")
    return np.asarray(prediction[0, ..., 0], dtype=np.int32)


def segment_cells_docker(
    dapi: np.ndarray,
    membrane: np.ndarray,
    image_mpp: float,
    work_dir: Path,
    docker_image: str,
) -> np.ndarray:
    """Run the pinned DeepCell Mesmer CLI without writing the access token to disk."""
    if not os.environ.get("DEEPCELL_ACCESS_TOKEN", "").strip():
        raise RuntimeError(
            "DEEPCELL_ACCESS_TOKEN is not set. Set it in the PowerShell session before running TLSProfiler."
        )
    docker = shutil.which("docker")
    if docker is None:
        raise RuntimeError("Docker was not found on PATH.")

    work_dir.mkdir(parents=True, exist_ok=True)
    input_path = work_dir / "deepcell_input.tiff"
    output_name = "whole_cell_mask.tiff"
    output_path = work_dir / output_name
    if output_path.exists():
        output_path.unlink()

    # The formal project workflow uses the raw DAPI image and the summed raw
    # membrane image in channel-first TIFF format.
    model_input = np.stack([dapi, membrane], axis=0).astype(np.float32, copy=False)
    tifffile.imwrite(input_path, model_input)

    command = [
        docker,
        "run",
        "--rm",
        "-i",
        "-e",
        "DEEPCELL_ACCESS_TOKEN",
        "-v",
        f"{work_dir.resolve()}:/data",
        docker_image,
        "mesmer",
        "--nuclear-image",
        "/data/deepcell_input.tiff",
        "--nuclear-channel",
        "0",
        "--membrane-image",
        "/data/deepcell_input.tiff",
        "--membrane-channel",
        "1",
        "--output-directory",
        "/data",
        "--output-name",
        output_name,
        "--compartment",
        "whole-cell",
        "--image-mpp",
        str(image_mpp),
        "--batch-size",
        "1",
    ]
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0 or not output_path.exists():
        raise RuntimeError(
            "DeepCell Mesmer failed. Check Docker is running, the image is available, "
            "and DEEPCELL_ACCESS_TOKEN is valid."
        )
    return np.asarray(tifffile.imread(output_path)).squeeze().astype(np.int32)


def quantify(mask: np.ndarray, channels: dict[str, np.ndarray], sample_id: str, cofactor: float) -> pd.DataFrame:
    if mask.ndim != 2:
        raise ValueError("The segmentation mask must be two-dimensional")
    labels = np.unique(mask)
    labels = labels[labels > 0]
    if len(labels) == 0:
        raise ValueError("The segmentation mask contains no cells")
    properties = pd.DataFrame(regionprops_table(mask, properties=("label", "area", "centroid")))
    properties = properties.rename(columns={
        "area": "cell_size",
        "centroid-0": "centroid-0",
        "centroid-1": "centroid-1",
    })
    properties.insert(0, "fov", sample_id)
    for marker, image in channels.items():
        mean_intensity = ndimage.mean(image, labels=mask, index=properties["label"].to_numpy(int))
        properties[marker] = np.arcsinh(np.asarray(mean_intensity, float) / cofactor)
    return properties


def build_cell_table(
    sample_id: str,
    channel_paths: dict[str, Path],
    output_dir: Path,
    mask_path: Path | None = None,
    cofactor: float = 0.01,
    image_mpp: float = 0.5,
    segmentation_mode: str = "docker",
    docker_image: str = DEFAULT_DEEPCELL_IMAGE,
) -> tuple[Path, Path]:
    channels = {marker: read_channel(path) for marker, path in channel_paths.items()}
    shapes = {image.shape for image in channels.values()}
    if len(shapes) != 1:
        raise ValueError(f"Channel dimensions differ: {sorted(shapes)}")
    if mask_path is None:
        membrane = np.sum(
            [channels["CD45"], channels["CD20"], channels["CD3"], channels["CD66B"]], axis=0
        )
        if segmentation_mode == "docker":
            mask = segment_cells_docker(
                channels["DAPI"], membrane, image_mpp, output_dir / "deepcell", docker_image
            )
        elif segmentation_mode == "python":
            mask = segment_cells_python(channels["DAPI"], membrane, image_mpp=image_mpp)
        else:
            raise ValueError(f"Unsupported segmentation mode: {segmentation_mode}")
    else:
        mask = np.asarray(tifffile.imread(mask_path)).squeeze().astype(np.int32)
    if mask.shape != next(iter(shapes)):
        raise ValueError("The segmentation mask and marker images have different dimensions")
    output_dir.mkdir(parents=True, exist_ok=True)
    mask_output = output_dir / "whole_cell_mask.tif"
    tifffile.imwrite(mask_output, mask, compression="zlib")
    cell_table = quantify(mask, {key: value for key, value in channels.items() if key != "DAPI"}, sample_id, cofactor)
    table_output = output_dir / "cell_table_arcsinh.csv.gz"
    cell_table.to_csv(table_output, index=False, compression="gzip")
    return table_output, mask_output
